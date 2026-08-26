/*
 * CAGRA-style graph search over node-level PQ windows with shared-memory ADC.
 *
 * Graph traversal and the open-addressed visited hashmap follow NVIDIA cuVS
 * CAGRA (Apache-2.0).  Distance evaluation is replaced with a per-query
 * lookup table in dynamic shared memory because cuVS does not accept a custom
 * device distance.
 *
 * SPDX-FileCopyrightText: Copyright (c) 2023-2026, NVIDIA CORPORATION & AFFILIATES.
 * SPDX-FileCopyrightText: Copyright (c) 2026, GINflow authors.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "cagra_pq_adc.hpp"

#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include <math_constants.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ginflow::pq_cagra {
namespace {

constexpr std::uint32_t kCheckedBit = 0x80000000u;
constexpr std::uint32_t kIndexMask = 0x7fffffffu;
constexpr std::uint32_t kEmpty = 0xffffffffu;
constexpr std::uint32_t kMaxItopk = 8192;
constexpr std::uint32_t kMaxGraphDegree = 128;
constexpr std::uint32_t kMaxIntermediate = 128;
constexpr std::uint32_t kMaxSearchWidth = 32;
constexpr std::uint32_t kMaxSmemBytes = 96u * 1024u;
constexpr std::uint32_t kBruteForceLimit = 16384;

#define PQ_CAGRA_CHECK(expr)                                                         \
    do {                                                                             \
        const cudaError_t _err = (expr);                                             \
        if (_err != cudaSuccess) {                                                   \
            throw std::runtime_error(std::string("CUDA error: ") +                   \
                                     cudaGetErrorString(_err) + " (" + #expr + ")"); \
        }                                                                            \
    } while (0)

void check_launch() {
    PQ_CAGRA_CHECK(cudaGetLastError());
    PQ_CAGRA_CHECK(cudaDeviceSynchronize());
}

void report_progress(const char* phase, std::uint32_t percent, const char* detail = nullptr) {
    if (detail && detail[0] != '\0') {
        std::fprintf(stdout,
                     "PQ_CAGRA_PROGRESS scope=graph_build phase=%s percent=%u %s\n",
                     phase,
                     percent,
                     detail);
    } else {
        std::fprintf(stdout,
                     "PQ_CAGRA_PROGRESS scope=graph_build phase=%s percent=%u\n",
                     phase,
                     percent);
    }
    std::fflush(stdout);
}

void report_iteration_progress(const char* phase,
                               std::uint32_t percent,
                               std::uint32_t iteration,
                               std::uint32_t total) {
    char detail[48];
    std::snprintf(detail, sizeof(detail), "iteration=%u/%u", iteration, total);
    report_progress(phase, percent, detail);
}

std::uint32_t round_up_pow2(std::uint32_t value) {
    std::uint32_t padded = 1;
    while (padded < value) {
        if (padded > (1u << 20)) throw std::runtime_error("bitonic padded length overflow");
        padded <<= 1;
    }
    return padded;
}

// Adapted from cuVS CAGRA hashmap (linear probing, empty key = ~0).
namespace hashmap {

__device__ __forceinline__ std::uint32_t size(std::uint32_t bitlen) { return 1u << bitlen; }

__device__ __forceinline__ std::uint32_t insert(std::uint32_t* table,
                                                std::uint32_t bitlen,
                                                std::uint32_t key) {
    const std::uint32_t n = size(bitlen);
    const std::uint32_t mask = n - 1u;
    std::uint32_t index = (key ^ (key >> bitlen)) & mask;
    for (std::uint32_t probe = 0; probe < n; ++probe) {
        const std::uint32_t old = atomicCAS(table + index, kEmpty, key);
        if (old == kEmpty) return 1;
        if (old == key) return 0;
        index = (index + 1u) & mask;
    }
    return 0;
}

}  // namespace hashmap

__device__ __forceinline__ float sdc_score(const std::uint8_t* a,
                                           const std::uint8_t* b,
                                           const float* similarity,
                                           std::uint32_t code_dim,
                                           std::uint32_t pq_m,
                                           std::uint32_t ksub) {
    float score = 0.0f;
    const std::uint32_t table_stride = ksub * ksub;
    const std::uint32_t window_size = code_dim / pq_m;
    for (std::uint32_t position = 0; position < window_size; ++position) {
        const std::uint8_t* a_node = a + position * pq_m;
        const std::uint8_t* b_node = b + position * pq_m;
        for (std::uint32_t sub = 0; sub < pq_m; ++sub) {
            const std::uint32_t ca = a_node[sub];
            const std::uint32_t cb = b_node[sub];
            score += similarity[sub * table_stride + ca * ksub + cb];
        }
    }
    return score;
}

__device__ __forceinline__ float warp_adc(const __half* lut,
                                          const std::uint8_t* codes,
                                          std::uint32_t code_dim,
                                          std::uint32_t ksub) {
    float sum = 0.0f;
    const unsigned lane = threadIdx.x & 31u;
    for (std::uint32_t i = lane; i < code_dim; i += 32u) {
        const std::uint32_t code = codes[i];
        sum += __half2float(lut[i * ksub + code]);
    }
#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1) {
        sum += __shfl_xor_sync(0xffffffffu, sum, offset);
    }
    return sum;
}

__device__ void bitonic_sort_kv(float* keys, std::uint32_t* vals, std::uint32_t n, std::uint32_t padded) {
    for (std::uint32_t i = threadIdx.x; i < padded; i += blockDim.x) {
        if (i >= n) {
            keys[i] = CUDART_INF_F;
            vals[i] = kEmpty;
        }
    }
    __syncthreads();
    for (std::uint32_t size = 2; size <= padded; size <<= 1) {
        for (std::uint32_t stride = size >> 1; stride > 0; stride >>= 1) {
            for (std::uint32_t i = threadIdx.x; i < padded; i += blockDim.x) {
                const std::uint32_t j = i ^ stride;
                if (j > i) {
                    const bool ascending = (i & size) == 0;
                    const bool should_swap =
                        ascending ? (keys[i] > keys[j]) : (keys[i] < keys[j]);
                    if (should_swap) {
                        const float tk = keys[i];
                        keys[i] = keys[j];
                        keys[j] = tk;
                        const std::uint32_t tv = vals[i];
                        vals[i] = vals[j];
                        vals[j] = tv;
                    }
                }
            }
            __syncthreads();
        }
    }
}

__global__ void knn_sdc_kernel(const std::uint8_t* __restrict__ codes,
                               const float* __restrict__ similarity,
                               std::uint32_t* __restrict__ knn_idx,
                               float* __restrict__ knn_sim,
                               std::uint32_t n,
                               std::uint32_t code_dim,
                               std::uint32_t pq_m,
                               std::uint32_t ksub,
                               std::uint32_t k) {
    const std::uint32_t row = blockIdx.x;
    if (row >= n) return;

    extern __shared__ char smem[];
    const std::uint32_t nwarps = blockDim.x / 32u;
    const std::uint32_t warp = threadIdx.x / 32u;
    const std::uint32_t lane = threadIdx.x & 31u;
    float* all_sim = reinterpret_cast<float*>(smem);
    std::uint32_t* all_idx = reinterpret_cast<std::uint32_t*>(all_sim + nwarps * k);
    int* all_fill = reinterpret_cast<int*>(all_idx + nwarps * k);
    float* warp_sim = all_sim + warp * k;
    std::uint32_t* warp_idx = all_idx + warp * k;
    int* warp_fill = all_fill + warp;

    if (lane == 0) warp_fill[0] = 0;
    __syncwarp();

    const std::uint8_t* row_codes = codes + static_cast<std::size_t>(row) * code_dim;
    const std::uint32_t stride = nwarps * 32u;
    for (std::uint32_t base = 0; base < n; base += stride) {
        const std::uint32_t col = base + warp * 32u + lane;
        float score = -CUDART_INF_F;
        std::uint32_t other = kEmpty;
        if (col < n && col != row) {
            score = sdc_score(row_codes,
                              codes + static_cast<std::size_t>(col) * code_dim,
                              similarity,
                              code_dim,
                              pq_m,
                              ksub);
            other = col;
        }
        float lane_scores[32];
        std::uint32_t lane_ids[32];
        for (int src = 0; src < 32; ++src) {
            lane_scores[src] = __shfl_sync(0xffffffffu, score, src);
            lane_ids[src] = __shfl_sync(0xffffffffu, other, src);
        }
        if (lane == 0) {
            for (int src = 0; src < 32; ++src) {
                if (lane_ids[src] == kEmpty) continue;
                const float s = lane_scores[src];
                int fill = warp_fill[0];
                if (fill < static_cast<int>(k)) {
                    warp_sim[fill] = s;
                    warp_idx[fill] = lane_ids[src];
                    warp_fill[0] = fill + 1;
                } else {
                    int worst = 0;
                    for (int t = 1; t < static_cast<int>(k); ++t) {
                        if (warp_sim[t] < warp_sim[worst]) worst = t;
                    }
                    if (s > warp_sim[worst]) {
                        warp_sim[worst] = s;
                        warp_idx[worst] = lane_ids[src];
                    }
                }
            }
        }
        __syncwarp();
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        float best_sim[kMaxIntermediate];
        std::uint32_t best_idx[kMaxIntermediate];
        int fill = 0;
        for (std::uint32_t w = 0; w < nwarps; ++w) {
            float* src_sim = all_sim + w * k;
            std::uint32_t* src_idx = all_idx + w * k;
            const int src_fill = all_fill[w];
            for (int t = 0; t < src_fill; ++t) {
                const float s = src_sim[t];
                const std::uint32_t id = src_idx[t];
                if (fill < static_cast<int>(k)) {
                    best_sim[fill] = s;
                    best_idx[fill] = id;
                    ++fill;
                } else {
                    int worst = 0;
                    for (int u = 1; u < static_cast<int>(k); ++u) {
                        if (best_sim[u] < best_sim[worst]) worst = u;
                    }
                    if (s > best_sim[worst]) {
                        best_sim[worst] = s;
                        best_idx[worst] = id;
                    }
                }
            }
        }
        for (std::uint32_t t = 0; t < k; ++t) {
            if (t < static_cast<std::uint32_t>(fill)) {
                knn_sim[row * k + t] = best_sim[t];
                knn_idx[row * k + t] = best_idx[t];
            } else {
                knn_sim[row * k + t] = -CUDART_INF_F;
                knn_idx[row * k + t] = kEmpty;
            }
        }
    }
}

__global__ void reverse_edges_kernel(const std::uint32_t* __restrict__ knn,
                                     std::uint32_t* __restrict__ rev,
                                     std::uint32_t* __restrict__ rev_count,
                                     std::uint32_t n,
                                     std::uint32_t k) {
    const std::uint32_t edge = blockIdx.x * blockDim.x + threadIdx.x;
    if (edge >= n * k) return;
    const std::uint32_t src = edge / k;
    const std::uint32_t dst = knn[edge];
    if (dst >= n) return;
    const std::uint32_t pos = atomicAdd(rev_count + dst, 1u);
    if (pos < k) rev[dst * k + pos] = src;
}

__global__ void prune_graph_kernel(const std::uint8_t* __restrict__ codes,
                                   const float* __restrict__ similarity,
                                   const std::uint32_t* __restrict__ knn,
                                   const std::uint32_t* __restrict__ rev,
                                   std::uint32_t* __restrict__ graph,
                                   std::uint32_t n,
                                   std::uint32_t code_dim,
                                   std::uint32_t pq_m,
                                   std::uint32_t ksub,
                                   std::uint32_t k_int,
                                   std::uint32_t degree) {
    const std::uint32_t row = blockIdx.x;
    if (row >= n) return;
    std::uint32_t candidates[2 * kMaxIntermediate];
    std::uint32_t n_cand = 0;
    auto push = [&](std::uint32_t id) {
        if (id >= n || id == row) return;
        for (std::uint32_t t = 0; t < n_cand; ++t) {
            if (candidates[t] == id) return;
        }
        if (n_cand < 2 * k_int) candidates[n_cand++] = id;
    };
    if (threadIdx.x == 0) {
        for (std::uint32_t t = 0; t < k_int; ++t) {
            push(knn[row * k_int + t]);
            push(rev[row * k_int + t]);
        }
        const std::uint8_t* row_codes = codes + static_cast<std::size_t>(row) * code_dim;
        float scores[2 * kMaxIntermediate];
        for (std::uint32_t t = 0; t < n_cand; ++t) {
            scores[t] = sdc_score(row_codes,
                                  codes + static_cast<std::size_t>(candidates[t]) * code_dim,
                                  similarity,
                                  code_dim,
                                  pq_m,
                                  ksub);
        }
        const std::uint32_t keep = degree < n_cand ? degree : n_cand;
        for (std::uint32_t t = 0; t < keep; ++t) {
            std::uint32_t best = t;
            for (std::uint32_t u = t + 1; u < n_cand; ++u) {
                if (scores[u] > scores[best]) best = u;
            }
            const float ts = scores[t];
            scores[t] = scores[best];
            scores[best] = ts;
            const std::uint32_t ti = candidates[t];
            candidates[t] = candidates[best];
            candidates[best] = ti;
            graph[row * degree + t] = candidates[t];
        }
        for (std::uint32_t t = keep; t < degree; ++t) graph[row * degree + t] = kEmpty;
    }
}

__global__ void random_graph_kernel(std::uint32_t* __restrict__ graph,
                                    std::uint32_t n,
                                    std::uint32_t degree,
                                    std::uint32_t seed) {
    const std::uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n) return;
    for (std::uint32_t k = 0; k < degree; ++k) {
        std::uint32_t x = (row + 1u) * 747796405u + (k + 1u) * 2891336453u + seed;
        x ^= x >> 16;
        x *= 0x7feb352du;
        x ^= x >> 15;
        x *= 0x846ca68bu;
        x ^= x >> 16;
        std::uint32_t neighbor = n > 1 ? (x % (n - 1)) : 0u;
        if (neighbor >= row) neighbor += 1u;
        for (std::uint32_t prev = 0; prev < k; ++prev) {
            if (graph[row * degree + prev] == neighbor) {
                neighbor = (neighbor + 1u) % n;
                if (neighbor == row) neighbor = (neighbor + 1u) % n;
            }
        }
        graph[row * degree + k] = neighbor;
    }
}

__global__ void nndescent_step_kernel(const std::uint8_t* __restrict__ codes,
                                      const float* __restrict__ similarity,
                                      const std::uint32_t* __restrict__ graph_in,
                                      std::uint32_t* __restrict__ graph_out,
                                      std::uint32_t n,
                                      std::uint32_t code_dim,
                                      std::uint32_t pq_m,
                                      std::uint32_t ksub,
                                      std::uint32_t degree) {
    const std::uint32_t row = blockIdx.x;
    if (row >= n) return;

    extern __shared__ char smem[];
    const std::uint32_t nwarps = blockDim.x / 32u;
    const std::uint32_t warp = threadIdx.x / 32u;
    const std::uint32_t lane = threadIdx.x & 31u;
    std::uint32_t* curr = reinterpret_cast<std::uint32_t*>(smem);
    float* all_sim = reinterpret_cast<float*>(curr + degree);
    std::uint32_t* all_idx = reinterpret_cast<std::uint32_t*>(all_sim + nwarps * degree);
    int* all_fill = reinterpret_cast<int*>(all_idx + nwarps * degree);
    float* warp_sim = all_sim + warp * degree;
    std::uint32_t* warp_idx = all_idx + warp * degree;
    int* warp_fill = all_fill + warp;

    for (std::uint32_t k = threadIdx.x; k < degree; k += blockDim.x) {
        curr[k] = graph_in[row * degree + k];
    }
    if (lane == 0) warp_fill[0] = 0;
    __syncthreads();

    const std::uint8_t* row_codes = codes + static_cast<std::size_t>(row) * code_dim;
    const std::uint32_t n_cand = degree + degree * degree;
    const std::uint32_t stride = nwarps * 32u;
    for (std::uint32_t base = 0; base < n_cand; base += stride) {
        const std::uint32_t e = base + warp * 32u + lane;
        float score = -CUDART_INF_F;
        std::uint32_t other = kEmpty;
        if (e < n_cand) {
            std::uint32_t candidate = kEmpty;
            if (e < degree) {
                candidate = curr[e];
            } else {
                const std::uint32_t e2 = e - degree;
                const std::uint32_t parent = curr[e2 / degree];
                if (parent < n) candidate = graph_in[parent * degree + (e2 % degree)];
            }
            if (candidate < n && candidate != row) {
                score = sdc_score(row_codes,
                                  codes + static_cast<std::size_t>(candidate) * code_dim,
                                  similarity,
                                  code_dim,
                                  pq_m,
                                  ksub);
                other = candidate;
            }
        }
        float lane_scores[32];
        std::uint32_t lane_ids[32];
        for (int src = 0; src < 32; ++src) {
            lane_scores[src] = __shfl_sync(0xffffffffu, score, src);
            lane_ids[src] = __shfl_sync(0xffffffffu, other, src);
        }
        if (lane == 0) {
            for (int src = 0; src < 32; ++src) {
                if (lane_ids[src] == kEmpty) continue;
                const float s = lane_scores[src];
                const std::uint32_t id = lane_ids[src];
                int fill = warp_fill[0];
                bool seen = false;
                for (int t = 0; t < fill; ++t) {
                    if (warp_idx[t] == id) {
                        seen = true;
                        break;
                    }
                }
                if (seen) continue;
                if (fill < static_cast<int>(degree)) {
                    warp_sim[fill] = s;
                    warp_idx[fill] = id;
                    warp_fill[0] = fill + 1;
                } else {
                    int worst = 0;
                    for (int t = 1; t < static_cast<int>(degree); ++t) {
                        if (warp_sim[t] < warp_sim[worst]) worst = t;
                    }
                    if (s > warp_sim[worst]) {
                        warp_sim[worst] = s;
                        warp_idx[worst] = id;
                    }
                }
            }
        }
        __syncwarp();
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        float best_sim[kMaxIntermediate];
        std::uint32_t best_idx[kMaxIntermediate];
        int fill = 0;
        for (std::uint32_t w = 0; w < nwarps; ++w) {
            float* src_sim = all_sim + w * degree;
            std::uint32_t* src_idx = all_idx + w * degree;
            const int src_fill = all_fill[w];
            for (int t = 0; t < src_fill; ++t) {
                const float s = src_sim[t];
                const std::uint32_t id = src_idx[t];
                bool seen = false;
                for (int u = 0; u < fill; ++u) {
                    if (best_idx[u] == id) {
                        seen = true;
                        break;
                    }
                }
                if (seen) continue;
                if (fill < static_cast<int>(degree)) {
                    best_sim[fill] = s;
                    best_idx[fill] = id;
                    ++fill;
                } else {
                    int worst = 0;
                    for (int u = 1; u < static_cast<int>(degree); ++u) {
                        if (best_sim[u] < best_sim[worst]) worst = u;
                    }
                    if (s > best_sim[worst]) {
                        best_sim[worst] = s;
                        best_idx[worst] = id;
                    }
                }
            }
        }
        for (std::uint32_t t = 0; t < degree; ++t) {
            graph_out[row * degree + t] = t < static_cast<std::uint32_t>(fill) ? best_idx[t] : kEmpty;
        }
    }
}

__global__ void search_kernel(const std::uint8_t* __restrict__ codes,
                              const std::uint32_t* __restrict__ graph,
                              const float* __restrict__ codebook,
                              const float* __restrict__ queries,
                              std::uint32_t* __restrict__ hashmap,
                              std::uint32_t* __restrict__ workspace_ids,
                              float* __restrict__ workspace_dists,
                              std::int32_t* __restrict__ out_labels,
                              float* __restrict__ out_distances,
                              std::uint32_t n,
                              std::uint32_t nq,
                              std::uint32_t window_size,
                              std::uint32_t pq_m,
                              std::uint32_t ksub,
                              std::uint32_t dsub,
                              std::uint32_t graph_degree,
                              std::uint32_t itopk,
                              std::uint32_t search_width,
                              std::uint32_t min_iterations,
                              std::uint32_t max_iterations,
                              std::uint32_t k,
                              std::uint32_t hash_bitlen,
                              std::uint32_t code_dim,
                              std::uint32_t result_padded) {
    const std::uint32_t query_id = blockIdx.x;
    if (query_id >= nq) return;

    extern __shared__ char smem[];
    __half* lut = reinterpret_cast<__half*>(smem);
    const std::uint32_t lut_elems = code_dim * ksub;
    const std::uint32_t lut_bytes = (lut_elems * static_cast<std::uint32_t>(sizeof(__half)) + 7u) & ~7u;
    auto* parents = reinterpret_cast<std::uint32_t*>(smem + lut_bytes);
    std::uint32_t* n_parents_ptr = parents + search_width;
    std::uint32_t* ids = workspace_ids + static_cast<std::size_t>(query_id) * result_padded;
    float* dists = workspace_dists + static_cast<std::size_t>(query_id) * result_padded;

    const std::uint32_t dim = pq_m * dsub;
    const float* query = queries + static_cast<std::size_t>(query_id) * window_size * dim;

    for (std::uint32_t idx = threadIdx.x; idx < lut_elems; idx += blockDim.x) {
        const std::uint32_t code = idx % ksub;
        const std::uint32_t sub = (idx / ksub) % pq_m;
        const std::uint32_t position = idx / (ksub * pq_m);
        const float* qsub = query + position * dim + sub * dsub;
        const float* centroid = codebook + (sub * ksub + code) * dsub;
        float acc = 0.0f;
        for (std::uint32_t d = 0; d < dsub; ++d) acc += qsub[d] * centroid[d];
        lut[idx] = __float2half(-acc);
    }

    std::uint32_t* visited = hashmap + static_cast<std::size_t>(query_id) * (1u << hash_bitlen);
    __syncthreads();

    for (std::uint32_t i = threadIdx.x; i < result_padded; i += blockDim.x) {
        ids[i] = kEmpty;
        dists[i] = CUDART_INF_F;
    }
    __syncthreads();

    // Random seeds, CAGRA-style xorshift pickup into the itopk slots.
    for (std::uint32_t slot = threadIdx.x; slot < itopk; slot += blockDim.x) {
        std::uint32_t seed = slot + query_id * 2654435761u;
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        const std::uint32_t index = (n == 0) ? 0u : (seed % n);
        if (hashmap::insert(visited, hash_bitlen, index)) {
            ids[slot] = index;
        } else {
            ids[slot] = kEmpty;
            dists[slot] = CUDART_INF_F;
        }
    }
    __syncthreads();

    const unsigned warp = threadIdx.x / 32u;
    const unsigned lane = threadIdx.x & 31u;
    const unsigned nwarps = blockDim.x / 32u;
    for (std::uint32_t slot = warp; slot < itopk; slot += nwarps) {
        const std::uint32_t index = ids[slot];
        float dist = CUDART_INF_F;
        if (index != kEmpty) {
            dist = warp_adc(lut, codes + static_cast<std::size_t>(index) * code_dim, code_dim, ksub);
        }
        if (lane == 0) dists[slot] = dist;
    }
    __syncthreads();
    bitonic_sort_kv(dists, ids, itopk, result_padded);
    __syncthreads();

    for (std::uint32_t iter = 0; iter < max_iterations; ++iter) {
        if (threadIdx.x == 0) {
            std::uint32_t count = 0;
            for (std::uint32_t i = 0; i < itopk && count < search_width; ++i) {
                const std::uint32_t raw = ids[i];
                if (raw == kEmpty) continue;
                if (raw & kCheckedBit) continue;
                parents[count++] = raw & kIndexMask;
                ids[i] = raw | kCheckedBit;
            }
            n_parents_ptr[0] = count;
        }
        __syncthreads();
        const std::uint32_t n_parents = n_parents_ptr[0];
        if (n_parents == 0 && iter >= min_iterations) break;
        if (n_parents == 0) {
            __syncthreads();
            continue;
        }

        const std::uint32_t n_expand = n_parents * graph_degree;
        for (std::uint32_t t = threadIdx.x; t < search_width * graph_degree; t += blockDim.x) {
            ids[itopk + t] = kEmpty;
            dists[itopk + t] = CUDART_INF_F;
        }
        __syncthreads();

        for (std::uint32_t t = threadIdx.x; t < n_expand; t += blockDim.x) {
            const std::uint32_t parent = parents[t / graph_degree];
            const std::uint32_t neigh = graph[parent * graph_degree + (t % graph_degree)];
            if (neigh >= n) continue;
            if (hashmap::insert(visited, hash_bitlen, neigh)) {
                ids[itopk + t] = neigh;
            }
        }
        __syncthreads();

        for (std::uint32_t t = warp; t < n_expand; t += nwarps) {
            const std::uint32_t index = ids[itopk + t];
            float dist = CUDART_INF_F;
            if (index != kEmpty) {
                dist = warp_adc(
                    lut, codes + static_cast<std::size_t>(index) * code_dim, code_dim, ksub);
            }
            if (lane == 0) dists[itopk + t] = dist;
        }
        __syncthreads();
        bitonic_sort_kv(dists, ids, itopk + n_expand, result_padded);
        __syncthreads();
    }

    for (std::uint32_t rank = threadIdx.x; rank < k; rank += blockDim.x) {
        if (rank < itopk && ids[rank] != kEmpty) {
            out_labels[query_id * k + rank] = static_cast<std::int32_t>(ids[rank] & kIndexMask);
            out_distances[query_id * k + rank] = dists[rank];
        } else {
            out_labels[query_id * k + rank] = -1;
            out_distances[query_id * k + rank] = CUDART_INF_F;
        }
    }
}

void validate_dims(std::uint32_t n,
                   std::uint32_t window_size,
                   std::uint32_t pq_m,
                   std::uint32_t nbits,
                   std::uint32_t dsub) {
    if (n < 2) throw std::runtime_error("CAGRA requires at least 2 vectors");
    if (window_size < 1 || pq_m < 1 || dsub < 1) {
        throw std::runtime_error("window_size, pq_m, and dsub must be positive");
    }
    if (nbits < 1 || nbits > 8) throw std::runtime_error("nbits must be in [1, 8]");
    if (n >= kIndexMask) throw std::runtime_error("dataset is too large for 31-bit ids");
}

}  // namespace

void build_similarity_table(const float* codebook,
                            std::uint32_t pq_m,
                            std::uint32_t ksub,
                            std::uint32_t dsub,
                            float* similarity) {
    for (std::uint32_t m = 0; m < pq_m; ++m) {
        const float* centroids = codebook + static_cast<std::size_t>(m) * ksub * dsub;
        float* table = similarity + static_cast<std::size_t>(m) * ksub * ksub;
        for (std::uint32_t i = 0; i < ksub; ++i) {
            for (std::uint32_t j = 0; j < ksub; ++j) {
                float acc = 0.0f;
                const float* a = centroids + i * dsub;
                const float* b = centroids + j * dsub;
                for (std::uint32_t d = 0; d < dsub; ++d) acc += a[d] * b[d];
                table[i * ksub + j] = acc;
            }
        }
    }
}

Index::Index(Index&& other) noexcept { *this = std::move(other); }

Index& Index::operator=(Index&& other) noexcept {
    if (this == &other) return *this;
    release();
    n_ = other.n_;
    window_size_ = other.window_size_;
    pq_m_ = other.pq_m_;
    nbits_ = other.nbits_;
    ksub_ = other.ksub_;
    dsub_ = other.dsub_;
    code_dim_ = other.code_dim_;
    graph_degree_ = other.graph_degree_;
    host_codes_ = std::move(other.host_codes_);
    host_codebook_ = std::move(other.host_codebook_);
    host_graph_ = std::move(other.host_graph_);
    d_codes_ = other.d_codes_;
    d_codebook_ = other.d_codebook_;
    d_graph_ = other.d_graph_;
    other.d_codes_ = nullptr;
    other.d_codebook_ = nullptr;
    other.d_graph_ = nullptr;
    other.n_ = 0;
    return *this;
}

Index::~Index() { release(); }

void Index::release() const noexcept {
    if (d_codes_) cudaFree(d_codes_);
    if (d_codebook_) cudaFree(d_codebook_);
    if (d_graph_) cudaFree(d_graph_);
    d_codes_ = nullptr;
    d_codebook_ = nullptr;
    d_graph_ = nullptr;
}

std::size_t Index::lut_bytes() const {
    return static_cast<std::size_t>(code_dim_) * ksub_ * sizeof(__half);
}

void Index::upload_from_host() const {
    if (d_codes_ && d_codebook_ && d_graph_) return;
    release();
    const std::size_t code_bytes = host_codes_.size() * sizeof(std::uint8_t);
    const std::size_t codebook_bytes = host_codebook_.size() * sizeof(float);
    const std::size_t graph_bytes = host_graph_.size() * sizeof(std::uint32_t);
    PQ_CAGRA_CHECK(cudaMalloc(&d_codes_, code_bytes));
    PQ_CAGRA_CHECK(cudaMalloc(&d_codebook_, codebook_bytes));
    PQ_CAGRA_CHECK(cudaMalloc(&d_graph_, graph_bytes));
    PQ_CAGRA_CHECK(cudaMemcpy(d_codes_, host_codes_.data(), code_bytes, cudaMemcpyHostToDevice));
    PQ_CAGRA_CHECK(
        cudaMemcpy(d_codebook_, host_codebook_.data(), codebook_bytes, cudaMemcpyHostToDevice));
    PQ_CAGRA_CHECK(cudaMemcpy(d_graph_, host_graph_.data(), graph_bytes, cudaMemcpyHostToDevice));
}

void Index::build(const std::uint8_t* codes,
                  const float* codebook,
                  std::uint32_t n,
                  std::uint32_t window_size,
                  std::uint32_t pq_m,
                  std::uint32_t nbits,
                  std::uint32_t dsub,
                  const BuildParams& params) {
    report_progress("start", 0);
    validate_dims(n, window_size, pq_m, nbits, dsub);
    if (params.graph_degree < 1 || params.intermediate_graph_degree < params.graph_degree) {
        throw std::runtime_error("graph degrees are invalid");
    }
    if (params.graph_degree > kMaxGraphDegree || params.intermediate_graph_degree > kMaxIntermediate) {
        throw std::runtime_error("graph degree exceeds compiled maximum");
    }
    const std::uint32_t ksub = 1u << nbits;
    const std::uint32_t code_dim = window_size * pq_m;
    const std::uint32_t k_int = std::min(params.intermediate_graph_degree, n - 1);
    const std::uint32_t degree = std::min(params.graph_degree, n - 1);

    n_ = n;
    window_size_ = window_size;
    pq_m_ = pq_m;
    nbits_ = nbits;
    ksub_ = ksub;
    dsub_ = dsub;
    code_dim_ = code_dim;
    graph_degree_ = degree;

    host_codes_.assign(codes, codes + static_cast<std::size_t>(n) * code_dim);
    host_codebook_.assign(codebook, codebook + static_cast<std::size_t>(pq_m) * ksub * dsub);
    host_graph_.assign(static_cast<std::size_t>(n) * degree, kEmpty);

    std::vector<float> similarity(static_cast<std::size_t>(pq_m) * ksub * ksub);
    build_similarity_table(codebook, pq_m, ksub, dsub, similarity.data());
    report_progress("host_data_ready", 10);

    std::uint8_t* d_codes = nullptr;
    float* d_sim = nullptr;
    std::uint32_t* d_knn = nullptr;
    float* d_knn_sim = nullptr;
    std::uint32_t* d_knn_alt = nullptr;
    std::uint32_t* d_rev = nullptr;
    std::uint32_t* d_rev_count = nullptr;
    std::uint32_t* d_graph = nullptr;
    const bool use_nndescent =
        params.nndescent_iterations > 0 || n > kBruteForceLimit;
    const std::uint32_t nd_iters =
        params.nndescent_iterations > 0 ? params.nndescent_iterations : 8u;
    try {
        PQ_CAGRA_CHECK(cudaMalloc(&d_codes, host_codes_.size()));
        PQ_CAGRA_CHECK(cudaMalloc(&d_sim, similarity.size() * sizeof(float)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_knn, static_cast<std::size_t>(n) * k_int * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_rev, static_cast<std::size_t>(n) * k_int * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_rev_count, static_cast<std::size_t>(n) * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_graph, static_cast<std::size_t>(n) * degree * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMemcpy(d_codes, host_codes_.data(), host_codes_.size(), cudaMemcpyHostToDevice));
        PQ_CAGRA_CHECK(cudaMemcpy(
            d_sim, similarity.data(), similarity.size() * sizeof(float), cudaMemcpyHostToDevice));
        PQ_CAGRA_CHECK(cudaMemset(d_rev, 0xff, static_cast<std::size_t>(n) * k_int * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMemset(d_rev_count, 0, static_cast<std::size_t>(n) * sizeof(std::uint32_t)));
        report_progress("device_data_ready", 15);

        constexpr std::uint32_t threads = 256;
        const std::uint32_t nwarps = threads / 32u;
        if (use_nndescent) {
            PQ_CAGRA_CHECK(
                cudaMalloc(&d_knn_alt, static_cast<std::size_t>(n) * k_int * sizeof(std::uint32_t)));
            random_graph_kernel<<<(n + threads - 1) / threads, threads>>>(d_knn, n, k_int, 1u);
            check_launch();
            report_progress("random_graph_ready", 20, "method=nndescent");
            const std::size_t nd_smem =
                k_int * sizeof(std::uint32_t) +
                nwarps * k_int * (sizeof(float) + sizeof(std::uint32_t)) + nwarps * sizeof(int);
            std::uint32_t* src = d_knn;
            std::uint32_t* dst = d_knn_alt;
            for (std::uint32_t iter = 0; iter < nd_iters; ++iter) {
                nndescent_step_kernel<<<n, threads, nd_smem>>>(
                    d_codes, d_sim, src, dst, n, code_dim, pq_m, ksub, k_int);
                check_launch();
                const auto percent =
                    20u + static_cast<std::uint32_t>(
                               (50ull * static_cast<std::uint64_t>(iter + 1)) / nd_iters);
                report_iteration_progress("nndescent", percent, iter + 1, nd_iters);
                std::uint32_t* tmp = src;
                src = dst;
                dst = tmp;
            }
            if (src != d_knn) {
                PQ_CAGRA_CHECK(cudaMemcpy(d_knn,
                                          src,
                                          static_cast<std::size_t>(n) * k_int * sizeof(std::uint32_t),
                                          cudaMemcpyDeviceToDevice));
            }
            report_progress("knn_ready", 70, "method=nndescent");
        } else {
            PQ_CAGRA_CHECK(
                cudaMalloc(&d_knn_sim, static_cast<std::size_t>(n) * k_int * sizeof(float)));
            const std::size_t knn_smem =
                nwarps * k_int * (sizeof(float) + sizeof(std::uint32_t)) + nwarps * sizeof(int);
            knn_sdc_kernel<<<n, threads, knn_smem>>>(
                d_codes, d_sim, d_knn, d_knn_sim, n, code_dim, pq_m, ksub, k_int);
            check_launch();
            report_progress("knn_ready", 70, "method=brute_force");
        }

        const std::uint32_t reverse_threads = 256;
        const std::uint32_t reverse_blocks =
            static_cast<std::uint32_t>((static_cast<std::uint64_t>(n) * k_int + reverse_threads - 1) /
                                       reverse_threads);
        reverse_edges_kernel<<<reverse_blocks, reverse_threads>>>(d_knn, d_rev, d_rev_count, n, k_int);
        check_launch();
        report_progress("reverse_edges", 80);

        prune_graph_kernel<<<n, 32>>>(
            d_codes, d_sim, d_knn, d_rev, d_graph, n, code_dim, pq_m, ksub, k_int, degree);
        check_launch();
        report_progress("prune_graph", 90);

        PQ_CAGRA_CHECK(cudaMemcpy(host_graph_.data(),
                                  d_graph,
                                  static_cast<std::size_t>(n) * degree * sizeof(std::uint32_t),
                                  cudaMemcpyDeviceToHost));
        report_progress("graph_copied_to_host", 95);
    } catch (...) {
        cudaFree(d_codes);
        cudaFree(d_sim);
        cudaFree(d_knn);
        cudaFree(d_knn_sim);
        cudaFree(d_knn_alt);
        cudaFree(d_rev);
        cudaFree(d_rev_count);
        cudaFree(d_graph);
        throw;
    }
    cudaFree(d_codes);
    cudaFree(d_sim);
    cudaFree(d_knn);
    cudaFree(d_knn_sim);
    cudaFree(d_knn_alt);
    cudaFree(d_rev);
    cudaFree(d_rev_count);
    cudaFree(d_graph);
    upload_from_host();
    report_progress("complete", 100);
}

void Index::search(const float* queries,
                   std::uint32_t nq,
                   const SearchParams& params,
                   std::int32_t* labels,
                   float* distances) const {
    if (host_codes_.empty() || host_graph_.empty() || host_codebook_.empty()) {
        throw std::runtime_error("index has no host data; call build() or load() first");
    }
    upload_from_host();
    if (nq < 1) throw std::runtime_error("search requires at least one query");
    if (params.k < 1 || params.itopk_size < params.k) {
        throw std::runtime_error("itopk_size must be >= k");
    }
    if (params.itopk_size > kMaxItopk) throw std::runtime_error("itopk_size exceeds compiled maximum");
    if (params.search_width < 1 || params.max_iterations < 1) {
        throw std::runtime_error("search_width and max_iterations must be positive");
    }
    if (params.hash_bitlen < 8 || params.hash_bitlen > 20) {
        throw std::runtime_error("hash_bitlen must be in [8, 20]");
    }
    if (params.search_width > kMaxSearchWidth) {
        throw std::runtime_error("search_width exceeds compiled maximum");
    }
    const std::uint32_t itopk = std::min(params.itopk_size, n_);
    const std::uint32_t k = std::min(params.k, n_);
    const std::uint32_t search_width = std::min(params.search_width, itopk);
    const std::uint32_t result_count = itopk + search_width * graph_degree_;
    const std::uint32_t result_padded = round_up_pow2(result_count);
    const std::uint32_t lut_elems = code_dim_ * ksub_;
    const std::uint32_t lut_bytes =
        (lut_elems * static_cast<std::uint32_t>(sizeof(__half)) + 7u) & ~7u;
    const std::uint32_t smem =
        lut_bytes + search_width * sizeof(std::uint32_t) + sizeof(std::uint32_t);
    if (smem > kMaxSmemBytes) {
        throw std::runtime_error("shared-memory LUT does not fit: " + std::to_string(smem) +
                                 " bytes required");
    }

    float* d_queries = nullptr;
    std::int32_t* d_labels = nullptr;
    float* d_distances = nullptr;
    std::uint32_t* d_hash = nullptr;
    std::uint32_t* d_ws_ids = nullptr;
    float* d_ws_dists = nullptr;
    const std::size_t query_elems = static_cast<std::size_t>(nq) * window_size_ * dim();
    const std::size_t hash_elems = static_cast<std::size_t>(nq) * (1u << params.hash_bitlen);
    const std::size_t workspace = static_cast<std::size_t>(nq) * result_padded;
    try {
        PQ_CAGRA_CHECK(cudaMalloc(&d_queries, query_elems * sizeof(float)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_labels, static_cast<std::size_t>(nq) * k * sizeof(std::int32_t)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_distances, static_cast<std::size_t>(nq) * k * sizeof(float)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_hash, hash_elems * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_ws_ids, workspace * sizeof(std::uint32_t)));
        PQ_CAGRA_CHECK(cudaMalloc(&d_ws_dists, workspace * sizeof(float)));
        PQ_CAGRA_CHECK(cudaMemcpy(
            d_queries, queries, query_elems * sizeof(float), cudaMemcpyHostToDevice));
        PQ_CAGRA_CHECK(cudaMemset(d_hash, 0xff, hash_elems * sizeof(std::uint32_t)));

        constexpr int threads = 256;
        PQ_CAGRA_CHECK(cudaFuncSetAttribute(
            search_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(smem)));
        search_kernel<<<nq, threads, smem>>>(d_codes_,
                                             d_graph_,
                                             d_codebook_,
                                             d_queries,
                                             d_hash,
                                             d_ws_ids,
                                             d_ws_dists,
                                             d_labels,
                                             d_distances,
                                             n_,
                                             nq,
                                             window_size_,
                                             pq_m_,
                                             ksub_,
                                             dsub_,
                                             graph_degree_,
                                             itopk,
                                             search_width,
                                             params.min_iterations,
                                             params.max_iterations,
                                             k,
                                             params.hash_bitlen,
                                             code_dim_,
                                             result_padded);
        check_launch();
        PQ_CAGRA_CHECK(cudaMemcpy(
            labels, d_labels, static_cast<std::size_t>(nq) * k * sizeof(std::int32_t), cudaMemcpyDeviceToHost));
        PQ_CAGRA_CHECK(cudaMemcpy(distances,
                                  d_distances,
                                  static_cast<std::size_t>(nq) * k * sizeof(float),
                                  cudaMemcpyDeviceToHost));
    } catch (...) {
        cudaFree(d_queries);
        cudaFree(d_labels);
        cudaFree(d_distances);
        cudaFree(d_hash);
        cudaFree(d_ws_ids);
        cudaFree(d_ws_dists);
        throw;
    }
    cudaFree(d_queries);
    cudaFree(d_labels);
    cudaFree(d_distances);
    cudaFree(d_hash);
    cudaFree(d_ws_ids);
    cudaFree(d_ws_dists);
}

}  // namespace ginflow::pq_cagra
