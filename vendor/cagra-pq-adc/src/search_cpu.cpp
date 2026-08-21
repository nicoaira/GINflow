/*
 * CPU beam search over a serialized PQ-CAGRA graph.
 *
 * The adjacency list is the GPU-built CAGRA graph.  Distances use the same
 * node-level PQ ADC as the CUDA kernel (query stays in original space).
 * This is the CPU counterpart of cuVS CAGRA-to-HNSW: the graph is a single
 * base layer, not an upstream hnswlib hierarchy.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "cagra_pq_adc.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ginflow::pq_cagra {
namespace {

constexpr std::uint32_t kCheckedBit = 0x80000000u;
constexpr std::uint32_t kIndexMask = 0x7fffffffu;
constexpr std::uint32_t kEmpty = 0xffffffffu;

struct HostHash {
    std::vector<std::uint32_t> table;
    std::uint32_t mask = 0;

    explicit HostHash(std::uint32_t bitlen)
        : table(std::size_t{1} << bitlen, kEmpty), mask((1u << bitlen) - 1u) {}

    bool insert(std::uint32_t key) {
        std::uint32_t index = (key ^ (key >> 16)) & mask;
        for (std::size_t probe = 0; probe < table.size(); ++probe) {
            if (table[index] == kEmpty) {
                table[index] = key;
                return true;
            }
            if (table[index] == key) return false;
            index = (index + 1u) & mask;
        }
        return false;
    }
};

void build_lut(const float* query,
               const float* codebook,
               float* lut,
               std::uint32_t window_size,
               std::uint32_t pq_m,
               std::uint32_t ksub,
               std::uint32_t dsub) {
    const std::uint32_t dim = pq_m * dsub;
    for (std::uint32_t position = 0; position < window_size; ++position) {
        for (std::uint32_t sub = 0; sub < pq_m; ++sub) {
            const float* qsub = query + position * dim + sub * dsub;
            for (std::uint32_t code = 0; code < ksub; ++code) {
                const float* centroid = codebook + (sub * ksub + code) * dsub;
                float acc = 0.0f;
                for (std::uint32_t d = 0; d < dsub; ++d) acc += qsub[d] * centroid[d];
                lut[(position * pq_m + sub) * ksub + code] = -acc;
            }
        }
    }
}

float adc_distance(const float* lut,
                   const std::uint8_t* codes,
                   std::uint32_t code_dim,
                   std::uint32_t ksub) {
    float sum = 0.0f;
    for (std::uint32_t i = 0; i < code_dim; ++i) {
        sum += lut[i * ksub + codes[i]];
    }
    return sum;
}

std::uint32_t xorshift(std::uint32_t seed) {
    seed ^= seed << 13;
    seed ^= seed >> 17;
    seed ^= seed << 5;
    return seed;
}

void search_one(const std::uint8_t* codes,
                const float* codebook,
                const std::uint32_t* graph,
                std::uint32_t n,
                std::uint32_t window_size,
                std::uint32_t pq_m,
                std::uint32_t ksub,
                std::uint32_t dsub,
                std::uint32_t code_dim,
                std::uint32_t graph_degree,
                const float* query,
                std::uint32_t query_id,
                const SearchParams& params,
                std::int32_t* labels,
                float* distances) {
    const std::uint32_t itopk = std::min(params.itopk_size, n);
    const std::uint32_t k = std::min(params.k, n);
    const std::uint32_t search_width = std::min(params.search_width, itopk);
    const std::uint32_t result_count = itopk + search_width * graph_degree;
    std::vector<float> lut(static_cast<std::size_t>(code_dim) * ksub);
    build_lut(query, codebook, lut.data(), window_size, pq_m, ksub, dsub);

    std::vector<float> dists(result_count, std::numeric_limits<float>::infinity());
    std::vector<std::uint32_t> ids(result_count, kEmpty);
    HostHash visited(params.hash_bitlen);

    for (std::uint32_t slot = 0; slot < itopk; ++slot) {
        std::uint32_t seed = slot + query_id * 2654435761u;
        seed = xorshift(seed);
        const std::uint32_t index = n == 0 ? 0u : (seed % n);
        if (visited.insert(index)) {
            ids[slot] = index;
            dists[slot] = adc_distance(lut.data(), codes + static_cast<std::size_t>(index) * code_dim, code_dim, ksub);
        }
    }
    const auto by_dist = [&](std::uint32_t a, std::uint32_t b) {
        if (dists[a] != dists[b]) return dists[a] < dists[b];
        return (ids[a] & kIndexMask) < (ids[b] & kIndexMask);
    };
    std::vector<std::uint32_t> order(itopk);
    const auto sort_prefix = [&](std::uint32_t count) {
        order.resize(count);
        for (std::uint32_t i = 0; i < count; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), by_dist);
        std::vector<float> nd(count);
        std::vector<std::uint32_t> ni(count);
        for (std::uint32_t i = 0; i < count; ++i) {
            nd[i] = dists[order[i]];
            ni[i] = ids[order[i]];
        }
        for (std::uint32_t i = 0; i < count; ++i) {
            dists[i] = nd[i];
            ids[i] = ni[i];
        }
    };
    sort_prefix(itopk);

    for (std::uint32_t iter = 0; iter < params.max_iterations; ++iter) {
        std::uint32_t parents[32];
        std::uint32_t n_parents = 0;
        for (std::uint32_t i = 0; i < itopk && n_parents < search_width; ++i) {
            const std::uint32_t raw = ids[i];
            if (raw == kEmpty) continue;
            if (raw & kCheckedBit) continue;
            parents[n_parents++] = raw & kIndexMask;
            ids[i] = raw | kCheckedBit;
        }
        if (n_parents == 0 && iter >= params.min_iterations) break;
        if (n_parents == 0) continue;

        const std::uint32_t n_expand = n_parents * graph_degree;
        for (std::uint32_t t = 0; t < n_expand; ++t) {
            ids[itopk + t] = kEmpty;
            dists[itopk + t] = std::numeric_limits<float>::infinity();
        }
        for (std::uint32_t t = 0; t < n_expand; ++t) {
            const std::uint32_t parent = parents[t / graph_degree];
            const std::uint32_t neigh = graph[parent * graph_degree + (t % graph_degree)];
            if (neigh >= n) continue;
            if (!visited.insert(neigh)) continue;
            ids[itopk + t] = neigh;
            dists[itopk + t] = adc_distance(
                lut.data(), codes + static_cast<std::size_t>(neigh) * code_dim, code_dim, ksub);
        }
        sort_prefix(itopk + n_expand);
    }

    for (std::uint32_t rank = 0; rank < k; ++rank) {
        if (rank < itopk && ids[rank] != kEmpty) {
            labels[rank] = static_cast<std::int32_t>(ids[rank] & kIndexMask);
            distances[rank] = dists[rank];
        } else {
            labels[rank] = -1;
            distances[rank] = std::numeric_limits<float>::infinity();
        }
    }
}

}  // namespace

void Index::search_cpu(const float* queries,
                       std::uint32_t nq,
                       const SearchParams& params,
                       std::int32_t* labels,
                       float* distances) const {
    if (host_codes_.empty() || host_graph_.empty() || host_codebook_.empty()) {
        throw std::runtime_error("index has no host data; call build() or load() first");
    }
    if (nq < 1) throw std::runtime_error("search requires at least one query");
    if (params.k < 1 || params.itopk_size < params.k) {
        throw std::runtime_error("itopk_size must be >= k");
    }
    if (params.search_width < 1 || params.search_width > 32 || params.max_iterations < 1) {
        throw std::runtime_error("search_width must be in [1, 32] and max_iterations positive");
    }
    if (params.hash_bitlen < 8 || params.hash_bitlen > 20) {
        throw std::runtime_error("hash_bitlen must be in [8, 20]");
    }
    const std::uint32_t k = std::min(params.k, n_);
    const std::uint32_t dim = pq_m_ * dsub_;
#ifdef _OPENMP
    const int threads = params.num_threads > 0 ? static_cast<int>(params.num_threads) : omp_get_max_threads();
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
    for (std::int64_t q = 0; q < static_cast<std::int64_t>(nq); ++q) {
        search_one(host_codes_.data(),
                   host_codebook_.data(),
                   host_graph_.data(),
                   n_,
                   window_size_,
                   pq_m_,
                   ksub_,
                   dsub_,
                   code_dim_,
                   graph_degree_,
                   queries + static_cast<std::size_t>(q) * window_size_ * dim,
                   static_cast<std::uint32_t>(q),
                   params,
                   labels + static_cast<std::size_t>(q) * k,
                   distances + static_cast<std::size_t>(q) * k);
    }
}

}  // namespace ginflow::pq_cagra
