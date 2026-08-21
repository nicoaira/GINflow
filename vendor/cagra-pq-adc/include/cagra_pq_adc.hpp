#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace ginflow::pq_cagra {

struct BuildParams {
    std::uint32_t graph_degree = 32;
    std::uint32_t intermediate_graph_degree = 64;
    // 0 = brute-force kNN when n is small, otherwise NN-Descent (8 iters).
    std::uint32_t nndescent_iterations = 0;
};

struct SearchParams {
    std::uint32_t k = 10;
    std::uint32_t itopk_size = 64;
    std::uint32_t search_width = 1;
    std::uint32_t min_iterations = 8;
    std::uint32_t max_iterations = 64;
    std::uint32_t hash_bitlen = 16;
    std::uint32_t num_threads = 0;
};

struct IndexHeader {
    char magic[8];
    std::uint32_t n;
    std::uint32_t window_size;
    std::uint32_t pq_m;
    std::uint32_t nbits;
    std::uint32_t dsub;
    std::uint32_t graph_degree;
    std::uint32_t code_dim;
    std::uint32_t reserved;
};

class Index {
public:
    Index() = default;
    Index(const Index&) = delete;
    Index& operator=(const Index&) = delete;
    Index(Index&&) noexcept;
    Index& operator=(Index&&) noexcept;
    ~Index();

    void build(const std::uint8_t* codes,
               const float* codebook,
               std::uint32_t n,
               std::uint32_t window_size,
               std::uint32_t pq_m,
               std::uint32_t nbits,
               std::uint32_t dsub,
               const BuildParams& params);

    void search(const float* queries,
                std::uint32_t nq,
                const SearchParams& params,
                std::int32_t* labels,
                float* distances) const;

    void search_cpu(const float* queries,
                    std::uint32_t nq,
                    const SearchParams& params,
                    std::int32_t* labels,
                    float* distances) const;

    void save(const std::string& path) const;
    static Index load(const std::string& path);

    std::uint32_t size() const { return n_; }
    std::uint32_t window_size() const { return window_size_; }
    std::uint32_t pq_m() const { return pq_m_; }
    std::uint32_t nbits() const { return nbits_; }
    std::uint32_t ksub() const { return ksub_; }
    std::uint32_t dsub() const { return dsub_; }
    std::uint32_t dim() const { return pq_m_ * dsub_; }
    std::uint32_t code_dim() const { return code_dim_; }
    std::uint32_t graph_degree() const { return graph_degree_; }

    const std::vector<std::uint8_t>& host_codes() const { return host_codes_; }
    const std::vector<float>& host_codebook() const { return host_codebook_; }
    const std::vector<std::uint32_t>& host_graph() const { return host_graph_; }

private:
    void release() const noexcept;
    void upload_from_host() const;
    std::size_t lut_bytes() const;

    std::uint32_t n_ = 0;
    std::uint32_t window_size_ = 0;
    std::uint32_t pq_m_ = 0;
    std::uint32_t nbits_ = 0;
    std::uint32_t ksub_ = 0;
    std::uint32_t dsub_ = 0;
    std::uint32_t code_dim_ = 0;
    std::uint32_t graph_degree_ = 0;

    std::vector<std::uint8_t> host_codes_;
    std::vector<float> host_codebook_;
    std::vector<std::uint32_t> host_graph_;

    mutable std::uint8_t* d_codes_ = nullptr;
    mutable float* d_codebook_ = nullptr;
    mutable std::uint32_t* d_graph_ = nullptr;
};

void build_similarity_table(const float* codebook,
                            std::uint32_t pq_m,
                            std::uint32_t ksub,
                            std::uint32_t dsub,
                            float* similarity);

}  // namespace ginflow::pq_cagra
