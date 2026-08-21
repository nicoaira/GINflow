/*
 * CPU-only stubs for GPU build/search. Load + search_cpu still work.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "cagra_pq_adc.hpp"

#include <stdexcept>
#include <utility>

namespace ginflow::pq_cagra {

Index::Index(Index&& other) noexcept { *this = std::move(other); }

Index& Index::operator=(Index&& other) noexcept {
    if (this == &other) return *this;
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
    d_codes_ = nullptr;
    d_codebook_ = nullptr;
    d_graph_ = nullptr;
    other.n_ = 0;
    other.d_codes_ = nullptr;
    other.d_codebook_ = nullptr;
    other.d_graph_ = nullptr;
    return *this;
}

Index::~Index() = default;

void Index::release() const noexcept {}

void Index::upload_from_host() const {}

std::size_t Index::lut_bytes() const { return 0; }

void Index::build(const std::uint8_t*,
                  const float*,
                  std::uint32_t,
                  std::uint32_t,
                  std::uint32_t,
                  std::uint32_t,
                  std::uint32_t,
                  const BuildParams&) {
    throw std::runtime_error(
        "PQ-CAGRA graph construction needs the GPU package nicolas.aira::pq-cagra-adc");
}

void Index::search(const float*,
                   std::uint32_t,
                   const SearchParams&,
                   std::int32_t*,
                   float*) const {
    throw std::runtime_error(
        "GPU PQ-CAGRA search needs nicolas.aira::pq-cagra-adc; use search_cpu() or --search_device cpu");
}

}  // namespace ginflow::pq_cagra
