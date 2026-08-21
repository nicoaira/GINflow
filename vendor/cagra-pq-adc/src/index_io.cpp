/*
 * Host-side save/load for the PQ-CAGRA graph. Shared by GPU and CPU builds.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "cagra_pq_adc.hpp"

#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ginflow::pq_cagra {

void Index::save(const std::string& path) const {
    if (host_codes_.empty() || host_graph_.empty()) {
        throw std::runtime_error("cannot save an empty index");
    }
    IndexHeader header{};
    std::memcpy(header.magic, "GINPCG01", 8);
    header.n = n_;
    header.window_size = window_size_;
    header.pq_m = pq_m_;
    header.nbits = nbits_;
    header.dsub = dsub_;
    header.graph_degree = graph_degree_;
    header.code_dim = code_dim_;
    header.reserved = 0;
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot create index file: " + path);
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));
    out.write(reinterpret_cast<const char*>(host_codebook_.data()),
              static_cast<std::streamsize>(host_codebook_.size() * sizeof(float)));
    out.write(reinterpret_cast<const char*>(host_codes_.data()),
              static_cast<std::streamsize>(host_codes_.size()));
    out.write(reinterpret_cast<const char*>(host_graph_.data()),
              static_cast<std::streamsize>(host_graph_.size() * sizeof(std::uint32_t)));
    if (!out) throw std::runtime_error("failed writing index file: " + path);
}

Index Index::load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open index file: " + path);
    IndexHeader header{};
    in.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!in || std::string(header.magic, header.magic + 8) != "GINPCG01") {
        throw std::runtime_error("invalid PQ-CAGRA index header");
    }
    Index index;
    index.n_ = header.n;
    index.window_size_ = header.window_size;
    index.pq_m_ = header.pq_m;
    index.nbits_ = header.nbits;
    index.ksub_ = 1u << header.nbits;
    index.dsub_ = header.dsub;
    index.graph_degree_ = header.graph_degree;
    index.code_dim_ = header.code_dim;
    index.host_codebook_.resize(static_cast<std::size_t>(index.pq_m_) * index.ksub_ * index.dsub_);
    index.host_codes_.resize(static_cast<std::size_t>(index.n_) * index.code_dim_);
    index.host_graph_.resize(static_cast<std::size_t>(index.n_) * index.graph_degree_);
    in.read(reinterpret_cast<char*>(index.host_codebook_.data()),
            static_cast<std::streamsize>(index.host_codebook_.size() * sizeof(float)));
    in.read(reinterpret_cast<char*>(index.host_codes_.data()),
            static_cast<std::streamsize>(index.host_codes_.size()));
    in.read(reinterpret_cast<char*>(index.host_graph_.data()),
            static_cast<std::streamsize>(index.host_graph_.size() * sizeof(std::uint32_t)));
    if (!in) throw std::runtime_error("truncated PQ-CAGRA index file");
    return index;
}

}  // namespace ginflow::pq_cagra
