/*
 * SPDX-License-Identifier: Apache-2.0
 */
#include "cagra_pq_adc.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <utility>

namespace py = pybind11;
using ginflow::pq_cagra::BuildParams;
using ginflow::pq_cagra::Index;
using ginflow::pq_cagra::SearchParams;

namespace {

py::array_t<std::uint8_t> as_c_uint8(const py::array& values) {
    return py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast>::ensure(values);
}

py::array_t<float> as_c_float(const py::array& values) {
    return py::array_t<float, py::array::c_style | py::array::forcecast>::ensure(values);
}

void build_index(Index& index,
                 const py::array& codes_in,
                 const py::array& codebook_in,
                 std::uint32_t window_size,
                 std::uint32_t pq_m,
                 std::uint32_t nbits,
                 std::uint32_t dsub,
                 std::uint32_t graph_degree,
                 std::uint32_t intermediate_graph_degree,
                 std::uint32_t nndescent_iterations) {
    auto codes = as_c_uint8(codes_in);
    auto codebook = as_c_float(codebook_in);
    if (codes.ndim() == 3) {
        if (window_size == 0) window_size = static_cast<std::uint32_t>(codes.shape(1));
        if (pq_m == 0) pq_m = static_cast<std::uint32_t>(codes.shape(2));
        if (codes.shape(1) != window_size || codes.shape(2) != pq_m) {
            throw std::runtime_error("codes shape must be (n, window_size, pq_m)");
        }
    } else if (codes.ndim() != 2) {
        throw std::runtime_error("codes must be 2-D (n, window_size * pq_m) or 3-D");
    }
    if (codebook.ndim() != 3) {
        throw std::runtime_error("codebook must have shape (pq_m, ksub, dsub)");
    }
    const std::uint32_t n = static_cast<std::uint32_t>(codes.shape(0));
    if (pq_m == 0) pq_m = static_cast<std::uint32_t>(codebook.shape(0));
    const std::uint32_t ksub = static_cast<std::uint32_t>(codebook.shape(1));
    if (dsub == 0) dsub = static_cast<std::uint32_t>(codebook.shape(2));
    if (nbits == 0) {
        std::uint32_t bits = 0;
        std::uint32_t value = ksub;
        while (value > 1) {
            if (value & 1u) throw std::runtime_error("ksub must be a power of two");
            value >>= 1;
            ++bits;
        }
        nbits = bits;
    }
    if ((1u << nbits) != ksub) throw std::runtime_error("codebook ksub does not match nbits");
    if (codebook.shape(0) != pq_m || codebook.shape(2) != dsub) {
        throw std::runtime_error("codebook shape does not match pq_m/dsub");
    }
    const std::uint32_t code_dim = window_size * pq_m;
    if (codes.ndim() == 2 && codes.shape(1) != code_dim) {
        throw std::runtime_error("codes second dimension must equal window_size * pq_m");
    }
    BuildParams params;
    params.graph_degree = graph_degree;
    params.intermediate_graph_degree = intermediate_graph_degree;
    params.nndescent_iterations = nndescent_iterations;
    py::gil_scoped_release release;
    index.build(codes.data(), codebook.data(), n, window_size, pq_m, nbits, dsub, params);
}

SearchParams make_search_params(std::uint32_t k,
                                std::uint32_t itopk_size,
                                std::uint32_t search_width,
                                std::uint32_t min_iterations,
                                std::uint32_t max_iterations,
                                std::uint32_t hash_bitlen,
                                std::uint32_t num_threads) {
    SearchParams params;
    params.k = k;
    params.itopk_size = itopk_size;
    params.search_width = search_width;
    params.min_iterations = min_iterations;
    params.max_iterations = max_iterations;
    params.hash_bitlen = hash_bitlen;
    params.num_threads = num_threads;
    return params;
}

py::array_t<float> require_query_windows(const Index& index, const py::array& queries_in) {
    auto queries = as_c_float(queries_in);
    if (queries.ndim() != 3) {
        throw std::runtime_error("queries must have shape (nq, window_size, dim)");
    }
    if (queries.shape(1) != index.window_size() || queries.shape(2) != index.dim()) {
        throw std::runtime_error("query window/dimension does not match the index");
    }
    return queries;
}

std::pair<py::array_t<std::int32_t>, py::array_t<float>> run_search(
    const Index& index,
    const py::array& queries_in,
    std::uint32_t k,
    std::uint32_t itopk_size,
    std::uint32_t search_width,
    std::uint32_t min_iterations,
    std::uint32_t max_iterations,
    std::uint32_t hash_bitlen,
    std::uint32_t num_threads,
    bool cpu) {
    auto queries = require_query_windows(index, queries_in);
    const std::uint32_t nq = static_cast<std::uint32_t>(queries.shape(0));
    SearchParams params = make_search_params(
        k, itopk_size, search_width, min_iterations, max_iterations, hash_bitlen, num_threads);
    const std::uint32_t result_k = k < index.size() ? k : index.size();
    py::array_t<std::int32_t> labels({static_cast<py::ssize_t>(nq), static_cast<py::ssize_t>(result_k)});
    py::array_t<float> distances({static_cast<py::ssize_t>(nq), static_cast<py::ssize_t>(result_k)});
    {
        py::gil_scoped_release release;
        if (cpu) {
            index.search_cpu(queries.data(), nq, params, labels.mutable_data(), distances.mutable_data());
        } else {
            index.search(queries.data(), nq, params, labels.mutable_data(), distances.mutable_data());
        }
    }
    return {std::move(labels), std::move(distances)};
}

std::pair<py::array_t<std::int32_t>, py::array_t<float>> search_index(const Index& index,
                                                                      const py::array& queries_in,
                                                                      std::uint32_t k,
                                                                      std::uint32_t itopk_size,
                                                                      std::uint32_t search_width,
                                                                      std::uint32_t min_iterations,
                                                                      std::uint32_t max_iterations,
                                                                      std::uint32_t hash_bitlen,
                                                                      std::uint32_t num_threads) {
    return run_search(index, queries_in, k, itopk_size, search_width, min_iterations,
                      max_iterations, hash_bitlen, num_threads, false);
}

std::pair<py::array_t<std::int32_t>, py::array_t<float>> search_cpu_index(
    const Index& index,
    const py::array& queries_in,
    std::uint32_t k,
    std::uint32_t itopk_size,
    std::uint32_t search_width,
    std::uint32_t min_iterations,
    std::uint32_t max_iterations,
    std::uint32_t hash_bitlen,
    std::uint32_t num_threads) {
    return run_search(index, queries_in, k, itopk_size, search_width, min_iterations,
                      max_iterations, hash_bitlen, num_threads, true);
}

}  // namespace

PYBIND11_MODULE(pq_cagra_adc_ext, m) {
    m.doc() = "CAGRA search over node-level PQ windows with shared-memory ADC";
    py::class_<Index>(m, "Index")
        .def(py::init<>())
        .def("build",
             &build_index,
             py::arg("codes"),
             py::arg("codebook"),
             py::arg("window_size") = 0,
             py::arg("pq_m") = 0,
             py::arg("nbits") = 0,
             py::arg("dsub") = 0,
             py::arg("graph_degree") = 32,
             py::arg("intermediate_graph_degree") = 64,
             py::arg("nndescent_iterations") = 0)
        .def("search",
             &search_index,
             py::arg("queries"),
             py::arg("k") = 10,
             py::arg("itopk_size") = 64,
             py::arg("search_width") = 1,
             py::arg("min_iterations") = 8,
             py::arg("max_iterations") = 64,
             py::arg("hash_bitlen") = 16,
             py::arg("num_threads") = 0)
        .def("search_cpu",
             &search_cpu_index,
             py::arg("queries"),
             py::arg("k") = 10,
             py::arg("itopk_size") = 64,
             py::arg("search_width") = 1,
             py::arg("min_iterations") = 8,
             py::arg("max_iterations") = 64,
             py::arg("hash_bitlen") = 16,
             py::arg("num_threads") = 0)
        .def("save", &Index::save, py::arg("path"))
        .def_static("load", &Index::load, py::arg("path"))
        .def_property_readonly("size", &Index::size)
        .def_property_readonly("window_size", &Index::window_size)
        .def_property_readonly("pq_m", &Index::pq_m)
        .def_property_readonly("nbits", &Index::nbits)
        .def_property_readonly("ksub", &Index::ksub)
        .def_property_readonly("dsub", &Index::dsub)
        .def_property_readonly("dim", &Index::dim)
        .def_property_readonly("code_dim", &Index::code_dim)
        .def_property_readonly("graph_degree", &Index::graph_degree);
}
