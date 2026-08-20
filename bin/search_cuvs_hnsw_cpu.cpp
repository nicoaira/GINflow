// Search a cuVS HNSW file from host memory without loading the CAGRA graph.
//
// The file is produced by cagra_to_hnsw_cpu.cpp. cuVS's HNSW wrapper is not
// compatible with the original hnswlib serialization format, but its search
// path is CPU-side and supports int8 vectors.

#include <cuda_runtime_api.h>

#include <cuvs/core/c_api.h>
#include <cuvs/distance/distance.h>
#include <cuvs/neighbors/hnsw.h>
#include <dlpack/dlpack.h>

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;

void check_cuvs(cuvsError_t status, const char* expression)
{
  if (status == CUVS_SUCCESS) { return; }
  const auto* detail = cuvsGetLastErrorText();
  throw std::runtime_error(
    std::string(expression) + " failed" + (detail == nullptr ? "" : ": " + std::string(detail)));
}

#define CHECK_CUVS(expression) check_cuvs((expression), #expression)

std::vector<int> parse_ints(const std::string& value)
{
  std::vector<int> values;
  std::stringstream stream(value);
  std::string item;
  while (std::getline(stream, item, ',')) {
    if (item.empty()) { continue; }
    const auto parsed = std::stoi(item);
    if (parsed < 1) { throw std::invalid_argument("integer values must be positive"); }
    values.push_back(parsed);
  }
  if (values.empty()) { throw std::invalid_argument("integer list is empty"); }
  return values;
}

template <typename T>
void read_binary(const std::string& path, std::vector<T>& values)
{
  std::ifstream input(path, std::ios::binary);
  if (!input) { throw std::runtime_error("cannot open query file: " + path); }
  input.seekg(0, std::ios::end);
  const auto bytes = input.tellg();
  input.seekg(0, std::ios::beg);
  const auto expected = static_cast<std::streamoff>(values.size() * sizeof(T));
  if (bytes != expected) { throw std::runtime_error("query file has unexpected size: " + path); }
  input.read(reinterpret_cast<char*>(values.data()), expected);
  if (!input) { throw std::runtime_error("failed to read query file: " + path); }
}

void write_binary(const std::string& path, const std::vector<std::uint64_t>& values)
{
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) { throw std::runtime_error("cannot open output file: " + path); }
  output.write(reinterpret_cast<const char*>(values.data()),
               static_cast<std::streamsize>(values.size() * sizeof(std::uint64_t)));
}

void write_binary(const std::string& path, const std::vector<float>& values)
{
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) { throw std::runtime_error("cannot open output file: " + path); }
  output.write(reinterpret_cast<const char*>(values.data()),
               static_cast<std::streamsize>(values.size() * sizeof(float)));
}

struct Tensor2D {
  DLManagedTensor tensor{};
  std::int64_t shape[2]{};

  Tensor2D(void* data, std::int64_t rows, std::int64_t columns, DLDataType dtype)
  {
    shape[0] = rows;
    shape[1] = columns;
    tensor.dl_tensor.data = data;
    tensor.dl_tensor.device = {kDLCPU, 0};
    tensor.dl_tensor.ndim = 2;
    tensor.dl_tensor.dtype = dtype;
    tensor.dl_tensor.shape = shape;
    tensor.dl_tensor.strides = nullptr;
    tensor.dl_tensor.byte_offset = 0;
  }
};

DLDataType dtype(std::uint8_t code, std::uint8_t bits)
{
  return {code, bits, 1};
}

std::string option(int argc, char** argv, const std::string& name)
{
  for (int i = 1; i + 1 < argc; ++i) {
    if (name == argv[i]) { return argv[i + 1]; }
  }
  throw std::invalid_argument("missing option: " + name);
}

void usage(const char* program)
{
  std::cerr << "usage: " << program
            << " --hnsw-index PATH --queries-int8 PATH --n-queries N"
               " --dimension D --k K --ef-values LIST --num-threads N"
               " --labels-prefix PATH --distances-prefix PATH --metrics PATH\n";
}

}  // namespace

int main(int argc, char** argv)
{
  try {
    if (argc < 19) {
      usage(argv[0]);
      return 2;
    }
    const auto hnsw_path        = option(argc, argv, "--hnsw-index");
    const auto query_path       = option(argc, argv, "--queries-int8");
    const auto n_queries        = static_cast<std::size_t>(std::stoull(option(argc, argv, "--n-queries")));
    const auto dimension        = std::stoi(option(argc, argv, "--dimension"));
    const auto k                = std::stoi(option(argc, argv, "--k"));
    const auto efs              = parse_ints(option(argc, argv, "--ef-values"));
    const auto threads          = std::stoi(option(argc, argv, "--num-threads"));
    const auto labels_prefix    = option(argc, argv, "--labels-prefix");
    const auto distances_prefix = option(argc, argv, "--distances-prefix");
    const auto metrics_path     = option(argc, argv, "--metrics");
    if (dimension < 1 || k < 1 || threads < 0) {
      throw std::invalid_argument("invalid dimension, k, or thread count");
    }
    for (const auto ef : efs) {
      if (ef < k) { throw std::invalid_argument("every ef value must be >= k"); }
    }

    std::vector<std::int8_t> queries(n_queries * static_cast<std::size_t>(dimension));
    read_binary(query_path, queries);
    Tensor2D query_tensor(queries.data(), n_queries, dimension, dtype(kDLInt, 8));

    cuvsResources_t resources = 0;
    cuvsHnswIndex_t index = nullptr;
    cuvsHnswSearchParams_t params = nullptr;
    CHECK_CUVS(cuvsResourcesCreate(&resources));
    CHECK_CUVS(cuvsHnswIndexCreate(&index));
    index->dtype = dtype(kDLInt, 8);

    const auto load_started = clock_type::now();
    CHECK_CUVS(cuvsHnswDeserialize(resources, hnsw_path.c_str(), dimension, L2Expanded, index));
    const auto load_seconds =
      std::chrono::duration<double>(clock_type::now() - load_started).count();
    CHECK_CUVS(cuvsHnswSearchParamsCreate(&params));

    std::vector<double> search_seconds;
    search_seconds.reserve(efs.size());
    for (const auto ef : efs) {
      std::vector<std::uint64_t> neighbors(n_queries * static_cast<std::size_t>(k));
      std::vector<float> distances(n_queries * static_cast<std::size_t>(k));
      Tensor2D neighbor_tensor(neighbors.data(), n_queries, k, dtype(kDLUInt, 64));
      Tensor2D distance_tensor(distances.data(), n_queries, k, dtype(kDLFloat, 32));
      params->ef = ef;
      params->numThreads = threads;
      CHECK_CUVS(cuvsHnswSearch(
        resources,
        params,
        index,
        &query_tensor.tensor,
        &neighbor_tensor.tensor,
        &distance_tensor.tensor));
      const auto started = clock_type::now();
      CHECK_CUVS(cuvsHnswSearch(
        resources,
        params,
        index,
        &query_tensor.tensor,
        &neighbor_tensor.tensor,
        &distance_tensor.tensor));
      search_seconds.push_back(
        std::chrono::duration<double>(clock_type::now() - started).count());
      write_binary(labels_prefix + std::to_string(ef) + ".u64", neighbors);
      write_binary(distances_prefix + std::to_string(ef) + ".f32", distances);
    }

    std::ofstream metrics(metrics_path, std::ios::trunc);
    if (!metrics) { throw std::runtime_error("cannot open metrics file: " + metrics_path); }
    metrics << "{\n";
    metrics << "  \"backend\": \"cuvs_hnsw_cpu_search_only\",\n";
    metrics << "  \"hnsw_index\": \"" << hnsw_path << "\",\n";
    metrics << "  \"n_queries\": " << n_queries << ",\n";
    metrics << "  \"dimension\": " << dimension << ",\n";
    metrics << "  \"k\": " << k << ",\n";
    metrics << "  \"num_threads\": " << threads << ",\n";
    metrics << "  \"hnsw_load_seconds\": " << load_seconds << ",\n";
    metrics << "  \"results\": [\n";
    for (std::size_t i = 0; i < efs.size(); ++i) {
      metrics << "    {\"ef_search\": " << efs[i]
              << ", \"search_seconds\": " << search_seconds[i]
              << ", \"labels_file\": \"" << labels_prefix << efs[i]
              << ".u64\", \"distances_file\": \"" << distances_prefix << efs[i]
              << ".f32\"}" << (i + 1 == efs.size() ? "\n" : ",\n");
    }
    metrics << "  ]\n}\n";

    CHECK_CUVS(cuvsHnswSearchParamsDestroy(params));
    CHECK_CUVS(cuvsHnswIndexDestroy(index));
    CHECK_CUVS(cuvsResourcesDestroy(resources));
    std::cout << "hnsw_load_seconds=" << load_seconds;
    for (std::size_t i = 0; i < efs.size(); ++i) {
      std::cout << " ef" << efs[i] << "_search_seconds=" << search_seconds[i];
    }
    std::cout << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }
}
