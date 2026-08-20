// Benchmark cuVS CAGRA-to-CPU-HNSW conversion and host-side search.
//
// This is intentionally a standalone research driver. It uses the cuVS C
// API because cuVS 24.10's Python hnsw wrapper does not accept the int8 dtype
// used by GINflow's GPU path, while the C API supports int8 HNSW search.

#include <cuda_runtime_api.h>

#include <cuvs/core/c_api.h>
#include <cuvs/distance/distance.h>
#include <cuvs/neighbors/cagra.h>
#include <cuvs/neighbors/hnsw.h>
#include <dlpack/dlpack.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;

struct Timings {
  double cagra_load_seconds = 0.0;
  double cagra_to_hnsw_seconds = 0.0;
  double hnsw_load_seconds = 0.0;
  std::vector<double> search_seconds;
  std::size_t baseline_vram = 0;
  std::size_t peak_vram = 0;
  std::size_t post_cagra_destroy_vram = 0;
};

void check_cuvs(cuvsError_t status, const char* expression)
{
  if (status == CUVS_SUCCESS) { return; }
  const auto* detail = cuvsGetLastErrorText();
  throw std::runtime_error(
    std::string(expression) + " failed" + (detail == nullptr ? "" : ": " + std::string(detail)));
}

#define CHECK_CUVS(expression) check_cuvs((expression), #expression)

std::size_t used_vram()
{
  std::size_t free_bytes = 0;
  std::size_t total_bytes = 0;
  if (cudaMemGetInfo(&free_bytes, &total_bytes) != cudaSuccess) { return 0; }
  return total_bytes - free_bytes;
}

class VramSampler {
 public:
  explicit VramSampler(std::size_t interval_millis = 10) : interval_millis_(interval_millis) {}
  ~VramSampler() { stop(); }

  void start()
  {
    baseline_ = used_vram();
    peak_.store(baseline_.load());
    running_ = true;
    thread_  = std::thread([this]() {
      while (running_) {
        peak_.store(std::max(peak_.load(), used_vram()));
        std::this_thread::sleep_for(std::chrono::milliseconds(interval_millis_));
      }
      peak_.store(std::max(peak_.load(), used_vram()));
    });
  }

  void stop()
  {
    running_ = false;
    if (thread_.joinable()) { thread_.join(); }
  }

  std::size_t baseline() const { return baseline_; }
  std::size_t peak() const { return peak_; }

 private:
  std::size_t interval_millis_;
  std::atomic<bool> running_{false};
  std::atomic<std::size_t> baseline_{0};
  std::atomic<std::size_t> peak_{0};
  std::thread thread_;
};

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
  if (bytes != expected) {
    throw std::runtime_error("query file has unexpected size: " + path);
  }
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
    tensor.manager_ctx = nullptr;
    tensor.deleter = nullptr;
  }
};

DLDataType int8_dtype()
{
  return {kDLInt, 8, 1};
}

DLDataType uint64_dtype()
{
  return {kDLUInt, 64, 1};
}

DLDataType float32_dtype()
{
  return {kDLFloat, 32, 1};
}

void write_json(
  const std::string& path,
  const Timings& timings,
  const std::vector<int>& efs,
  std::size_t n_queries,
  int dimension,
  int k,
  int threads,
  const std::string& hnsw_path,
  const std::string& labels_prefix,
  const std::string& distances_prefix)
{
  std::ofstream output(path, std::ios::trunc);
  if (!output) { throw std::runtime_error("cannot open metrics file: " + path); }
  output << "{\n";
  output << "  \"backend\": \"cuvs_c_api_cagra_to_hnsw_cpu\",\n";
  output << "  \"n_queries\": " << n_queries << ",\n";
  output << "  \"dimension\": " << dimension << ",\n";
  output << "  \"k\": " << k << ",\n";
  output << "  \"num_threads\": " << threads << ",\n";
  output << "  \"hnsw_index\": \"" << hnsw_path << "\",\n";
  output << "  \"cagra_load_seconds\": " << timings.cagra_load_seconds << ",\n";
  output << "  \"cagra_to_hnsw_conversion_seconds\": "
          << timings.cagra_to_hnsw_seconds << ",\n";
  output << "  \"hnsw_load_seconds\": " << timings.hnsw_load_seconds << ",\n";
  output << "  \"baseline_vram_bytes\": " << timings.baseline_vram << ",\n";
  output << "  \"peak_vram_bytes\": " << timings.peak_vram << ",\n";
  output << "  \"peak_vram_delta_bytes\": "
          << (timings.peak_vram >= timings.baseline_vram
                ? timings.peak_vram - timings.baseline_vram
                : 0)
          << ",\n";
  output << "  \"post_cagra_destroy_vram_bytes\": "
          << timings.post_cagra_destroy_vram << ",\n";
  output << "  \"results\": [\n";
  for (std::size_t i = 0; i < efs.size(); ++i) {
    output << "    {\n";
    output << "      \"ef_search\": " << efs[i] << ",\n";
    output << "      \"search_seconds\": " << timings.search_seconds[i] << ",\n";
    output << "      \"labels_file\": \"" << labels_prefix << efs[i] << ".u64\",\n";
    output << "      \"distances_file\": \"" << distances_prefix << efs[i]
            << ".f32\"\n";
    output << "    }" << (i + 1 == efs.size() ? "\n" : ",\n");
  }
  output << "  ]\n";
  output << "}\n";
}

void usage(const char* program)
{
  std::cerr << "usage: " << program
            << " --cagra-index PATH --hnsw-index PATH --queries-int8 PATH"
               " --n-queries N --dimension D --k K --ef-values LIST"
               " --num-threads N --labels-prefix PATH --distances-prefix PATH"
               " --metrics PATH\n";
}

std::string option(int argc, char** argv, const std::string& name)
{
  for (int i = 1; i + 1 < argc; ++i) {
    if (name == argv[i]) { return argv[i + 1]; }
  }
  throw std::invalid_argument("missing option: " + name);
}

}  // namespace

int main(int argc, char** argv)
{
  try {
    if (argc < 21) {
      usage(argv[0]);
      return 2;
    }
    const auto cagra_path       = option(argc, argv, "--cagra-index");
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
    Tensor2D query_tensor(queries.data(), n_queries, dimension, int8_dtype());

    cuvsResources_t resources = 0;
    cuvsCagraIndex_t cagra_index = nullptr;
    cuvsHnswIndex_t hnsw_index = nullptr;
    cuvsHnswSearchParams_t search_params = nullptr;
    CHECK_CUVS(cuvsResourcesCreate(&resources));

    Timings timings;
    VramSampler sampler;
    sampler.start();

    CHECK_CUVS(cuvsCagraIndexCreate(&cagra_index));
    auto started = clock_type::now();
    CHECK_CUVS(cuvsCagraDeserialize(resources, cagra_path.c_str(), cagra_index));
    CHECK_CUVS(cuvsStreamSync(resources));
    timings.cagra_load_seconds =
      std::chrono::duration<double>(clock_type::now() - started).count();
    const auto dtype = cagra_index->dtype;
    int loaded_dimension = 0;
    CHECK_CUVS(cuvsCagraIndexGetDims(cagra_index, &loaded_dimension));
    if (loaded_dimension != dimension) {
      throw std::runtime_error("CAGRA index dimension does not match query dimension");
    }

    started = clock_type::now();
    CHECK_CUVS(cuvsCagraSerializeToHnswlib(resources, hnsw_path.c_str(), cagra_index));
    CHECK_CUVS(cuvsStreamSync(resources));
    timings.cagra_to_hnsw_seconds =
      std::chrono::duration<double>(clock_type::now() - started).count();
    sampler.stop();
    timings.baseline_vram = sampler.baseline();
    timings.peak_vram = sampler.peak();

    CHECK_CUVS(cuvsCagraIndexDestroy(cagra_index));
    cagra_index = nullptr;
    CHECK_CUVS(cuvsStreamSync(resources));
    timings.post_cagra_destroy_vram = used_vram();

    CHECK_CUVS(cuvsHnswIndexCreate(&hnsw_index));
    hnsw_index->dtype = dtype;
    started = clock_type::now();
    CHECK_CUVS(cuvsHnswDeserialize(
      resources, hnsw_path.c_str(), dimension, L2Expanded, hnsw_index));
    timings.hnsw_load_seconds =
      std::chrono::duration<double>(clock_type::now() - started).count();
    CHECK_CUVS(cuvsHnswSearchParamsCreate(&search_params));

    timings.search_seconds.reserve(efs.size());
    for (const auto ef : efs) {
      std::vector<std::uint64_t> neighbors(n_queries * static_cast<std::size_t>(k));
      std::vector<float> distances(n_queries * static_cast<std::size_t>(k));
      Tensor2D neighbor_tensor(neighbors.data(), n_queries, k, uint64_dtype());
      Tensor2D distance_tensor(distances.data(), n_queries, k, float32_dtype());
      search_params->ef = ef;
      search_params->numThreads = threads;

      CHECK_CUVS(cuvsHnswSearch(
        resources,
        search_params,
        hnsw_index,
        &query_tensor.tensor,
        &neighbor_tensor.tensor,
        &distance_tensor.tensor));
      auto search_started = clock_type::now();
      CHECK_CUVS(cuvsHnswSearch(
        resources,
        search_params,
        hnsw_index,
        &query_tensor.tensor,
        &neighbor_tensor.tensor,
        &distance_tensor.tensor));
      timings.search_seconds.push_back(
        std::chrono::duration<double>(clock_type::now() - search_started).count());
      write_binary(labels_prefix + std::to_string(ef) + ".u64", neighbors);
      write_binary(distances_prefix + std::to_string(ef) + ".f32", distances);
    }

    write_json(
      metrics_path,
      timings,
      efs,
      n_queries,
      dimension,
      k,
      threads,
      hnsw_path,
      labels_prefix,
      distances_prefix);

    CHECK_CUVS(cuvsHnswSearchParamsDestroy(search_params));
    search_params = nullptr;
    CHECK_CUVS(cuvsHnswIndexDestroy(hnsw_index));
    hnsw_index = nullptr;
    CHECK_CUVS(cuvsResourcesDestroy(resources));
    resources = 0;
    std::cout << "cagra_load_seconds=" << timings.cagra_load_seconds
              << " cagra_to_hnsw_seconds=" << timings.cagra_to_hnsw_seconds
              << " hnsw_load_seconds=" << timings.hnsw_load_seconds
              << " peak_vram_delta_bytes="
              << (timings.peak_vram >= timings.baseline_vram
                    ? timings.peak_vram - timings.baseline_vram
                    : 0)
              << " post_cagra_destroy_vram_bytes=" << timings.post_cagra_destroy_vram << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }
}
