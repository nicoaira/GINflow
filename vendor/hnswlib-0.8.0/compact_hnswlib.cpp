// Compact hnswlib driver for fixed-width uint16 centroid-code windows.
// The bundled hnswlib headers are pinned to v0.8.0.
#include "hnswlib/hnswlib.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Options {
    std::string command;
    std::string codes;
    std::string similarity;
    std::string index;
    std::string labels_out;
    std::string distances_out;
    std::size_t count = 0;
    std::size_t query_count = 0;
    std::size_t window_size = 0;
    std::size_t n_centroids = 0;
    std::size_t k = 0;
    std::size_t m = 32;
    std::size_t ef_construction = 200;
    std::size_t ef_search = 100;
    std::size_t threads = 0;
    std::size_t random_seed = 1;
};

std::string require_value(const std::string &key, const std::string &value) {
    if (value.empty()) {
        throw std::runtime_error("missing value for " + key);
    }
    return value;
}

std::size_t parse_size(const std::string &key, const std::string &value) {
    const std::string checked = require_value(key, value);
    std::size_t consumed = 0;
    unsigned long long parsed = 0;
    try {
        parsed = std::stoull(checked, &consumed);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid integer for " + key + ": " + checked);
    }
    if (consumed != checked.size()) {
        throw std::runtime_error("invalid integer for " + key + ": " + checked);
    }
    if (parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("integer is too large for " + key);
    }
    return static_cast<std::size_t>(parsed);
}

Options parse_options(int argc, char **argv) {
    if (argc < 2) {
        throw std::runtime_error("usage: compact_hnswlib build|search [options]");
    }
    Options options;
    options.command = argv[1];
    for (int i = 2; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument.rfind("--", 0) != 0) {
            throw std::runtime_error("unexpected argument: " + argument);
        }
        const std::size_t equals = argument.find('=');
        std::string key = equals == std::string::npos ? argument : argument.substr(0, equals);
        std::string value;
        if (equals == std::string::npos) {
            if (i + 1 >= argc) {
                throw std::runtime_error("missing value for " + key);
            }
            value = argv[++i];
        } else {
            value = argument.substr(equals + 1);
        }
        if (key == "--codes") options.codes = require_value(key, value);
        else if (key == "--similarity") options.similarity = require_value(key, value);
        else if (key == "--index") options.index = require_value(key, value);
        else if (key == "--labels-out") options.labels_out = require_value(key, value);
        else if (key == "--distances-out") options.distances_out = require_value(key, value);
        else if (key == "--count") options.count = parse_size(key, value);
        else if (key == "--query-count") options.query_count = parse_size(key, value);
        else if (key == "--window-size") options.window_size = parse_size(key, value);
        else if (key == "--n-centroids") options.n_centroids = parse_size(key, value);
        else if (key == "--k") options.k = parse_size(key, value);
        else if (key == "--m") options.m = parse_size(key, value);
        else if (key == "--ef-construction") options.ef_construction = parse_size(key, value);
        else if (key == "--ef-search") options.ef_search = parse_size(key, value);
        else if (key == "--threads") options.threads = parse_size(key, value);
        else if (key == "--random-seed") options.random_seed = parse_size(key, value);
        else throw std::runtime_error("unknown option: " + key);
    }
    return options;
}

void require_common(const Options &options) {
    if (options.window_size < 1 || options.n_centroids < 1) {
        throw std::runtime_error("window size and centroid count must be positive");
    }
    if (options.n_centroids > std::numeric_limits<std::uint16_t>::max() + 1ULL) {
        throw std::runtime_error("compact code windows support at most 65536 centroids");
    }
    if (options.m < 2 || options.ef_construction < 1 || options.ef_search < 1) {
        throw std::runtime_error("M, ef-construction, and ef-search must be positive (M >= 2)");
    }
}

std::vector<std::uint16_t> read_codes(const std::string &path, std::size_t count, std::size_t window_size) {
    if (count < 1) {
        throw std::runtime_error("code count must be positive");
    }
    const std::uint64_t expected = static_cast<std::uint64_t>(count) * window_size * sizeof(std::uint16_t);
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("cannot open code file: " + path);
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::uint64_t>(size) != expected) {
        throw std::runtime_error("code file has unexpected size: " + path);
    }
    input.seekg(0, std::ios::beg);
    std::vector<std::uint16_t> codes(static_cast<std::size_t>(count * window_size));
    input.read(reinterpret_cast<char *>(codes.data()), static_cast<std::streamsize>(expected));
    if (!input) throw std::runtime_error("failed reading code file: " + path);
    return codes;
}

std::vector<float> read_similarity(const std::string &path, std::size_t n_centroids) {
    const std::uint64_t expected = static_cast<std::uint64_t>(n_centroids) * n_centroids * sizeof(float);
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("cannot open similarity file: " + path);
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::uint64_t>(size) != expected) {
        throw std::runtime_error("similarity file has unexpected size: " + path);
    }
    input.seekg(0, std::ios::beg);
    std::vector<float> similarity(static_cast<std::size_t>(n_centroids * n_centroids));
    input.read(reinterpret_cast<char *>(similarity.data()), static_cast<std::streamsize>(expected));
    if (!input) throw std::runtime_error("failed reading similarity file: " + path);
    return similarity;
}

struct CodeDistanceParams {
    const float *similarity;
    std::size_t n_centroids;
    std::size_t window_size;
};

float code_distance(const void *left, const void *right, const void *raw_params) {
    const auto *params = static_cast<const CodeDistanceParams *>(raw_params);
    const auto *left_codes = static_cast<const std::uint16_t *>(left);
    const auto *right_codes = static_cast<const std::uint16_t *>(right);
    float score = 0.0f;
    for (std::size_t position = 0; position < params->window_size; ++position) {
        const std::size_t left_code = left_codes[position];
        const std::size_t right_code = right_codes[position];
        if (left_code >= params->n_centroids || right_code >= params->n_centroids) {
            return std::numeric_limits<float>::infinity();
        }
        score += params->similarity[left_code * params->n_centroids + right_code];
    }
    // hnswlib ranks smaller distances first; negating preserves the desired score order.
    return -score;
}

class CodeSpace final : public hnswlib::SpaceInterface<float> {
public:
    CodeSpace(std::size_t n_centroids, std::size_t window_size, const float *similarity)
        : params_{similarity, n_centroids, window_size}, data_size_(window_size * sizeof(std::uint16_t)) {}

    std::size_t get_data_size() override { return data_size_; }
    hnswlib::DISTFUNC<float> get_dist_func() override { return code_distance; }
    void *get_dist_func_param() override { return &params_; }

private:
    CodeDistanceParams params_;
    std::size_t data_size_;
};

void configure_threads(std::size_t threads) {
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(static_cast<int>(threads));
#else
    (void)threads;
#endif
}

void build_index(const Options &options) {
    require_common(options);
    if (options.count < 1 || options.index.empty() || options.codes.empty() || options.similarity.empty()) {
        throw std::runtime_error("build requires --codes, --similarity, --count, and --index");
    }
    const std::vector<std::uint16_t> codes = read_codes(options.codes, options.count, options.window_size);
    const std::vector<float> similarity = read_similarity(options.similarity, options.n_centroids);
    CodeSpace space(options.n_centroids, options.window_size, similarity.data());
    hnswlib::HierarchicalNSW<float> index(
        &space,
        options.count,
        options.m,
        options.ef_construction,
        options.random_seed
    );
    index.setEf(options.ef_search);
    configure_threads(options.threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) if (options.threads > 1)
#endif
    for (std::int64_t row = 0; row < static_cast<std::int64_t>(options.count); ++row) {
        index.addPoint(codes.data() + static_cast<std::size_t>(row) * options.window_size,
                       static_cast<hnswlib::labeltype>(row));
    }
    index.saveIndex(options.index);
    std::cout << "built compact hnswlib index: elements=" << options.count
              << " code_bytes=" << options.window_size * sizeof(std::uint16_t)
              << " threads=" << options.threads << std::endl;
}

void write_results(
    const std::string &labels_path,
    const std::string &distances_path,
    const std::vector<std::uint64_t> &labels,
    const std::vector<float> &distances
) {
    std::ofstream labels_out(labels_path, std::ios::binary);
    std::ofstream distances_out(distances_path, std::ios::binary);
    if (!labels_out || !distances_out) throw std::runtime_error("cannot create search result files");
    labels_out.write(reinterpret_cast<const char *>(labels.data()),
                     static_cast<std::streamsize>(labels.size() * sizeof(std::uint64_t)));
    distances_out.write(reinterpret_cast<const char *>(distances.data()),
                        static_cast<std::streamsize>(distances.size() * sizeof(float)));
    if (!labels_out || !distances_out) throw std::runtime_error("failed writing search result files");
}

void search_index(const Options &options) {
    require_common(options);
    if (options.query_count < 1 || options.k < 1 || options.index.empty() || options.codes.empty()
        || options.similarity.empty() || options.labels_out.empty() || options.distances_out.empty()) {
        throw std::runtime_error("search requires code/index/similarity paths, query count, k, and output paths");
    }
    const std::vector<std::uint16_t> codes = read_codes(options.codes, options.query_count, options.window_size);
    const std::vector<float> similarity = read_similarity(options.similarity, options.n_centroids);
    CodeSpace space(options.n_centroids, options.window_size, similarity.data());
    hnswlib::HierarchicalNSW<float> index(&space, options.index);
    index.setEf(std::max(options.ef_search, options.k));
    configure_threads(options.threads);
    const std::size_t available = index.getCurrentElementCount();
    const std::size_t result_count = std::min(options.k, available);
    std::vector<std::uint64_t> labels(options.query_count * options.k, std::numeric_limits<std::uint64_t>::max());
    std::vector<float> distances(options.query_count * options.k, std::numeric_limits<float>::infinity());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) if (options.threads > 1)
#endif
    for (std::int64_t row = 0; row < static_cast<std::int64_t>(options.query_count); ++row) {
        const auto result = index.searchKnnCloserFirst(
            codes.data() + static_cast<std::size_t>(row) * options.window_size,
            result_count
        );
        for (std::size_t rank = 0; rank < result.size(); ++rank) {
            labels[static_cast<std::size_t>(row) * options.k + rank] = result[rank].second;
            distances[static_cast<std::size_t>(row) * options.k + rank] = result[rank].first;
        }
    }
    write_results(options.labels_out, options.distances_out, labels, distances);
    std::cout << "searched compact hnswlib index: queries=" << options.query_count
              << " k=" << options.k << " threads=" << options.threads << std::endl;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.command == "build") build_index(options);
        else if (options.command == "search") search_index(options);
        else throw std::runtime_error("command must be build or search");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << std::endl;
        return 1;
    }
}
