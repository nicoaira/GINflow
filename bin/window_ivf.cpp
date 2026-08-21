// Standalone IVF prototype for node-level PQ windows.
//
// Every window is a sequence of packed node-PQ codes.  Coarse assignment,
// nprobe routing, and list scanning all use the registered positional lookup
// score: sum(position, subquantizer) similarity[query_code, target_code].
// The original float16 residue embeddings remain outside this candidate index
// and can be used by the exact reranker after search.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
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
    std::size_t pq_m = 0;
    std::size_t nbits = 0;
    std::size_t k = 0;
    std::size_t nlist = 256;
    std::size_t nprobe = 8;
    std::size_t niter = 1;
    std::size_t threads = 0;
    std::size_t random_seed = 1;
};

std::string require_value(const std::string &key, const std::string &value) {
    if (value.empty()) throw std::runtime_error("missing value for " + key);
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
    if (consumed != checked.size() || parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("invalid integer for " + key + ": " + checked);
    }
    return static_cast<std::size_t>(parsed);
}

Options parse_options(int argc, char **argv) {
    if (argc < 2) throw std::runtime_error("usage: window_ivf build|search [options]");
    Options options;
    options.command = argv[1];
    for (int i = 2; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument.rfind("--", 0) != 0) throw std::runtime_error("unexpected argument: " + argument);
        const std::size_t equals = argument.find('=');
        const std::string key = equals == std::string::npos ? argument : argument.substr(0, equals);
        std::string value;
        if (equals == std::string::npos) {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + key);
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
        else if (key == "--pq-m") options.pq_m = parse_size(key, value);
        else if (key == "--nbits") options.nbits = parse_size(key, value);
        else if (key == "--k") options.k = parse_size(key, value);
        else if (key == "--nlist") options.nlist = parse_size(key, value);
        else if (key == "--nprobe") options.nprobe = parse_size(key, value);
        else if (key == "--niter") options.niter = parse_size(key, value);
        else if (key == "--threads") options.threads = parse_size(key, value);
        else if (key == "--random-seed") options.random_seed = parse_size(key, value);
        else throw std::runtime_error("unknown option: " + key);
    }
    return options;
}

void require_common(const Options &options) {
    if (options.window_size < 1 || options.pq_m < 1 || options.nbits < 1 || options.nbits > 8) {
        throw std::runtime_error("window size, PQ subquantizers, and nbits are invalid");
    }
    if (options.nlist < 1 || options.nprobe < 1 || options.niter < 1) {
        throw std::runtime_error("nlist, nprobe, and niter must be positive");
    }
}

std::size_t ksub(const Options &options) {
    return static_cast<std::size_t>(1) << options.nbits;
}

std::size_t node_code_bytes(const Options &options) {
    return (options.pq_m * options.nbits + 7) / 8;
}

std::size_t code_bytes(const Options &options) {
    if (node_code_bytes(options) > std::numeric_limits<std::size_t>::max() / options.window_size) {
        throw std::runtime_error("code size overflows size_t");
    }
    return options.window_size * node_code_bytes(options);
}

std::vector<std::uint8_t> read_codes(const std::string &path, std::size_t count, std::size_t bytes_per_window) {
    if (count < 1) throw std::runtime_error("code count must be positive");
    const std::uint64_t expected = static_cast<std::uint64_t>(count) * bytes_per_window;
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("cannot open code file: " + path);
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::uint64_t>(size) != expected) {
        throw std::runtime_error("code file has unexpected size: " + path);
    }
    input.seekg(0, std::ios::beg);
    std::vector<std::uint8_t> codes(static_cast<std::size_t>(expected));
    input.read(reinterpret_cast<char *>(codes.data()), static_cast<std::streamsize>(expected));
    if (!input) throw std::runtime_error("failed reading code file: " + path);
    return codes;
}

std::vector<float> read_similarity(const std::string &path, std::size_t pq_m, std::size_t ksub_value) {
    const std::uint64_t table_size = static_cast<std::uint64_t>(pq_m) * ksub_value * ksub_value;
    const std::uint64_t expected = table_size * sizeof(float);
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("cannot open similarity file: " + path);
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::uint64_t>(size) != expected) {
        throw std::runtime_error("similarity file has unexpected size: " + path);
    }
    input.seekg(0, std::ios::beg);
    std::vector<float> similarity(static_cast<std::size_t>(table_size));
    input.read(reinterpret_cast<char *>(similarity.data()), static_cast<std::streamsize>(expected));
    if (!input) throw std::runtime_error("failed reading similarity file: " + path);
    return similarity;
}

std::size_t unpack_code(const std::uint8_t *node_codes, std::size_t subquantizer, std::size_t nbits) {
    const std::size_t bit_offset = subquantizer * nbits;
    const std::size_t byte_offset = bit_offset / 8;
    const std::size_t shift = bit_offset % 8;
    std::uint16_t packed = node_codes[byte_offset];
    if (shift + nbits > 8) packed |= static_cast<std::uint16_t>(node_codes[byte_offset + 1]) << 8;
    const std::uint16_t mask = static_cast<std::uint16_t>((static_cast<std::size_t>(1) << nbits) - 1);
    return static_cast<std::size_t>((packed >> shift) & mask);
}

float window_score(
    const std::uint8_t *left,
    const std::uint8_t *right,
    const Options &options,
    const std::vector<float> &similarity
) {
    const std::size_t node_bytes = node_code_bytes(options);
    const std::size_t table_stride = ksub(options) * ksub(options);
    float score = 0.0f;
    for (std::size_t position = 0; position < options.window_size; ++position) {
        const std::uint8_t *left_node = left + position * node_bytes;
        const std::uint8_t *right_node = right + position * node_bytes;
        for (std::size_t subquantizer = 0; subquantizer < options.pq_m; ++subquantizer) {
            const std::size_t left_code = unpack_code(left_node, subquantizer, options.nbits);
            const std::size_t right_code = unpack_code(right_node, subquantizer, options.nbits);
            if (left_code >= ksub(options) || right_code >= ksub(options)) {
                return -std::numeric_limits<float>::infinity();
            }
            const float *table = similarity.data() + subquantizer * table_stride;
            score += table[left_code * ksub(options) + right_code];
        }
    }
    return score;
}

void configure_threads(std::size_t threads) {
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(static_cast<int>(threads));
#else
    (void)threads;
#endif
}

struct IndexHeader {
    char magic[8];
    std::uint64_t count;
    std::uint64_t code_bytes;
    std::uint64_t window_size;
    std::uint64_t pq_m;
    std::uint64_t nbits;
    std::uint64_t nlist;
};

struct IndexData {
    IndexHeader header{};
    std::vector<std::uint8_t> prototypes;
    std::vector<std::uint64_t> offsets;
    std::vector<std::uint8_t> codes;
    std::vector<std::uint64_t> labels;
};

void write_index(const std::string &path, const IndexData &data) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("cannot create IVF index: " + path);
    output.write(reinterpret_cast<const char *>(&data.header), sizeof(data.header));
    output.write(reinterpret_cast<const char *>(data.prototypes.data()), static_cast<std::streamsize>(data.prototypes.size()));
    output.write(reinterpret_cast<const char *>(data.offsets.data()), static_cast<std::streamsize>(data.offsets.size() * sizeof(std::uint64_t)));
    output.write(reinterpret_cast<const char *>(data.codes.data()), static_cast<std::streamsize>(data.codes.size()));
    output.write(reinterpret_cast<const char *>(data.labels.data()), static_cast<std::streamsize>(data.labels.size() * sizeof(std::uint64_t)));
    if (!output) throw std::runtime_error("failed writing IVF index: " + path);
}

IndexData read_index(const std::string &path, const Options &options) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open IVF index: " + path);
    IndexData data;
    input.read(reinterpret_cast<char *>(&data.header), sizeof(data.header));
    if (!input || std::string(data.header.magic, data.header.magic + 8) != "GINIVF01") {
        throw std::runtime_error("invalid WindowIVF index header");
    }
    if (data.header.code_bytes != code_bytes(options) || data.header.window_size != options.window_size ||
        data.header.pq_m != options.pq_m || data.header.nbits != options.nbits) {
        throw std::runtime_error("WindowIVF index parameters do not match the query codes");
    }
    const std::size_t nlist = static_cast<std::size_t>(data.header.nlist);
    data.prototypes.resize(nlist * static_cast<std::size_t>(data.header.code_bytes));
    data.offsets.resize(nlist + 1);
    data.codes.resize(static_cast<std::size_t>(data.header.count) * static_cast<std::size_t>(data.header.code_bytes));
    data.labels.resize(static_cast<std::size_t>(data.header.count));
    input.read(reinterpret_cast<char *>(data.prototypes.data()), static_cast<std::streamsize>(data.prototypes.size()));
    input.read(reinterpret_cast<char *>(data.offsets.data()), static_cast<std::streamsize>(data.offsets.size() * sizeof(std::uint64_t)));
    input.read(reinterpret_cast<char *>(data.codes.data()), static_cast<std::streamsize>(data.codes.size()));
    input.read(reinterpret_cast<char *>(data.labels.data()), static_cast<std::streamsize>(data.labels.size() * sizeof(std::uint64_t)));
    if (!input) throw std::runtime_error("truncated WindowIVF index");
    if (data.offsets.front() != 0 || data.offsets.back() != data.header.count) {
        throw std::runtime_error("invalid WindowIVF list offsets");
    }
    return data;
}

void build_index(const Options &options) {
    require_common(options);
    if (options.count < 1 || options.index.empty() || options.codes.empty() || options.similarity.empty()) {
        throw std::runtime_error("build requires --codes, --similarity, --count, and --index");
    }
    if (options.nlist > options.count) throw std::runtime_error("nlist cannot exceed count");
    const std::size_t bytes = code_bytes(options);
    const std::vector<std::uint8_t> codes = read_codes(options.codes, options.count, bytes);
    const std::vector<float> similarity = read_similarity(options.similarity, options.pq_m, ksub(options));
    configure_threads(options.threads);

    std::vector<std::uint8_t> prototypes(options.nlist * bytes);
    for (std::size_t list = 0; list < options.nlist; ++list) {
        const std::size_t source = (list * options.count) / options.nlist;
        std::copy_n(codes.data() + source * bytes, bytes, prototypes.data() + list * bytes);
    }
    std::vector<std::size_t> assignments(options.count, 0);
    // This is an intentionally simple coarse quantizer: prototypes are
    // deterministic sampled representatives and assignment uses the custom score.
    // It is a transparent baseline for later FAISS InvertedLists integration.
    for (std::size_t iteration = 0; iteration < options.niter; ++iteration) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) if (options.threads > 1)
#endif
        for (std::int64_t row = 0; row < static_cast<std::int64_t>(options.count); ++row) {
            float best_score = -std::numeric_limits<float>::infinity();
            std::size_t best_list = 0;
            const std::uint8_t *query = codes.data() + static_cast<std::size_t>(row) * bytes;
            for (std::size_t list = 0; list < options.nlist; ++list) {
                const float score = window_score(query, prototypes.data() + list * bytes, options, similarity);
                if (score > best_score) {
                    best_score = score;
                    best_list = list;
                }
            }
            assignments[static_cast<std::size_t>(row)] = best_list;
        }
    }

    std::vector<std::uint64_t> counts(options.nlist, 0);
    for (std::size_t assignment : assignments) ++counts[assignment];
    std::vector<std::uint64_t> offsets(options.nlist + 1, 0);
    for (std::size_t list = 0; list < options.nlist; ++list) offsets[list + 1] = offsets[list] + counts[list];
    std::vector<std::uint64_t> cursors = offsets;
    std::vector<std::uint8_t> reordered_codes(codes.size());
    std::vector<std::uint64_t> reordered_labels(options.count);
    for (std::size_t row = 0; row < options.count; ++row) {
        const std::size_t list = assignments[row];
        const std::size_t destination = static_cast<std::size_t>(cursors[list]++);
        std::copy_n(codes.data() + row * bytes, bytes, reordered_codes.data() + destination * bytes);
        reordered_labels[destination] = static_cast<std::uint64_t>(row);
    }

    IndexData data;
    std::copy_n("GINIVF01", 8, data.header.magic);
    data.header.count = options.count;
    data.header.code_bytes = bytes;
    data.header.window_size = options.window_size;
    data.header.pq_m = options.pq_m;
    data.header.nbits = options.nbits;
    data.header.nlist = options.nlist;
    data.prototypes = std::move(prototypes);
    data.offsets = std::move(offsets);
    data.codes = std::move(reordered_codes);
    data.labels = std::move(reordered_labels);
    write_index(options.index, data);
    std::cout << "built WindowIVF: elements=" << options.count
              << " nlist=" << options.nlist
              << " code_bytes=" << bytes
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
    if (!labels_out || !distances_out) throw std::runtime_error("cannot create IVF result files");
    labels_out.write(reinterpret_cast<const char *>(labels.data()), static_cast<std::streamsize>(labels.size() * sizeof(std::uint64_t)));
    distances_out.write(reinterpret_cast<const char *>(distances.data()), static_cast<std::streamsize>(distances.size() * sizeof(float)));
    if (!labels_out || !distances_out) throw std::runtime_error("failed writing IVF result files");
}

void search_index(const Options &options) {
    require_common(options);
    if (options.query_count < 1 || options.k < 1 || options.index.empty() || options.codes.empty() ||
        options.similarity.empty() || options.labels_out.empty() || options.distances_out.empty()) {
        throw std::runtime_error("search requires code/index/similarity paths, query count, k, and output paths");
    }
    const std::size_t bytes = code_bytes(options);
    const std::vector<std::uint8_t> queries = read_codes(options.codes, options.query_count, bytes);
    const std::vector<float> similarity = read_similarity(options.similarity, options.pq_m, ksub(options));
    const IndexData data = read_index(options.index, options);
    const std::size_t nlist = static_cast<std::size_t>(data.header.nlist);
    const std::size_t nprobe = std::min(options.nprobe, nlist);
    const std::size_t result_k = std::min(options.k, static_cast<std::size_t>(data.header.count));
    std::vector<std::uint64_t> labels(options.query_count * options.k, std::numeric_limits<std::uint64_t>::max());
    std::vector<float> distances(options.query_count * options.k, std::numeric_limits<float>::infinity());
    configure_threads(options.threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) if (options.threads > 1)
#endif
    for (std::int64_t row = 0; row < static_cast<std::int64_t>(options.query_count); ++row) {
        const std::uint8_t *query = queries.data() + static_cast<std::size_t>(row) * bytes;
        std::vector<std::pair<float, std::size_t>> coarse;
        coarse.reserve(nlist);
        for (std::size_t list = 0; list < nlist; ++list) {
            coarse.emplace_back(window_score(query, data.prototypes.data() + list * bytes, options, similarity), list);
        }
        std::partial_sort(coarse.begin(), coarse.begin() + nprobe, coarse.end(),
                          [](const auto &left, const auto &right) { return left.first > right.first; });
        std::vector<std::pair<float, std::uint64_t>> candidates;
        for (std::size_t probe = 0; probe < nprobe; ++probe) {
            const std::size_t list = coarse[probe].second;
            const std::size_t begin = static_cast<std::size_t>(data.offsets[list]);
            const std::size_t end = static_cast<std::size_t>(data.offsets[list + 1]);
            for (std::size_t position = begin; position < end; ++position) {
                const float score = window_score(query, data.codes.data() + position * bytes, options, similarity);
                candidates.emplace_back(score, data.labels[position]);
            }
        }
        const std::size_t keep = std::min(result_k, candidates.size());
        std::partial_sort(candidates.begin(), candidates.begin() + keep, candidates.end(),
                          [](const auto &left, const auto &right) { return left.first > right.first; });
        for (std::size_t rank = 0; rank < keep; ++rank) {
            labels[static_cast<std::size_t>(row) * options.k + rank] = candidates[rank].second;
            distances[static_cast<std::size_t>(row) * options.k + rank] = -candidates[rank].first;
        }
    }
    write_results(options.labels_out, options.distances_out, labels, distances);
    std::cout << "searched WindowIVF: queries=" << options.query_count
              << " k=" << options.k
              << " nprobe=" << nprobe
              << " threads=" << options.threads << std::endl;
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
