# Window indexes

GINflow searches sliding windows of GINFINITY embeddings. The `--index` option
selects the library; the library-specific option selects the concrete index
type. A published database contains the index plus `windows.tsv`, `meta.json`,
packed residue embeddings, and records.

The tables below separate parameters used while constructing a database from
parameters that affect queries against it. A parameter marked **fixed** is part
of the current adapter and is not a Nextflow flag. Values shown as `auto` are
resolved from the number of indexed windows and stored in `meta.json`.

## Choose a library: `--index`

| Library | Index-type option | Index types | GPU | On disk |
|---|---|---|---|---|
| `faiss` (default) | `--faiss_index` | `flatip`, `flatl2`, `ivfflat`, `hnsw`, `ivfpq`, `ivfsq`, `pq`, `sq`, `lsh`, `ivfpqr` | Optional for `flatip`, `flatl2`, `ivfflat`, `ivfpq`, and `ivfsq` | `index.faiss` |
| `scann` | none | ScaNN brute-force, asymmetric hashing, or tree + asymmetric hashing plan | No | `scann/` |
| `ngt` | `--ngt_index` | `ngt`, `qg`, `qbg` | No | `ngt/` |
| `cuvs` | `--cuvs_index` | `cagra`, `ivf`, `ivf-pq` | Required | `cuvs/` |
| `hnswlib` | HNSW parameters below | Quantized-node HNSWLIB graph, or optional GPU CAGRA companion | Optional (`--hnswlib_gpu`) | `index.bin` or `cagra/` |

`--database` is a path to an existing published database. On query-only runs,
GINflow reads `meta.json` and detects the library and index type. Passing an
explicit `--index` that disagrees with the database is an error.

By default, windows are 1408-dimensional (`--window_size 11`) and
L2-normalized. FAISS and ScaNN use inner-product scores directly. NGT and cuVS
return distances that the adapter converts to the cosine-like score used by
`--seed_min_similarity`.

The HNSWLIB backend is different: it assigns each 128-dimensional node
embedding to a centroid, persists the full float16 centroid vectors and their
float32 `k x k` similarity matrix, and stores only uint16 centroid codes in
the custom C++ HNSW elements. For a pair of code windows `A` and `B`, its raw
score is exactly `sum(S[A[p], B[p]])`. The public `--seed_min_similarity`
threshold remains on the existing per-position scale, so a non-reranked HNSW
search applies it internally as `window_size * seed_min_similarity`; with
`--hnswlib_rerank true`, the original float16 windows provide the final score
and threshold. The original float16 node embeddings are still packed into the
database and remain the only embeddings used by `ALIGN_CLUSTERS` and
`DRAW_SW`.

## Shared parameters

These controls apply across libraries and are not repeated in every type
table.

### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--index` | `faiss` | Index library used to build and search the database. |
| `--window_size` | `11` | Nucleotides per window; each window concatenates this many 128-dimensional residue embeddings. |
| `--window_stride` | `1` | Step between consecutive window starts. |
| `-profile gpu` | not set | Required by cuVS, by FAISS when `--faiss_gpu` is enabled, and by HNSWLIB when `--hnswlib_gpu true`. It selects the GPU task image and NVIDIA runtime. |

### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Maximum neighbours retained per query window before thresholding. It is also the number of neighbours requested from ScaNN and cuVS CAGRA. |
| `--seed_min_similarity` | `0.8` | Minimum cosine-like similarity required to keep a seed. |
| `--search_shard_size` | `--shard_size` | Query records processed by one search task. |

## Quantized-node HNSWLIB: `--index hnswlib`

HNSWLIB provides approximate candidate selection over centroid-code windows.
The [Python binding](https://github.com/nmslib/hnswlib#python-bindings) exposes
only built-in `ip`, `l2`, and `cosine` spaces. GINflow therefore uses the
official C++ `SpaceInterface`/`DISTFUNC` API instead: the vendored custom space
stores each window as `w` uint16 centroid codes and evaluates the registered
sum of pairwise centroid similarities directly. The index does not contain
expanded `w x 128` vectors.

`faiss/quantization/centroids.npy` stores the `k x 128` float16 codebook and
`similarity.npy` stores its float32 `k x k` lookup table. Those are side
artifacts used for assignment and the custom distance; `faiss/index.bin`
stores the quantized code windows. `faiss/embeddings.npz` remains the original
128-dimensional float16 residue data used by SW and alignment. With
`--hnswlib_rerank true`, HNSW selects a candidate pool and the original windows
are used for the final scores.

The implementation plan, cache layout, validation record, and Rouskin
FlatIP-recall benchmark are maintained in
[quantized-hnsw-research.md](quantized-hnsw-research.md).

The HNSWLIB environment contains the pinned hnswlib 0.8.0 headers/API and a C++
compiler. Use `-profile conda` or `-profile wave` for this backend; the existing
FAISS/ScaNN Docker images do not provide the complete custom C++ build path.
That CPU requirement applies when `--hnswlib_gpu false`; the CAGRA companion
uses the cuVS GPU image with `-profile docker,gpu`. The database build also
copies the original embeddings and graph records needed by the later
SW/alignment steps.

### GPU companion: `--hnswlib_gpu true`

The GPU mode is a companion implementation for the HNSWLIB workflow, not a
custom-distance extension to the hnswlib Python binding. It uses the pinned
cuVS 24.10 CAGRA graph over the original normalized window vectors after a
single global int8 scale. CAGRA's squared-L2 search is only candidate
selection. The candidate labels are exact-reranked with the preserved original
128-dimensional float16 node embeddings and the normal full-window cosine
score before seeds are written.

This design keeps the requested quantized-node artifacts available for the CPU
custom-distance path while avoiding expansion of every candidate into a dense
float vector on the GPU. It also means that `--hnswlib_gpu_candidate_k` must be
at least `--seed_k`, and `--hnswlib_gpu_itopk_size` must be at least the
candidate count. A GPU database is stored under `faiss/cagra/index.bin` and is
identified by `meta.json` as `HNSWLIB_GPU_CAGRA`. It requires `-profile gpu`
for both construction and query-only search.

| Parameter | Default | Description |
|---|---:|---|
| `--hnswlib_gpu` | `false` | Select the CAGRA companion for HNSWLIB build/search. |
| `--hnswlib_gpu_candidate_k` | `50` | Candidates retained before exact original-vector reranking. |
| `--hnswlib_gpu_itopk_size` | `256` | CAGRA intermediate beam. Must be at least `candidate_k`. |
| `--hnswlib_gpu_search_batch_size` | `512` | Query windows per GPU search batch. Lower it for smaller VRAM. |
| `--hnswlib_gpu_intermediate_graph_degree` | `128` | CAGRA build-time graph degree. |
| `--hnswlib_gpu_graph_degree` | `64` | Retained CAGRA graph degree. |
| `--hnswlib_gpu_build_algo` | `nn_descent` | Pinned cuVS graph builder (`nn_descent` or `ivf_pq`). |
| `--hnswlib_gpu_int8_scale` | auto | Optional non-clipping global scale; auto uses the validated 850 scale when safe and lowers it only when needed to avoid clipping. |

The GPU search does not replace the original embeddings: `faiss/embeddings.npz`
remains the SW/alignment source, and `faiss/quantization/` remains the fitted
node codebook/similarity artifact. See the measured Rouskin comparison in
[quantized-hnsw-research.md](quantized-hnsw-research.md).

### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--index` | `hnswlib` | Selects centroid-code quantization plus HNSWLIB. `hnsw` is accepted as an alias. |
| `--node_quantization_k` | `2048` | Requested spherical k-means centroid count. Small inputs use the effective count available. |
| `--node_quantization_sample_size` | `500000` | Maximum node vectors used to fit k-means. |
| `--node_quantization_niter` | `25` | Spherical k-means iterations. |
| `--node_quantization_seed` | `1` | Centroid-fitting seed. |
| `--hnswlib_m` | `32` | HNSW graph degree. |
| `--hnswlib_ef_construction` | `200` | Build-time exploration factor. |
| `--hnswlib_ef_search` | `100` | Default search-time exploration factor, persisted in `meta.json`. |
| `--hnswlib_random_seed` | `1` | HNSW graph construction seed. |
| `--hnswlib_num_threads` | `0` | HNSWLIB worker threads; zero keeps the library default. |
| `--hnswlib_candidate_k` | `50` | Compact code-HNSW candidates retrieved before output or optional original-vector reranking. Must be at least `--seed_k`. |
| `--hnswlib_rerank` | `false` | Rerank candidates with the original float16 windows and use the normal full-window cosine score. |

### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Maximum HNSW neighbours retained per query window before thresholding. |
| `--seed_min_similarity` | `0.8` | With reranking, the full-window cosine threshold; without reranking, the equivalent raw centroid-score sum is used internally. |
| `--hnswlib_ef_search` | `100` | Query-time HNSW exploration factor. |

For the Rouskin high-recall profile, use `--node_quantization_k 4096`,
`--hnswlib_m 32`, `--hnswlib_ef_construction 200`,
`--hnswlib_ef_search 5000`, `--hnswlib_candidate_k 5000`, and
`--hnswlib_rerank true`. The requested default remains `k=2048`; benchmark
trade-offs are recorded in [quantized-hnsw-research.md](quantized-hnsw-research.md).

The published HNSW artifacts are `faiss/index.bin`, `faiss/quantization/`,
`faiss/windows.tsv`, `faiss/embeddings.npz`, `faiss/records.tsv`, and
`faiss/meta.json`. The separate published `quantization/` and
`windows_quantized/` directories make the research representation inspectable;
they are disposable rebuild products, not substitutes for the original
`embeddings/` artifacts.

The following sections contain one build-parameter table and one
search-parameter table for every concrete index type. Parameters that do not
apply to a selected type are ignored; see [unused-parameter warnings](#unused-parameter-warnings).

## NGT: `--index ngt`

[NGT](https://github.com/NGT-labs/NGT) is CPU-only. The adapter exposes the
regular graph (`ngt`), quantized graph (`qg`), and quantized blob graph
(`qbg`). The NGT conda package includes the `qbg` command-line tool used by
the quantized variants.

For the complete upstream command reference, see the [NGT `ngt` README](https://github.com/NGT-labs/NGT/blob/main/bin/ngt/README.md)
and [NGT `qbg` README](https://github.com/NGT-labs/NGT/blob/main/bin/qbg/README.md).
GINflow exposes the important knobs that are supported by its pinned Python
and `qbg` interfaces; the fixed cosine/float representation is intentional
because the window vectors are already L2-normalized.

### NGT regular graph: `--ngt_index ngt`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--ngt_index` | `ngt` | Selects the regular NGT graph. |
| `--ngt_edge_size_for_creation` | `10` | Initial graph edges per node during construction. Higher values can improve recall at build-time and index-size cost. |
| `--ngt_edge_size_for_search` | `40` | Graph edges explored per node during search and stored as the regular graph's search edge setting. |
| `--ngt_num_threads` | `8` | Threads used by `ngtpy` batch insertion/index construction. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Maximum neighbours requested from NGT for each query window. |
| `--ngt_num_of_search_objects` | `20` | Initial candidate objects explored by NGT's two-stage search. Higher values generally improve recall at search-time cost. |
| `--ngt_max_no_of_edges` | unset | Maximum NGT edges explored at search time; when unset, `--ngt_edge_size_for_search` supplies the edge budget. |
| `--ngt_search_range_coefficient` | binding default (`0.1`) | NGT search epsilon/search-range coefficient. |
| `--ngt_search_radius` | unlimited | Optional maximum NGT distance/radius. |
| `--ngt_edge_size_for_search` | `40` | Per-node edge budget used by the Python search setter. |
| `--seed_min_similarity` | `0.8` | Cosine-like threshold after converting NGT cosine distance. |

### NGT quantized graph: `--ngt_index qg`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--ngt_index` | `qg` | Builds a regular graph and quantizes it with `qbg create-qg` and `qbg build-qg`. |
| `--ngt_edge_size_for_creation` | `10` | Initial graph edges per node in the regular graph before quantization. |
| `--ngt_edge_size_for_search` | `40` | Graph edges explored per node by the quantized graph search. |
| `--ngt_num_threads` | `8` | Threads used by the regular graph build. |
| `--ngt_max_no_of_edges` | unset (qbg default) | Maximum edges passed to `qbg build-qg -E` when set. |
| `--ngt_qg_subvector_dimensions` | unset (qbg default) | Subvector dimension passed to `qbg create-qg -Q`; controls the quantization trade-off. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Maximum neighbours requested from the quantized graph for each query window. |
| `--ngt_num_of_search_objects` | `20` | Initial candidate objects explored by the quantized search. |
| `--ngt_max_no_of_edges` | unset | Maximum quantized-graph edges explored at search time; when unset, `--ngt_edge_size_for_search` supplies the edge budget. |
| `--ngt_search_range_coefficient` | binding default (`0.02`) | Quantized graph search epsilon. |
| `--ngt_result_expansion` | binding default (`3.0`) | Expansion ratio for approximate inner-search candidates. Higher values can improve recall at search-time cost. |
| `--ngt_edge_size_for_search` | `40` | Per-node edge budget used by the quantized graph search setter. |
| `--ngt_search_radius` | unlimited | Optional maximum distance/radius. |
| `--seed_min_similarity` | `0.8` | Cosine-like threshold after converting NGT cosine distance. |

### NGT quantized blob graph: `--ngt_index qbg`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--ngt_index` | `qbg` | Builds an NGT quantized blob graph with the `qbg create`, `append`, and `build` commands. |
| `dimension` and `distance` *(fixed)* | window dimension / L2 (`-D 2`) | The adapter creates the QBG store with the input dimension and L2 distance so its distances can be converted back to the cosine-like seed score. |
| `build output count` *(fixed)* | `n_windows` | The adapter passes the number of indexed windows to `qbg build -o`. |
| `--ngt_max_no_of_edges` | unset (`128` at load) | Passes the maximum edge count to `qbg build -E` and to the QBG Python loader's `max_no_of_edges`; unset keeps the qbg build default and the current loader fallback of `128`. |
| `--ngt_qbg_subvectors` | unset (qbg default) | Local quantizer division/subvector count passed to `qbg create -N`. |
| `--ngt_qbg_cluster_data_type` | `pq4` | QBG cluster/codebook storage: `pq4` is compact; `sq8` can improve quantization fidelity at higher storage cost. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Maximum neighbours requested from QBG for each query window. |
| `--seed_min_similarity` | `0.8` | L2 distance is converted to a cosine-like score before thresholding. |
| `--ngt_max_no_of_edges` | `128` at load | Maximum QBG edges allocated by the Python loader and explored at search time. |
| `--ngt_num_of_search_objects` | `20` | Initial QBG search-object count. |
| `--ngt_edge_size_for_search` | `40` | QBG edge exploration budget used by the Python search setter. |
| `--ngt_search_range_coefficient` | binding default | QBG graph-stage search epsilon. |
| `--ngt_blob_search_range_coefficient` | binding default | QBG blob-stage search epsilon. |
| `--ngt_result_expansion` | `0.0` | QBG approximate-result expansion factor. |
| `--ngt_exploration_size` | `256` | Number of QBG nodes explored by the blob search. |
| `--ngt_exact_result_expansion` | `0.0` | QBG exact-result expansion factor. |
| `--ngt_num_of_probes` | `1` | QBG probe count; one probe avoids empty results on small collections with current NGT releases. |
| `--ngt_search_radius` | unlimited | Optional QBG search radius. |

NGT databases are published under `results/faiss/ngt/` for compatibility with
the rest of the pipeline. A later `--query --database results/faiss` run
detects NGT from `meta.json`.

## cuVS: `--index cuvs`

[cuVS](https://docs.nvidia.com/cuvs/) is NVIDIA's GPU-optimized vector-search
library. It requires `-profile gpu` and a visible NVIDIA GPU for both building
and searching. GINflow pins cuVS `24.10.00` with CUDA 11.8 and CuPy 13.0.0.

The [NVIDIA CAGRA guide](https://docs.nvidia.com/cuvs/user-guide/api-guides/indexing-guide/cagra)
lists the full upstream build and search API. GINflow exposes the important
graph-degree, graph-builder, and intermediate-search controls below. The
adapter fixes the metric to cosine and serializes the dataset with the index
attached so a later search-only run can load the database.

### cuVS CAGRA: `--cuvs_index cagra`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--cuvs_index` | `cagra` | Selects the GPU graph index. |
| `--cuvs_intermediate_graph_degree` | `128` | Number of neighbours in the initial graph before pruning. Higher values can improve recall but increase build time and memory. Values are capped at `n_windows - 1`. |
| `--cuvs_graph_degree` | `64` | Number of neighbours retained in the final graph. Higher values improve recall at graph-memory and search cost. It cannot exceed the intermediate degree. |
| `--cuvs_build_algo` | `nn_descent` | Initial graph builder. Choices are `ivf_pq`, `nn_descent`, `iterative_cagra_search`, and `ace`. |
| `--cuvs_itopk_size` | `64` | Intermediate candidate count persisted with the CAGRA database for later searches. The runtime raises it to at least `--seed_k`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--cuvs_itopk_size` | `64` | Intermediate candidate count for CAGRA search. It is stored in `meta.json` at build time and is raised to at least `--seed_k` for an actual query. Increase it for recall at throughput cost. |
| `--seed_k` | `50` | Final neighbours requested from CAGRA per query window. |
| `--seed_min_similarity` | `0.8` | Cosine-like threshold after converting cuVS cosine distance. |

The upstream CAGRA API also has options such as `metric_arg`, `compression`,
`graph_build_params`, and `guarantee_connectivity`. They are not silently
translated into GINflow flags: the current adapter fixes or omits them, as
shown above. Extending those options requires a matching cuVS adapter and
schema change.

### cuVS IVF-Flat: `--cuvs_index ivf`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--cuvs_index` | `ivf` | Selects cuVS IVF-Flat. |
| `--cuvs_n_lists` | `min(4 * sqrt(n), n)` | Number of coarse IVF lists. The resolved value is capped at the number of indexed windows. |
| `--cuvs_n_probes` | `min(20, n_lists)` | Lists searched per query by default. The resolved value is stored in `meta.json` and can be overridden for a later search. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--cuvs_n_probes` | stored build value | Number of IVF lists searched per query. Higher values usually improve recall and reduce throughput. It is capped at `n_lists`. |
| `--seed_k` | `50` | Maximum neighbours requested from IVF-Flat per query window. |
| `--seed_min_similarity` | `0.8` | Cosine-like threshold after converting cuVS cosine distance. |

### cuVS IVF-PQ: `--cuvs_index ivf-pq`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--cuvs_index` | `ivf-pq` | Selects cuVS IVF with product quantization. |
| `--cuvs_n_lists` | `min(4 * sqrt(n), n)` | Number of coarse IVF lists. The resolved value is capped at the number of indexed windows. |
| `--cuvs_n_probes` | `min(20, n_lists)` | Default lists searched per query. The resolved value is stored in `meta.json`. |
| `--cuvs_pq_bits` | `8` | Bits per product-quantization code. Supported values are `4`, `5`, `6`, `7`, and `8`. |
| `--cuvs_pq_dim` | `0` | Compressed dimension. `0` delegates the choice to cuVS; otherwise choose a dimension compatible with the window dimension. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--cuvs_n_probes` | stored build value | Number of IVF lists searched per query. Higher values usually improve recall and reduce throughput. |
| `--seed_k` | `50` | Maximum neighbours requested from IVF-PQ per query window. |
| `--seed_min_similarity` | `0.8` | Cosine-like threshold after converting cuVS cosine distance. |

cuVS databases are stored in `results/faiss/cuvs/index.bin`. Search-only runs
detect `backend: "cuvs"` from `meta.json`, but still require `-profile gpu`.

## FAISS: `--index faiss`

[FAISS index summary](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes)
is a summary of the methods and provides further documentation. GINflow pins
FAISS 1.10.0. The on-disk `index.faiss` is always a CPU index; `--faiss_gpu`
uses the GPU only inside the build or search task.

### FAISS type overview

| Type | Metric | Exact? | GPU | Typical use |
|---|---|---:|---:|---|
| `flatip` | Inner product | Yes | Yes | Default exact cosine search. |
| `flatl2` | Squared L2 | Yes | Yes | Exact search with L2-to-cosine conversion. |
| `ivfflat` | Inner product over probed IVF lists | No | Yes | First approximation to try for large databases. |
| `hnsw` | Inner product graph walk | No | No | High recall without IVF training. |
| `ivfpq` | IVF + product quantization | No | Yes | Lower serving memory for large databases. |
| `ivfsq` | IVF + scalar quantization | No | Yes | Simpler scalar compression with IVF search. |
| `pq` | Product quantization | No | No | Compression without IVF partitioning. |
| `sq` | Scalar quantization | No | No | Light compression without IVF partitioning. |
| `lsh` | Hamming distance | No | No | Cheap, low-memory filtering with lower recall. |
| `ivfpqr` | IVF + PQ + refinement | No | No | Extra PQ refinement; L2-only in pinned FAISS. |

`pq`, `ivfpq`, and `ivfpqr` need at least `2^n` training windows, where `n` is
the value of `--faiss_pq_nbits`. With the default `n=8`, that is 256 windows.

### FAISS `flatip`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `flatip` | Selects exact inner-product search. |
| `--faiss_gpu` | `false` | Build on GPU when `true`; requires `-profile gpu`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_gpu` | `false` | Load and search the CPU index on GPU; requires `-profile gpu`. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |

### FAISS `flatl2`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `flatl2` | Selects exact squared-L2 search. |
| `--faiss_gpu` | `false` | Build on GPU when `true`; requires `-profile gpu`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_gpu` | `false` | Load and search the CPU index on GPU; requires `-profile gpu`. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine-like similarity after L2 conversion. |

### FAISS `ivfflat`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `ivfflat` | Selects IVF with uncompressed vectors. |
| `--faiss_nlist` | `max(1, min(floor(4 * sqrt(n)), n, floor(n / 39) if n >= 39 else n))` | Number of coarse IVF centroids. The automatic value is capped to leave enough training vectors per centroid. |
| `--faiss_nprobe` | `min(8, nlist)` | Initial number of IVF lists searched; stored in `meta.json`. |
| `--faiss_gpu` | `false` | Build on GPU when `true`; requires `-profile gpu`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_nprobe` | stored build value | IVF lists searched per query. Higher values usually improve recall and reduce throughput. |
| `--faiss_gpu` | `false` | Search on GPU when `true`; requires `-profile gpu`. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |

### FAISS `hnsw`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `hnsw` | Selects an inner-product HNSW graph. |
| `--faiss_hnsw_m` | `32` | Maximum graph neighbours per node. Higher values improve recall and increase memory. |
| `--faiss_hnsw_ef_construction` | `40` | Build-time graph exploration depth. Higher values can improve graph quality at build cost. |
| `--faiss_hnsw_ef_search` | `16` | Initial search exploration depth saved in the index. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_hnsw_ef_search` | stored build value | Search exploration depth. Higher values usually improve recall and reduce throughput. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |

### FAISS `ivfpq`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `ivfpq` | Selects IVF with product-quantized vectors. |
| `--faiss_nlist` | `max(1, min(floor(4 * sqrt(n)), n, floor(n / 39) if n >= 39 else n))` | Number of coarse IVF centroids. |
| `--faiss_nprobe` | `min(8, nlist)` | Initial IVF probe count; stored in `meta.json`. |
| `--faiss_pq_m` | `16` | Number of PQ subvectors; must divide the window dimension. |
| `--faiss_pq_nbits` | `8` | Bits per PQ subquantizer: `4`, `8`, `12`, or `16`. |
| `--faiss_gpu` | `false` | Build on GPU when `true`; requires `-profile gpu`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_nprobe` | stored build value | IVF lists searched per query. |
| `--faiss_gpu` | `false` | Search on GPU when `true`; requires `-profile gpu`. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |

### FAISS `ivfsq`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `ivfsq` | Selects IVF with scalar-quantized vectors. |
| `--faiss_nlist` | `max(1, min(floor(4 * sqrt(n)), n, floor(n / 39) if n >= 39 else n))` | Number of coarse IVF centroids. |
| `--faiss_nprobe` | `min(8, nlist)` | Initial IVF probe count; stored in `meta.json`. |
| `--faiss_sq_type` | `8bit` | Scalar quantizer width: `8bit`, `6bit`, `4bit`, or `fp16`. |
| `--faiss_gpu` | `false` | Build on GPU when `true`; requires `-profile gpu`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_nprobe` | stored build value | IVF lists searched per query. |
| `--faiss_gpu` | `false` | Search on GPU when `true`; requires `-profile gpu`. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |

### FAISS `pq`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `pq` | Selects standalone product quantization without IVF. |
| `--faiss_pq_m` | `16` | Number of PQ subvectors; must divide the window dimension. |
| `--faiss_pq_nbits` | `8` | Bits per PQ subquantizer: `4`, `8`, `12`, or `16`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |
| Index-specific search tuning | none | FAISS PQ has no separate GINflow search override; search uses the built codebooks. |

### FAISS `sq`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `sq` | Selects standalone scalar quantization without IVF. |
| `--faiss_sq_type` | `8bit` | Scalar quantizer width: `8bit`, `6bit`, `4bit`, or `fp16`. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity to emit as a seed. |
| Index-specific search tuning | none | FAISS SQ has no separate GINflow search override; search uses the built quantizer. |

### FAISS `lsh`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `lsh` | Selects FAISS locality-sensitive hashing. |
| `--faiss_lsh_nbits` | `2 * window_dim` | Hamming code length. Larger codes increase memory and can change recall. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum converted similarity to emit as a seed. |
| Index-specific search tuning | none | FAISS LSH has no separate GINflow search override. |

### FAISS `ivfpqr`

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_index` | `ivfpqr` | Selects IVF-PQ with a refinement code. |
| `--faiss_nlist` | `max(1, min(floor(4 * sqrt(n)), n, floor(n / 39) if n >= 39 else n))` | Number of coarse IVF centroids. |
| `--faiss_nprobe` | `min(8, nlist)` | Initial IVF probe count; stored in `meta.json`. |
| `--faiss_pq_m` | `16` | Number of primary PQ subvectors; must divide the window dimension. |
| `--faiss_pq_nbits` | `8` | Bits per primary PQ subquantizer: `4`, `8`, `12`, or `16`. |
| `--faiss_pq_m_refine` | `4` | Number of extra PQ subvectors used for refinement. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--faiss_nprobe` | stored build value | IVF lists searched per query. |
| `--seed_k` | `50` | Neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum converted similarity to emit as a seed. |
| Index-specific search tuning | none | Refinement uses the codebooks and `--faiss_pq_m_refine` selected at build time. |

Example GPU FAISS build:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --index faiss \
    --faiss_index ivfflat \
    --faiss_gpu \
    --input structures.tsv \
    --outdir results
```

## ScaNN: `--index scann`

[ScaNN](https://github.com/google-research/google-research/tree/master/scann)
(`scann==1.4.2`) is CPU-only and uses AVX/FMA on x86. ScaNN has no separate
index-type flag in GINflow: the builder selects a brute-force, asymmetric
hashing, or partitioned tree + asymmetric hashing plan from the database size
and `--scann_leaves`.

| Database size | Effective plan when `--scann_leaves` is unset |
|---|---|
| `< 20,000` windows | Brute-force exact scoring |
| `20,000`–`99,999` windows | Asymmetric hashing plus exact reorder |
| `>= 100,000` windows | Spherical tree, asymmetric hashing, and exact reorder |
| Any size with `--scann_leaves` | Spherical tree, asymmetric hashing, and exact reorder |

### ScaNN searcher

#### Build parameters

| Parameter | Default | Description |
|---|---:|---|
| `--index` | `scann` | Selects the ScaNN backend. |
| `--scann_leaves` | auto `round(sqrt(n))` for `n >= 100,000` | Number of tree partitions. Setting it forces the tree plan, including for small databases. |
| `--scann_ah_dim` | `2` | Dimensions per asymmetric-hash block. The ScaNN guidance recommends `2`. |
| `--scann_anisotropic` | `0.2` | Anisotropic quantization threshold used by asymmetric hashing. |
| `--scann_soar` | `false` | Enables SOAR spilling when a tree plan is built. |
| `--scann_leaves_to_search` | auto `min(8, leaves)` | Number of tree partitions searched by the persisted tree plan; a query-time value can override it. |
| `--scann_reorder` | `max(100, --seed_k)` | Exact-reorder pool after asymmetric hashing; the effective value is raised to at least `--seed_k`. |
| `--seed_k` | `50` | Number of neighbours passed to the ScaNN builder. |

#### Search parameters

| Parameter | Default | Description |
|---|---:|---|
| `--scann_leaves_to_search` | auto `min(8, leaves)` | Tree partitions scored per query. Only applies to a partitioned plan. |
| `--scann_reorder` | stored build value | Exact-reorder pool; a query-time value overrides the stored value. |
| `--seed_k` | `50` | Final neighbours requested per query window. |
| `--seed_min_similarity` | `0.8` | Minimum inner-product similarity to emit as a seed. |

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --index scann \
    --scann_leaves 4096 \
    --scann_leaves_to_search 32 \
    --scann_reorder 200 \
    --input structures.tsv \
    --query queries.tsv \
    --outdir results
```

Smoke-test and other small databases stay brute-force unless
`--scann_leaves` is set.

## Unused-parameter warnings

At launch, GINflow checks parameters explicitly supplied on the command line
or changed from their schema defaults. Library-inapplicable parameters and
FAISS-type-inapplicable parameters produce a warning and are ignored. For
cuVS, only the parameters consumed by the selected cuVS type are meaningful;
the current adapter records the applicable values and ignores the rest.

```
WARN: Unused index parameters for --index faiss / --faiss_index flatip (ignored): --faiss_nlist, --scann_reorder. See docs/indexes.md.
```

Hard errors (not warnings):

- `--index` not `faiss`, `scann`, `ngt`, or `cuvs`
- `--faiss_index scann` (use `--index scann`)
- `--index` disagrees with an existing `--database`
- `--faiss_gpu` without `-profile gpu`
- `--faiss_gpu` with `--index scann`, `ngt`, `cuvs`, or a CPU-only FAISS type
- `--index cuvs` without `-profile gpu`

## Capacity planning and tuning guide

This section is deliberately explicit about evidence. **Implementation fact**
means that the current GINflow source enforces it. **Source-backed** means an
upstream project's documentation or original paper says it. **Derived** is
arithmetic from the window count and dimension. **Planning heuristic** is a
conservative starting point to validate with a measurement; it is not a
benchmark result or a vendor hardware guarantee.

### Version scope

These rules apply to the images pinned by this repository, rather than to the
latest release mentioned by an upstream web page:

| Library | GINflow version | Scope of this guide |
|---|---:|---|
| FAISS CPU / GPU | 1.10.0 | The classes and GPU support exposed by `--faiss_index` / `--faiss_gpu`. |
| ScaNN | 1.4.2 | The Python builder path used here: dot product, spherical tree, AH, optional exact reorder. |
| NGT | 2.3.12 | The `ngtpy` and `qbg` paths used by the NGT adapter. |
| cuVS | 24.10.00 | The CUDA 11 Python API invoked by the cuVS adapter. Newer cuVS documentation may use different defaults or APIs. |

Treat post-24.10 cuVS documentation as conceptual background only until a
setting has been checked against the pinned adapter and a small build. In
particular, do not copy a current cuVS example's undocumented default into a
24.10 experiment. The pinned 24.10 branch identifies this release as the
first fully featured migration of the ANN implementations from RAFT
([NVIDIA cuVS 24.10 branch](https://github.com/NVIDIA/cuvs/tree/branch-24.10)).

### Dataset model and raw-vector lower bound

**Implementation fact:** the standard window is 1,408 float values (11
residues × 128) and the index builder materializes them as contiguous
`float32`. Let `N` be the final window count, `D = 1408`, and `B` the bytes
per component. The raw vector lower bound is:

```
raw_bytes = N × D × B
float32: B = 4                 float16: B = 2
```

This excludes the index topology, IDs, centroids/codebooks, metadata, query
buffers, Python objects, the packed residue embeddings, and transient build
buffers. It is therefore a *lower bound*, not a machine requirement.

The `float16` column below is a storage comparator only. The current GINflow
builder materializes the common input matrix as `float32` for every backend.

The planned test inputs, with the default 11-nt window and stride 1, project
to the following counts. They are **derived estimates** from the supplied TSVs;
the completed build manifest is the source of truth if short records are
skipped or windowing settings change.

| Input | Projected `N` | Raw `float32` | Raw `float16` | `sqrt(N)` | GINflow auto IVF `4 × sqrt(N)` |
|---|---:|---:|---:|---:|---:|
| `tests/data/rouskin_sample_6k.tsv` | 839,188 | 4.402 GiB (4.726 GB) | 2.201 GiB | 916 | 3,664 |
| `tests/data/rouskin_sample_30k.tsv` | 4,115,576 | 21.587 GiB (23.179 GB) | 10.794 GiB | 2,029 | 8,114 |

**Implementation fact:** GINflow's current builder reads shard arrays, keeps
them in a Python list, then concatenates them into one `xb` matrix before
calling a backend. It is not a streaming index builder. During a large build,
allow at least two raw CPU copies (the shard arrays plus `xb`), then add the
index and process overhead. A **planning heuristic** of `3 × raw_float32`
before index-specific overhead is a safer first reservation for this
implementation.

That yields about 13.2 GiB before index overhead for 6k and 64.8 GiB before
index overhead for 30k. The window mapping is also held as Python tuples, so
the 30k case is not a sensible first run on a nominal 64 GiB machine. Prefer
96–128 GiB system RAM for full 30k builds, or measure a reduced subset first.

### Hardware planning by backend

The table is intentionally conservative. “Recommended” is a **planning
heuristic**, not a claim that a smaller machine cannot work. CPU RAM includes
the non-streaming builder behaviour above; VRAM applies to a single visible
GPU and excludes other workloads. Disk should also have space for the raw
window cache, a database, work directories, and repeated benchmark runs.

| Backend / structure | 6k build recommendation | 30k build recommendation | Important constraint |
|---|---|---|---|
| FAISS CPU Flat / IVF-Flat / HNSW | 32 GiB RAM | 96–128 GiB RAM | Flat-class storage starts at the raw `float32` size; HNSW adds graph links. |
| FAISS CPU IVF-PQ / PQ / SQ | 32 GiB RAM | 96 GiB RAM | Compression reduces the *final index*, not the current builder's raw training/input buffers. |
| FAISS GPU | 32 GiB RAM + 12–16 GiB VRAM | 96 GiB RAM + 32 GiB+ VRAM | GPU build/search copies the index/data to device; 8 GiB has too little robust headroom for 6k Flat/IVF and cannot hold the 30k raw matrix. |
| ScaNN | 32 GiB RAM, AVX/FMA CPU | 96–128 GiB RAM, AVX/FMA CPU | Current builder receives the full `float32` matrix and exact reorder needs raw-vector access. CPU-only. |
| NGT / QG / QBG | 32 GiB RAM | 96–128 GiB RAM | Final graph/quantized artifact size is data- and implementation-dependent; do not assume NGT memory mapping because GINflow does not enable it. |
| cuVS CAGRA / IVF / IVF-PQ | 32 GiB RAM + 16 GiB VRAM (24 GiB preferred) | 96 GiB RAM + 32 GiB VRAM minimum (48 GiB preferred) | The adapter uploads the full `float32` dataset with `cupy.asarray` even for IVF-PQ; graph/training workspace is extra. |

For the 8 GiB GPU used by the benchmark host, the 6k raw matrix alone is
4.40 GiB and the 30k raw matrix alone is 21.59 GiB. **Derived conclusion:**
the full 30k cuVS database, and GPU-FAISS Flat or IVF-Flat databases, cannot
fit in 8 GiB. This calculation alone does not determine feasibility for
compressed GPU-FAISS IVFPQ or IVFSQ, whose device-resident structures are
smaller; measure those configurations separately. The 6k cuVS/GPU case may
also fail during construction because build workspace and the index are
additional allocations; treat it as an intentionally infeasible test unless a
small pilot proves otherwise. A compressed persisted IVF-PQ index does not
remove cuVS's full-dataset build-time upload requirement in the current
adapter.

**Source-backed:** FAISS documents all of its indexes as RAM-resident and
lists the underlying bytes/vector formulas; its GPU path is a separate copy
of the CPU index ([FAISS index summary](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes),
[FAISS on GPU](https://github.com/facebookresearch/faiss/wiki/Faiss-on-the-GPU)).
NGT advertises a memory-mapped option for objects larger than memory, but that
is an upstream capability rather than a setting used by this pipeline
([NGT README](https://github.com/NGT-labs/NGT)).

### FAISS memory model

The following are lower bounds from FAISS' index summary; the scalar-code
widths, including SQ6, are also verified against the pinned FAISS 1.10 build.
Centroids/codebooks and container/process overhead are omitted. They describe
the final stored code, not peak GINflow build memory.

| FAISS family | Approximate bytes/vector | Consequence at `D = 1408` |
|---|---:|---|
| `flatip`, `flatl2` | `4 × D` | 5,632 B; 4.402 GiB / 21.587 GiB for the two datasets. |
| `ivfflat` | `4 × D + 8` | Raw vector plus 64-bit ID; nearly Flat-sized. |
| `hnsw` | `4 × D + x × M × 2 × 4` | Raw vector plus graph links; `x` depends on the graph layers. Higher `M` costs more memory. |
| `sq` | `D` for SQ8, `ceil(6D / 8)` for SQ6, `ceil(D / 2)` for SQ4, `2D` for fp16 | Standalone scalar-quantizer codes; SQ8 is 1,408 B/vector. |
| `ivfsq` | scalar-quantizer code + 8 | The same code widths as `sq`, plus a 64-bit IVF ID. |
| `pq` | `ceil(M × nbits / 8)` | `M=16, nbits=8` is 16 B/vector; `M=32` is 32 B/vector; `M=64` is 64 B/vector. |
| `ivfpq` | `ceil(M × nbits / 8) + 8` | PQ code plus ID, with centroids/codebooks in addition. |
| `ivfpqr` | code + refinement code + ID | More storage than IVFPQ for reranking. |

For example, 16-byte PQ codes occupy only 12.8 MiB for 6k or 62.8 MiB for
30k, before IDs and codebooks. This is why compression can dramatically
reduce *serving* RAM, but it does not prove that the full build will fit in
the same amount of RAM.

### General tuning protocol

1. Build a deterministic exact `flatip` reference on the same normalized
   vectors. It is the ground truth for recall@100, not an ANN contender.
2. Use a representative, fixed query set and request `k=100` from every
   backend. Measure recall as overlap with the exact top-100, and record
   build time, index bytes, peak host RAM, peak VRAM, warm query QPS, and
   p50/p95 latency separately.
3. Change one query-time knob at a time after a fixed build whenever possible.
   This makes the recall/QPS curve attributable: `nprobe`,
   `leaves_to_search`, and `n_probes` are query-time ladders.
4. Stop a ladder once it has enough points above recall 0.80 and no longer
   buys useful QPS. Keep failures/infeasible configurations in the raw table,
   but do not place them on the requested recall > 0.80 plot.
5. Never compare cold-load latency with warm search throughput. Record the
   batch size and whether transfer/load was excluded from timing.

The grids below are **starting grids**, not reported results. They use only
flags accepted by the current schema and adapters. For the full 30k database,
prune a grid early after measuring capacity; a clean infeasibility record is
more useful than an out-of-memory retry loop.

## Backend-specific tuning

### FAISS: choose precision, then the query-time ladder

**Source-backed:** FAISS says an IVF search scans roughly `nprobe / nlist` of
the collection (often an underestimate due to uneven lists), recommends an
`nlist` proportional to `sqrt(n)`, and exposes `nprobe` at search time for the
speed/accuracy trade-off. It describes higher HNSW `M` as more accurate and
more memory hungry; `efConstruction` and `efSearch` control construction and
search exploration ([FAISS index summary](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes)).

**Implementation fact:** GINflow's automatic list count is
`max(1, min(floor(4 × sqrt(n)), floor(n / 39) when n >= 39 else n, n))`,
yielding 3,664 and 8,114 above. The automatic `nprobe` is `min(8, nlist)`.
`--faiss_pq_m` must divide the actual window dimension; 16, 32, and 64 all
divide 1,408. Only FlatIP, FlatL2, IVFFlat, IVFPQ, and IVFSQ variants selected
by the table above can use `--faiss_gpu`.

| Goal | First choice | Tune this first | What increases / decreases |
|---|---|---|---|
| Ground truth / safest recall | `flatip` | none | Exact; cost rises linearly with `N`. |
| Fast approximate, retain full vectors | `ivfflat` | `--faiss_nprobe` | Higher probe count tends to increase recall and latency. |
| Fit a much smaller serving index | `ivfpq` | `--faiss_nprobe`, then `--faiss_pq_m` / `--faiss_pq_nbits` | More probes help recall; more code bytes usually help approximation quality and use more memory. |
| High recall CPU graph | `hnsw` | `--faiss_hnsw_ef_search`, then `--faiss_hnsw_m` and `--faiss_hnsw_ef_construction` | Larger values use more memory/build time or search time. |
| Lightweight compression | `ivfsq` or `sq` | `--faiss_sq_type`, then `nprobe` for IVF | Lower-bit scalar codes reduce storage with a recall trade-off. |

Suggested FAISS benchmark ladders (with `k=100`):

| Dataset | Build candidates | Query-time ladder |
|---|---|---|
| 6k | IVFFlat `nlist=2048,4096`; IVFPQ with `nlist=2048,4096`, `pq_m=16,32,64`, `pq_nbits=4,8`; HNSW `M=16,32,48` | IVF `nprobe=8,16,32,64,128`; HNSW `ef_search=16,32,64,128` (start `ef_construction=80`, then test 40/160 only if needed). |
| 30k | IVFFlat `nlist=4096,8192`; IVFPQ with `nlist=4096,8192`, `pq_m=16,32,64`, `pq_nbits=4,8`; HNSW only if RAM permits | IVF `nprobe=16,32,64,128,256`; HNSW `ef_search=32,64,128,256`. |

Keep the `pq_m` / bit sweep small: a full cross-product with every probe
value is expensive and obscures the curve. Start with `pq_m=32, nbits=8`, then
change compression only after locating a useful probe range. Current FAISS
allows 4, 8, 12, or 16 PQ bits; the grid avoids 12/16 initially because their
larger codebooks raise training cost without being necessary to establish a
recall > 0.80 frontier.

### ScaNN: tree coverage and exact reorder

**Source-backed:** ScaNN's official rules say brute-force below 20k points,
AH plus reorder below 100k, and partitioning + AH + reorder above 100k. It
recommends `dimensions_per_block=2`, approximately `sqrt(n)` leaves, tuning
the `leaves_to_search` ratio to the recall target, and an exact-reorder pool
larger than final `k` ([ScaNN algorithms and configuration](https://github.com/google-research/google-research/blob/master/scann/docs/algorithms.md)).
The anisotropic quantization motivation is described in the original
[ScaNN paper](https://arxiv.org/abs/1908.10396).

**Implementation fact:** both supplied datasets are above 100k, so GINflow
uses a spherical tree, asymmetric hashing (AH), and exact reorder by default.
It derives `round(sqrt(n))` leaves: 916 for 6k and 2,029 for 30k. The default
`--scann_leaves_to_search` is 8, `--scann_reorder` is 100, and the reorder
pool is raised to at least `--seed_k`. The build stores initial
`leaves_to_search` and reorder values; `SEARCH_SCANN` can override either on a
later query run.

| Dataset | Stable build baseline | Query-time ladder |
|---|---|---|
| 6k | `--scann_leaves 1024 --scann_ah_dim 2 --scann_anisotropic 0.2` | `--scann_leaves_to_search 8,16,32,64,128`; then `--scann_reorder 100,200,400,800`. |
| 30k | `--scann_leaves 2048 --scann_ah_dim 2 --scann_anisotropic 0.2` | `--scann_leaves_to_search 16,32,64,128,256`; then `--scann_reorder 100,200,400,800`. |

Use leaves 512/2,048 as a 6k build sensitivity and 1,024/4,096 as a 30k
sensitivity only after the baseline ladder. Test `--scann_soar` as a separate
build variant, not mixed into every condition. Keep `--scann_ah_dim 2` fixed
for the primary comparison because that is the upstream rule of thumb; a
small `--scann_anisotropic 0.1,0.2,0.3` sensitivity is optional after the
main curve is complete.

### NGT: compare structures and exposed knobs

**Implementation fact:** GINflow exposes the important NGT construction and
search controls documented in the NGT section above, including graph edge
sizes, search exploration, QBG probing, and the quantizer edge/subvector
settings. The default QBG search configuration remains deliberately
small-collection-safe (`num_of_search_objects=20`, `exploration_size=256`,
`num_of_probes=1`). A useful baseline ladder is:

| Dataset | Valid benchmark variants | Interpretation |
|---|---|---|
| 6k and 30k | `--ngt_index ngt`, `qg`, `qbg` | Compare index type, index size, build time, recall@100, and QPS with the documented defaults. |
| 6k and 30k | `--ngt_edge_size_for_creation 10,20`, `--ngt_edge_size_for_search 40,80` | Measure graph-quality and search-recall sensitivity for regular NGT and QG. |
| 6k and 30k | `--ngt_max_no_of_edges 64,128,256` | Measure QG/QBG graph compactness and search recall; keep the other parameters fixed. |
| 6k and 30k | `--ngt_num_of_search_objects 20,50,100`, `--ngt_exploration_size 128,256,512` | Measure search recall/QPS sensitivity; `--ngt_exploration_size` applies to QBG. |

**Source-backed:** NGT identifies these as graph/tree (NGT), quantized graph
(QG), and quantized blob graph (QBG) methods ([NGT README](https://github.com/NGT-labs/NGT)).
Their compactness and recall must be measured on the RNA vectors; do not
extrapolate a memory ranking from another dataset. If a desired recall target
cannot be reached, the correct action is to report the limitation or add a
deliberate new pipeline parameter in a separate change—not to pass an
undocumented flag to the benchmark script.

### cuVS: GPU capacity first, then graph or IVF coverage

**Implementation fact:** GINflow supports CAGRA, IVF-Flat, and IVF-PQ only.
For CAGRA it forwards `--cuvs_intermediate_graph_degree`,
`--cuvs_graph_degree`, and `--cuvs_build_algo` to the build adapter; final
graph degree cannot exceed intermediate degree. `--cuvs_itopk_size` is also a
build-pipeline setting: `BUILD_CUVS_INDEX` resolves and stores it in
`meta.json`, then `SEARCH_CUVS` reloads that value. It has no search-process
override, so changing it through GINflow requires rebuilding the database.
For IVF variants, `--cuvs_n_lists`, `--cuvs_pq_bits`, and `--cuvs_pq_dim` are
build-only. `--cuvs_n_probes` seeds metadata at build time and is the one cuVS
setting that `SEARCH_CUVS` can override later. The current adapter's automatic
list count is `min(n, floor(4 × sqrt(n)))`; automatic probes are
`min(20, n_lists)`.

**Version caveat:** the options below are valid GINflow parameter names, but
the exact quality/capacity effect must be established with cuVS 24.10, not
inferred from a newer release. NVIDIA describes cuVS as GPU ANN software and
publishes CAGRA as its graph-search research basis
([cuVS 24.10 branch](https://github.com/NVIDIA/cuvs/tree/branch-24.10),
[CAGRA paper](https://arxiv.org/abs/2308.15136)).

| Structure | Capacity gate | Build variants (each requires rebuilding) | Query-time ladder |
|---|---|---|---|
| CAGRA | Full raw dataset + graph + workspace must fit VRAM | `intermediate/graph=128/64`, `256/64`, `256/128`; `--cuvs_itopk_size 100,128,256`; begin with `--cuvs_build_algo nn_descent` | No CAGRA query override is exposed by the current pipeline. |
| IVF-Flat | Full raw dataset + IVF workspace must fit VRAM | `--cuvs_n_lists 2048,4096` (6k) or `4096,8192` (30k) | `--cuvs_n_probes 8,16,32,64,128`. |
| IVF-PQ | Build still uploads raw data; final index may be smaller | Same lists; `--cuvs_pq_bits 4,6,8 --cuvs_pq_dim 0` | Same probe ladder. Try explicit `pq_dim` only after a 24.10 smoke build accepts it. |

For the required `k=100`, the adapter creates an effective CAGRA search
parameter of at least 100 even if the stored `itopk_size` is lower. Thus a
stored `64` versus `100` comparison is not meaningful for a top-100 benchmark;
start at 100. Values above 100 must be built into separate GINflow databases,
because the search process has no `--cuvs_itopk_size` flag. Keep `nn_descent`
as the main CAGRA build algorithm because it is the pipeline default.
Treat the other schema-listed build algorithms as optional 24.10 compatibility
experiments, not required benchmark points.

### Sources used for these rules

- [FAISS index summary and bytes/vector formulas](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes)
- [FAISS selection guidelines](https://github.com/facebookresearch/faiss/wiki/Guidelines-to-choose-an-index)
- [FAISS GPU documentation](https://github.com/facebookresearch/faiss/wiki/Faiss-on-the-GPU)
- [ScaNN algorithms and configuration](https://github.com/google-research/google-research/blob/master/scann/docs/algorithms.md)
- [ScaNN anisotropic quantization paper](https://arxiv.org/abs/1908.10396)
- [NGT `ngt` command reference](https://github.com/NGT-labs/NGT/blob/main/bin/ngt/README.md)
- [NGT `qbg` command reference](https://github.com/NGT-labs/NGT/blob/main/bin/qbg/README.md)
- [NVIDIA cuVS 24.10 source branch](https://github.com/NVIDIA/cuvs/tree/branch-24.10)
- [CAGRA original paper](https://arxiv.org/abs/2308.15136)
