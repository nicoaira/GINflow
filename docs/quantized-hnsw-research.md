# Quantized-node HNSW research plan, implementation, and results

Status: implemented and validated on the smoke pipeline; the CPU compact path
exceeds 97% FlatIP recall with the recommended candidate-rerank profile, the
optional GPU CAGRA companion reaches 99.0% R@50, a full CAGRA-to-cuVS-CPU-
HNSW conversion/search experiment has completed on the Rouskin 6k set, and
FAISS's true scalar-quantized HNSW path has been measured against the same
int8 representation, and a node-level product-quantization HNSW sweep has
been measured against the same FlatIP labels. The complete documented nf-test regression suite and
the full default k=2048 test-corpus build have passed. Larger downstream
comparisons and whole-pipeline resource measurements remain follow-up work.

This document is deliberately comprehensive. It records the design, the
temporary benchmark layout, the measured trade-offs, and the cleanup TODOs so
that large research caches under `/mnt/ssd_samsung` can be removed or rebuilt
without reconstructing the investigation.

## Objective and invariant

The HNSWLIB backend adds a compact candidate selector for windows of GINFINITY
node embeddings:

1. fit spherical k-means centroids over the node embeddings;
2. persist the complete `k x 128` centroid matrix as float16;
3. persist the `k x k` float32 centroid-similarity matrix;
4. replace each node with a centroid code and form fixed-width code windows;
5. search those code windows with a custom C++ hnswlib distance;
6. optionally rerank the candidate pool with the original window vectors; and
7. compare the result with exact FlatIP on `queries_rouskin_6k.tsv`.

The non-negotiable data contract is that quantization is candidate selection
only. The original residue embeddings remain full 128-dimensional float16
arrays and are still the inputs to GINFINITY-SW, alignment, and plotting. The
new representation is additive:

```text
original embeddings (float16, n_nodes x 128)
        |
        +--> original windows (float32 search view; float16 source preserved)
        |       `faiss/embeddings.npz` remains the SW/alignment source
        |
        +--> centroid assignment --> uint16 node codes
                                  --> uint16 code windows
                                  --> compact custom-distance HNSW candidates
                                  --> optional exact original-vector rerank
```

The optional GPU branch is a separate candidate representation:

```text
original windows (normalized float32 search view)
        |
        +--> global non-clipping int8 scale
                --> cuVS CAGRA candidates on GPU
                --> exact rerank with original residue embeddings
```

It does not alter the CPU custom-distance index or replace the preserved
float16 node embeddings. `HNSWLIB_GPU_CAGRA` is the metadata name for this
companion index; it is deliberately explicit because cuVS CAGRA does not
provide hnswlib's custom `SpaceInterface`.

## C++ custom distance

The Python hnswlib binding exposes only the built-in `l2`, `ip`, and `cosine`
spaces; it does not expose a Python callback for a new distance. The official
C++ API does expose `SpaceInterface`, `DISTFUNC`, and a raw data-size hook. This
implementation uses that API rather than expanding every code window into a
`w * 128` float vector.

The pinned vendor bundle is:

```text
vendor/hnswlib-0.8.0/
  hnswlib/*.h                 # upstream header-only implementation
  compact_hnswlib.cpp         # GINflow build/search driver and custom space
  LICENSE                     # upstream Apache-2.0 license
```

`compact_hnswlib.cpp` registers a `CodeSpace` with:

```cpp
data_size = window_size * sizeof(uint16_t)
distance(A, B) = -sum(S[A[p], B[p]] for p=0..window_size-1)
```

hnswlib ranks smaller distances first, hence the negative sign. `S` is the
persisted centroid-similarity matrix. At index build time hnswlib stores the
raw uint16 code window; at search time it receives another raw code window and
calls the registered C++ function. The index does not store expanded centroid
vectors.

There are therefore four distinct stored artifacts:

| Artifact | Stored data | Purpose |
|---|---|---|
| `quantization/centroids.npy` | `k x 128`, float16 | Assign query nodes and retain the fitted codebook |
| `quantization/similarity.npy` | `k x k`, float32 | Lookup table used by the custom distance |
| `faiss/index.bin` | HNSW graph plus `w` uint16 codes per element | Compact candidate selection |
| `faiss/embeddings.npz` | original arrays, float16, each `n x 128` | SW, alignment, and plots |

For the default window length 11, the custom index data payload is only 22
bytes per window. The 128-dimensional centroids are side data for assignment
and audit; they are not substituted for the original node embeddings and are
not copied into each HNSW element.

## Metric definition

Node vectors are L2-normalized before spherical k-means. Centroids are
renormalized and written as float16. Similarity is calculated in float32 from
the stored float16 centroids:

```text
S[i, j] = dot(float32(centroid[i]), float32(centroid[j]))
```

For code windows `A` and `B` of length `w`:

```text
code_score(A, B) = sum(S[A[p], B[p]] for p=0..w-1)
```

Without reranking, the exported seed score is this raw sum. Since the public
`seed_min_similarity` parameter is a per-position cosine-like threshold, the
compact search converts it to `seed_min_similarity * window_size` internally.

With `--hnswlib_rerank true`, the HNSW code search only chooses candidates. Each
candidate is then scored with the original normalized window formed from the
preserved float16 residue embeddings. The exported score and threshold return
to the normal full-window cosine convention. This is the recommended mode when
recall matters. The production reranker groups candidates by source embedding
array and scores each group in vectorized batches; short hnswlib result rows
are represented as missing candidates and do not abort the search.

## Pipeline shape

The Nextflow path is:

```text
BUILD_RNA_GRAPHS
        |
EMBED_RNA_GRAPHS ----------------------------+
        |                                      |
        |                               original embeddings
        |                                      |
        +--> FIT_NODE_QUANTIZER                +--> GENERATE_WINDOWS
                  |                                  |
                  +--> centroids + S + node codes   +--> original windows
                                      |
                         GENERATE_QUANTIZED_WINDOWS
                                      |
                         BUILD_HNSWLIB_INDEX
                                      |
                         SEARCH_HNSWLIB
                                      |
                         CLUSTER_SEEDS -> ALIGN_CLUSTERS -> DRAW_SW
```

For a query-only run, query nodes are assigned to the database's persisted
centroids. Query fitting is never performed independently. The search process
receives both the quantized query windows and all original query embedding
shards when reranking is enabled.

When `--hnswlib_gpu true` is selected, the build and search processes also
stage the original window arrays and their manifests. The GPU search quantizes
those normalized window arrays to int8 only for CAGRA candidate traversal and
then uses the original embedding shards for exact scoring.

## Parameters and recommended profiles

The requested defaults remain conservative and compatible with the existing
pipeline:

| Parameter | Default | Role |
|---|---:|---|
| `--node_quantization_k` | `2048` | Requested spherical k-means centroid count |
| `--node_quantization_sample_size` | `500000` | Maximum node vectors used for fitting |
| `--node_quantization_niter` | `25` | Spherical k-means iterations |
| `--node_quantization_seed` | `1` | Reproducible centroid fitting |
| `--hnswlib_m` | `32` | HNSW graph degree |
| `--hnswlib_ef_construction` | `200` | Build exploration |
| `--hnswlib_ef_search` | `100` | Search exploration when no override is supplied |
| `--hnswlib_random_seed` | `1` | HNSW construction seed |
| `--hnswlib_num_threads` | `0` | C++ OpenMP worker setting; zero uses the default |
| `--hnswlib_candidate_k` | `50` | Code-HNSW candidates retrieved before output/rerank |
| `--hnswlib_rerank` | `false` | Score the candidate pool with original float16 windows |

GPU companion defaults:

| Parameter | Default | Role |
|---|---:|---|
| `--hnswlib_gpu` | `false` | Select cuVS CAGRA for HNSWLIB candidate search; requires `-profile gpu`. |
| `--hnswlib_gpu_candidate_k` | `50` | Candidates retained before exact reranking. |
| `--hnswlib_gpu_itopk_size` | `256` | CAGRA intermediate beam; must be at least candidate_k. |
| `--hnswlib_gpu_search_batch_size` | `512` | Query windows per GPU batch. |
| `--hnswlib_gpu_intermediate_graph_degree` | `128` | CAGRA build-time degree. |
| `--hnswlib_gpu_graph_degree` | `64` | Retained CAGRA degree. |
| `--hnswlib_gpu_build_algo` | `nn_descent` | CAGRA build algorithm in the pinned cuVS release. |
| `--hnswlib_gpu_int8_scale` | auto | Optional global non-clipping scale for original window vectors. |

The high-recall Rouskin profile is:

```bash
--index hnswlib \
--node_quantization_k 4096 \
--hnswlib_m 32 \
--hnswlib_ef_construction 200 \
--hnswlib_ef_search 5000 \
--hnswlib_candidate_k 5000 \
--hnswlib_rerank true
```

The 2048-centroid value remains the requested default, but it is not the
highest-recall setting measured here. Increasing centroid count improves the
quantization fidelity; increasing `M` or `ef` alone does not recover the
information lost by a small codebook. The recommended profile is intentionally
explicit rather than silently changing the user-requested default.

Small inputs use an effective centroid count no larger than the available
training sample. Compact uint16 codes support up to 65,536 centroids.

## Published artifacts and temporary cache

Published database artifacts under `faiss/` are:

```text
faiss/index.bin
faiss/windows.tsv
faiss/meta.json
faiss/quantization/centroids.npy
faiss/quantization/similarity.npy
faiss/quantization/quantization.json
faiss/embeddings.npz
faiss/records.tsv
```

The workflow also publishes the inspectable derived artifacts at the output
root:

```text
quantization/
windows_quantized/
```

The latter are not inputs to SW. The database metadata records
`index_type=HNSWLIB_COMPACT`, `candidate_representation=uint16_centroid_codes`,
`index_data_dtype=uint16`, the code bytes per element, centroid dimensions,
HNSW settings, and the original embedding dtype/preservation flag.

For the GPU companion, the database additionally contains:

```text
faiss/cagra/index.bin             # cuVS CAGRA graph and serialized int8 dataset
faiss/meta.json                   # HNSWLIB_GPU_CAGRA, scale, graph, and batch settings
```

The GPU metadata reports
`candidate_representation=original_normalized_window_int8`,
`gpu_metric=sqeuclidean`, and `original_embeddings_preserved=true`. The
copied node quantization artifacts remain available for the CPU path and
inspection; they are not used as the GPU CAGRA metric.

Large, disposable research artifacts live under:

```text
/mnt/ssd_samsung/ginflow-hnsw-research/
```

The important subdirectories used in this investigation are:

```text
env-benchmark/                         # Python 3.12, NumPy, FAISS, hnswlib
compact_hnswlib                        # compiled C++ driver
benchmark-compact-rouskin-6k-v2/       # compact no-rerank k/M/ef sweep
benchmark-compact-highk/                # k=8192/16384 and M=32/64
benchmark-compact-rerank-k4096/        # candidate_k=1000
benchmark-compact-rerank-k4096-c5000/  # candidate_k=5000
benchmark-compact-rerank-final/        # fresh paired final measurement
benchmark-compact-rerank-thresholds/   # fresh paired threshold sweep
benchmark-compact-candidate-sweep/     # candidate_k=1000/2000/3000/5000
benchmark-compact-ef-rerank/           # ef-search sweep with reranking
benchmark-compact-k-default-rerank/    # k=2048 versus k=4096 rerank comparison
benchmark-compact-rerank-w7/           # window-size comparison
benchmark-compact-rerank-w15/
real-compact-rerank-vectorized-smoke/  # real end-to-end smoke output
real-compact-rerank-query/             # real query-only output
cagra-original-int8-scale850-cuvs2410/ # full Rouskin CAGRA source benchmark
cagra-to-hnsw-full/                    # converted index, CPU labels, metrics
cagra-to-hnsw-smoke/                   # small conversion smoke check
faiss-hnsw-int8-full/                  # FAISS IndexHNSWSQ, IP, full Rouskin
faiss-hnsw-flat-int8-full/             # float-backed int8-coordinate control
faiss-hnsw-int8-l2-full/               # FAISS IndexHNSWSQ, L2 comparison
faiss-hnsw-int8-ip-k100/               # loaded-index candidate-width sweep
faiss-hnsw-int8-ip-k200/               # loaded-index candidate-width sweep
faiss-hnsw-int8-ip-k500/               # loaded-index candidate-width sweep
faiss-hnsw-int8-ip-k1000/              # loaded-index candidate-width sweep
```

Do not commit these generated directories. The benchmark scripts are kept in
`bin/` because they are useful for reproducing the measurements; the large
arrays and indexes remain disposable.

## Benchmark method

The benchmark uses:

- `tests/data/rouskin_sample_6k.tsv`: 5,840 records, 897,588 nodes;
- `tests/data/queries_rouskin_6k.tsv`: 512 query rows;
- 839,188 database windows at `window_size=11`, `stride=1`; and
- original node embeddings cast to the pipeline's stored float16 dtype before
  forming normalized FlatIP windows.

The query TSV contains stale `vector_id` and `record_ordinal` values from a
different window ordering. Stable `transcript_id + window_offset` resolves all
512 rows and is the authoritative query identity. Only 3 rows match the
current `record_ordinal`, so the benchmark does not use that field to select
queries.

For a pipeline-native query run, convert those selections into a structures TSV
with `bin/convert_query_selections.py`. The generated
`tests/data/queries_rouskin_6k_structures.tsv` contains the source molecule
sequence and structure plus one 0-based half-open `start`/`end` slice per
selection. It is directly ingestible by `--query`; the original selection TSV
itself is not, because it has no sequence or secondary-structure columns.

If the intended experiment is ordinary full-molecule querying rather than
reproducing the 512 selected benchmark windows, pass `--full-molecules`. This
creates `tests/data/queries_rouskin_6k_full_structures.tsv` with one
deduplicated full molecule per selected transcript. The pipeline then searches
all sliding windows from each query molecule, so its output is not directly
comparable to the one-window-per-row benchmark recall.

```bash
python3 bin/convert_query_selections.py \
  --structures tests/data/rouskin_sample_6k.tsv \
  --selections tests/data/queries_rouskin_6k.tsv \
  --output tests/data/queries_rouskin_6k_structures.tsv \
  --window-size 11
```

For full-molecule query records:

```bash
python3 bin/convert_query_selections.py \
  --structures tests/data/rouskin_sample_6k.tsv \
  --selections tests/data/queries_rouskin_6k.tsv \
  --output tests/data/queries_rouskin_6k_full_structures.tsv \
  --full-molecules
```

The exact reference is `faiss.IndexFlatIP` over the original full windows. The
compact result is measured twice: first as raw code-HNSW candidates, then after
optional exact reranking against the same original-window vectors. This
separates graph recall from quantization/rerank recall.

The GPU experiment uses the pinned cuVS 24.10.00/CuPy 13.0.0 image and the
same normalized original windows. cuVS 24.10's CAGRA path supports the
sqeuclidean metric for signed int8 data, so the benchmark uses a global scale
of 850.0. The exact original-window matrix remains the reference; the int8
matrix is never used for final scores.

Reproduction command:

```bash
env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV -u CONDA_SHLVL \
  PYTHONNOUSERSITE=1 \
  /mnt/ssd_samsung/ginflow-hnsw-research/env-benchmark/bin/python \
  bin/benchmark_compact_hnswlib.py \
  --embeddings /home/nicolas/programs/GINFINITY/experiments/rouskin_sample_6k_quantization/embeddings.float32.npy \
  --structures tests/data/rouskin_sample_6k.tsv \
  --queries tests/data/queries_rouskin_6k.tsv \
  --centroids-dir /home/nicolas/programs/GINFINITY/experiments/rouskin_sample_6k_quantizer_comparison \
  --outdir /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds \
  --executable /mnt/ssd_samsung/ginflow-hnsw-research/compact_hnswlib \
  --k-values 4096 --m-values 32 --ef-search-values 5000 \
  --ef-construction 200 --num-threads 8 --batch-size 32768 \
  --top-k 50 --candidate-k 5000 --rerank --rerank-batch-size 4 \
  --thresholds 0.70,0.75,0.80,0.85,0.90,0.95
```

The C++ driver can be rebuilt from the repository with:

```bash
g++ -O3 -std=c++11 -fopenmp \
  -I vendor/hnswlib-0.8.0 \
  vendor/hnswlib-0.8.0/compact_hnswlib.cpp \
  -o /mnt/ssd_samsung/ginflow-hnsw-research/compact_hnswlib
```

The candidate-pool experiment reuses the same harness with
`--outdir .../benchmark-compact-candidate-sweep`,
`--candidate-k-values 1000,2000,3000,5000`, and
`--ef-search-values 5000`. The ef experiment uses
`--outdir .../benchmark-compact-ef-rerank`,
`--candidate-k 5000`, and `--ef-search-values 1000,3000,5000`.

The standalone GPU sweep helper is `bin/benchmark_cagra_gpu.py`. The Rouskin
measurement used the pinned cuVS 24.10 image and was reproduced with:

```bash
docker run --rm --gpus all --ipc=host --ulimit memlock=-1 \
  -v "$PWD":/work -v /mnt/ssd_samsung:/mnt/ssd_samsung -w /work \
  community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4 \
  python3 bin/benchmark_cagra_gpu.py \
    --data /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flat_windows_float16_w11_s1.npy \
    --queries /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/gpu_original_queries.float32.npy \
    --reference /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flatip_reference_labels.int64.npy \
    --outdir /mnt/ssd_samsung/ginflow-hnsw-research/cagra-original-int8-scale850-cuvs2410 \
    --dtype int8 --int8-scale 850 --k 50 --itopk-size 256 \
    --intermediate-graph-degree 128 --graph-degree 64 \
    --search-batch-size 512
```

The production adapter is `bin/hnswlib_gpu.py`; it adds manifest validation,
database packaging, query batching, exact reranking, seed thresholding, and
the `HNSWLIB_GPU_CAGRA` metadata contract around the same CAGRA representation.

The CAGRA-to-CPU-HNSW experiment is research-only and does not change the
production Nextflow path. The two C++ drivers use the cuVS C API because the
pinned 24.10 Python binding rejects int8 in `hnsw.from_cagra`, and its Python
search wrapper expects a CAGRA object even after conversion. The C API can
serialize the int8 CAGRA graph to cuVS's HNSW wrapper and search that wrapper
from host memory.

The generated HNSW file is not compatible with upstream hnswlib despite the
historical `hnswlib` name in cuVS. It is a cuVS-specific immutable base-layer
format. `bin/search_cuvs_hnsw_cpu.cpp` loads that file without loading the
CAGRA graph and accepts host `int8` query windows.

The runtime image is deliberately kept unchanged. For this one-off research
build, a compiler and development headers were installed in a disposable
container. The exact source files are:

```text
bin/cagra_to_hnsw_cpu.cpp       # CAGRA load, serialize, VRAM sampling
bin/search_cuvs_hnsw_cpu.cpp    # HNSW-only CPU search
bin/benchmark_cagra_hnsw_cpu.py # FlatIP recall and original-window rerank
```

For the pinned runtime image used in this investigation, the temporary build
command was equivalent to:

```bash
docker run --rm --gpus all \
  -v "$PWD":/work \
  -v /usr/local/cuda-12.9/targets/x86_64-linux/include:/host-cuda:ro \
  -w /work \
  community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4 \
  bash -lc '
    apt-get update -qq && apt-get install -y -qq g++ >/dev/null
    mkdir -p /tmp/dlpack/dlpack
    ln -sf /opt/conda/lib/python3.12/site-packages/cupy/_core/include/cupy/_dlpack/dlpack.h /tmp/dlpack/dlpack/dlpack.h
    for source in bin/cagra_to_hnsw_cpu.cpp bin/search_cuvs_hnsw_cpu.cpp; do
      output=/work/bin/$(basename "$source" .cpp)
      g++ -std=c++17 -O3 -I/opt/conda/include -I/tmp/dlpack -I/host-cuda \
        -pthread "$source" -L/opt/conda/lib -Wl,-rpath,/opt/conda/lib \
        -lcuvs_c -lcuvs -lraft -lcudart -o "$output"
    done
  '
```

The installed compiler is disposable; it is not part of the production GPU
environment or the Nextflow module.

After compiling both C++ drivers with the cuVS libraries, quantize the 512
normalized query windows with the same scale as the CAGRA database:

```bash
python3 - <<'PY'
import numpy as np

source = np.load(
    "/mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/"
    "gpu_original_queries.float32.npy"
)
np.rint(source * 850.0).astype(np.int8).tofile(
    "/mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/queries.int8"
)
PY
```

Convert the full Rouskin CAGRA index once:

```bash
docker run --rm --gpus all \
  -v "$PWD":/work -v /mnt/ssd_samsung:/mnt/ssd_samsung -w /work \
  community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4 \
  /work/bin/cagra_to_hnsw_cpu \
    --cagra-index /mnt/ssd_samsung/ginflow-hnsw-research/cagra-original-int8-scale850-cuvs2410/cagra.index \
    --hnsw-index /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/hnsw.index \
    --queries-int8 /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/queries.int8 \
    --n-queries 512 --dimension 1408 --k 50 --ef-values 50,100,200,400,800 \
    --num-threads 8 \
    --labels-prefix /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/labels- \
    --distances-prefix /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/distances- \
    --metrics /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/metrics.json
```

Then run the persisted HNSW file in a container with no GPU device. The
resulting labels were byte-for-byte identical to the labels from the
GPU-enabled conversion process:

```bash
docker run --rm -e CUDA_VISIBLE_DEVICES='' \
  -v "$PWD":/work -v /mnt/ssd_samsung:/mnt/ssd_samsung -w /work \
  community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4 \
  /work/bin/search_cuvs_hnsw_cpu \
    --hnsw-index /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/hnsw.index \
    --queries-int8 /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/queries.int8 \
    --n-queries 512 --dimension 1408 --k 50 --ef-values 50,100,200,400,800 \
    --num-threads 16 \
    --labels-prefix /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/nogpu-labels- \
    --distances-prefix /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/nogpu-distances- \
    --metrics /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/nogpu-metrics.json
```

Evaluate candidate recall and exact rerank against the FlatIP labels:

```bash
python3 bin/benchmark_cagra_hnsw_cpu.py \
  --database-windows /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flat_windows_float16_w11_s1.npy \
  --queries /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/gpu_original_queries.float32.npy \
  --reference-labels /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flatip_reference_labels.int64.npy \
  --labels-prefix /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/nogpu-labels- \
  --search-metrics /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/nogpu-metrics.json \
  --candidate-k 50 --ef-search-values 50,100,200,400,800 \
  --output /mnt/ssd_samsung/ginflow-hnsw-research/cagra-to-hnsw-full/nogpu-evaluation.json
```

## Results

### Exact FlatIP reference versus the compact rerank profile

The fresh paired threshold run is `benchmark-compact-rerank-thresholds`:

| Metric | FlatIP | Compact code-HNSW + original-vector rerank |
|---|---:|---:|
| Database windows | 839,188 | 839,188 |
| Query windows | 512 | 512 |
| Recall@1 | reference | 0.9746 |
| Recall@5 | reference | 0.9852 |
| Recall@10 | reference | 0.9848 |
| Recall@50 | reference | **0.9755** |
| Self-hit@1 | 0.9727 | 0.9785 |
| Index bytes | 4,726,306,861 | 250,282,332 |
| Build seconds | 4.59 | 27.29 |
| Query seconds | 24.67 | 0.79 C++ search + 3.01 rerank = 3.80 |

The exact FlatIP index is approximately 4.73 GB; the compact HNSW index is
approximately 250 MB, about 18.9 times smaller. The custom index stores 22
bytes of code payload per w=11 window; graph links account for most of the
remaining index size. Independent fresh runs with the same profile obtained
R@50 values between `0.9724` and `0.9780`; every fresh c=5000 run cleared 97%.
Parallel HNSW insertion is not bit-for-bit deterministic even with a fixed
seed, which explains the small run-to-run variation.

The approximate query timing includes the original-vector rerank but excludes
one-time query code packing and Python process startup. FlatIP timing likewise
measures the FAISS search call after its query matrix is available. These are
index-search timings, not total Nextflow wall time.

### Requested default k versus high-recall k

The requested default `k=2048` remains usable, but the Rouskin data shows why
the recommended high-recall profile uses `k=4096`. Both rows use the same
M=32, ef-construction=200, candidate_k=5000, and original-vector rerank:

| Centroids | Recall@1 | Recall@5 | Recall@10 | Recall@50 | Threshold-0.8 micro-recall |
|---:|---:|---:|---:|---:|---:|
| 2,048 | 0.9766 | 0.9816 | 0.9760 | 0.9629 | 0.9705 |
| 4,096 | 0.9727 | 0.9844 | 0.9844 | **0.9766** | **0.9768** |

This is the evidence behind keeping 2048 as the compatibility default while
documenting 4096 as the optimized profile.

### Compact code-only sweep

These are raw code-HNSW results with `top_k=50` and no original-vector rerank.
They show the quantization ceiling clearly:

| k | M | ef search | Recall@1 | Recall@50 | Self-hit@1 | Index bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 1,024 | 32 | 100 | 0.8477 | 0.3672 | 0.8516 | 250,282,200 |
| 1,024 | 32 | 400 | 0.8535 | 0.3682 | 0.8574 | 250,282,200 |
| 2,048 | 32 | 100 | 0.8242 | 0.4124 | 0.8281 | 250,282,332 |
| 2,048 | 32 | 400 | 0.8301 | 0.4138 | 0.8340 | 250,282,332 |
| 4,096 | 32 | 100 | 0.8457 | 0.4699 | 0.8496 | 250,282,200 |
| 4,096 | 32 | 400 | 0.8457 | 0.4714 | 0.8496 | 250,282,200 |
| 8,192 | 32 | 1,000 | 0.8613 | 0.4971 | 0.8652 | 250,282,200 |
| 8,192 | 64 | 200 | 0.8613 | 0.4983 | 0.8652 | 465,082,216 |
| 16,384 | 32 | 1,000 | 0.8711 | 0.5421 | 0.8750 | 250,282,200 |
| 16,384 | 64 | 200 | 0.8711 | 0.5440 | 0.8750 | 465,082,216 |

Increasing `k` is the dominant code-only improvement. M=64 approximately
doubles graph-link storage and did not materially improve recall. Increasing
`ef_search` from 200 to 1,000 helped only modestly. Code-only HNSW should
therefore be treated as a compact candidate generator, not as the final
high-recall ranking method.

### Candidate pool, thresholds, and window size

At `k=4096`, `M=32`, `ef_construction=200`, `ef_search=5000`, and
`window_size=11`, the fresh candidate-pool sweep was:

| Candidate pool | ef search | R@1 | R@5 | R@10 | R@50 |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 5,000 | 0.9668 | 0.9590 | 0.9457 | 0.9020 |
| 2,000 | 5,000 | 0.9668 | 0.9723 | 0.9645 | 0.9438 |
| 3,000 | 5,000 | 0.9668 | 0.9793 | 0.9758 | 0.9608 |
| 5,000 | 5,000 | 0.9707 | 0.9801 | 0.9797 | **0.9724** |

The candidate-only R@50 was only 0.4738 even for the 5,000 pool; the >97%
result comes from scoring those candidates with the original embeddings. The
5,000 pool is the smallest tested pool that clears the target. It also means
that the C++ driver must use an effective `ef >= 5000`, because hnswlib cannot
return 5,000 neighbors with a smaller exploration queue; requested ef values
below `candidate_k` are raised to `candidate_k` internally.

The reranked ef sweep at candidate_k=5000 produced effectively identical
rankings for requested ef=1000, 3000, and 5000 (R@50 0.9780 in that run),
because all three searches use the effective ef=5000. This makes the candidate
pool, rather than a separately tunable lower ef, the dominant cost/recall
control at the recommended profile.

Threshold analysis compares the exact FlatIP top-50 set with the emitted
reranked top-50 set in original-window cosine units. `thresholds.json` records
the complete table; the representative fresh run was:

| Threshold | Exact hits/query | Approx hits/query | Micro recall | Micro precision |
|---:|---:|---:|---:|---:|
| 0.70 | 24.516 | 24.490 | 0.9830 | 0.9840 |
| 0.75 | 18.697 | 18.693 | 0.9789 | 0.9791 |
| 0.80 | 14.510 | 14.510 | 0.9732 | 0.9732 |
| 0.85 | 10.598 | 10.598 | 0.9633 | 0.9633 |
| 0.90 | 6.729 | 6.729 | 0.9422 | 0.9422 |
| 0.95 | 3.758 | 3.758 | 0.8966 | 0.8966 |

At the pipeline default threshold 0.8, this fresh run emitted exactly the same
mean number of seeds per query as FlatIP. Independent graph builds can change
the count by a few thousandths of a hit per query. Thresholds are applied after
reranking; the code-only path uses the equivalent centroid raw sum divided by
window size for this analysis.

Window-size exploration used the same k/M/ef/candidate settings. The w=15
run has 484 valid query rows because 28 selected offsets do not fit in a
15-residue window; it is not a direct 512-row comparison.

| Window size | Valid queries | Windows | Recall@1 | Recall@50 |
|---:|---:|---:|---:|---:|
| 7 | 512 | 862,548 | 0.9668 | 0.9602 |
| 11 | 512 | 839,188 | 0.9746 | 0.9755 |
| 15 | 484 | 815,828 | 0.9814 | 0.9873 |

The current public default remains w=11 because it preserves the existing
query-window contract and provides a complete comparison set. w=15 is a
promising dataset-specific tuning direction, not a silent global default
change.

### GPU CAGRA companion versus FlatIP and compact CPU HNSW

The GPU benchmark uses 839,188 database windows, 512 query windows, `k=50`,
`itopk_size=256`, CAGRA degrees 128/64, and `nn_descent`. The database is the
same normalized original-window matrix used by the exact FlatIP reference.

| Path | Candidate representation | Query/search seconds | R@1 | R@5 | R@10 | R@50 |
|---|---|---:|---:|---:|---:|---:|
| Exact FlatIP | float32 normalized windows | 24.67 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| Compact CPU HNSW + exact rerank | uint16 centroid-code windows | 0.79 candidate + 3.01 rerank = 3.80 | 0.9746 | 0.9852 | 0.9848 | 0.9755 |
| GPU CAGRA candidates | original windows, int8 scale 850 | **0.0298 CAGRA search** | 0.9688 | 0.9813 | 0.9820 | 0.9900 |
| GPU CAGRA + exact rerank | same candidates, original embeddings for score | ~0.105 scoring/search kernel | **0.9727** | **0.9875** | **0.9898** | **0.9900** |

The CAGRA search kernel is approximately 26.5x faster than the recorded
0.79-second compact CPU candidate traversal, and the GPU search-plus-rerank
kernel is about 36x faster than the 3.80-second compact CPU search/rerank
measurement. A serialized GPU index load took 0.813 seconds in the same
benchmark, so a cold per-process end-to-end query is about 0.92 seconds before
pipeline startup and I/O; the more-than-10x claim applies to the query search
kernel after the index is loaded. The GPU CAGRA build took 39.0 seconds and
the serialized index was about 1.4 GiB, so this is a query-throughput option,
not a smaller-index option.

The exact rerank is essential: raw int8 CAGRA labels are only approximate
candidates, while final scores and thresholding use the original embeddings.
At the default candidate count of 50, the measured GPU final R@50 is already
above the CPU compact result. Increase `--hnswlib_gpu_candidate_k` and
`--hnswlib_gpu_itopk_size` together if a new dataset shows lower recall; lower
`--hnswlib_gpu_search_batch_size` when VRAM is constrained.

### CAGRA converted to CPU HNSW

This experiment answers a different question from the GPU CAGRA benchmark:
can a built CAGRA graph be converted once and then searched on the CPU? The
answer is yes for the pinned cuVS 24.10 build, through its C API. The CPU
HNSW path does not use centroid codes or the custom centroid-similarity
function. It searches the original normalized window vectors after the same
global `int8` scale of 850.0 used by CAGRA, with squared L2 as the candidate
metric. Exact reranking still uses the original normalized float32 view of the
preserved float16 node embeddings.

The full run used 839,188 windows and 512 queries. Recall is measured against
the exact FlatIP top-50 labels. The CPU row uses the HNSW-only process with no
GPU device exposed, `ef_search=200`, and 8 search threads; the exact rerank is
included in the reported total.

| Path | Stored index bytes | One-time build/convert | Search + exact rerank for 512 queries | R@1 | R@5 | R@10 | R@50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Exact FlatIP | 4,726,306,861 | 4.59 s build | 24.672 s | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| CAGRA GPU | 1,396,410,545 | 39.037 s build | 0.0297 s GPU search | 0.9688 | 0.9813 | 0.9820 | 0.9900 |
| Converted cuVS CPU HNSW, candidate_k=50 | 1,409,835,936 | 15.034 s conversion | 0.210 + 0.022 = **0.233 s** | 0.9688 | 0.9871 | 0.9898 | 0.9895 |

The converted file is about 1.31 GiB, only slightly larger than the 1.30 GiB
CAGRA file. Its one-time HNSW reload took 0.728 s. The CPU search is about
105x faster than the recorded FlatIP query call at the candidate_k=50 profile,
while losing about 1.05 percentage points of R@50. Increasing the candidate
width lets exact reranking recover more of the FlatIP set:

| Candidate width | ef_search | CPU search | Exact rerank | Total | Final R@50 |
|---:|---:|---:|---:|---:|---:|
| 50 | 200 | 0.210 s | 0.022 s | 0.233 s | 0.9895 |
| 100 | 400 | 0.333 s | 0.097 s | 0.430 s | 0.9941 |
| 200 | 800 | 0.559 s | 0.131 s | 0.690 s | 0.9943 |
| 500 | 2,000 | 1.094 s | 0.297 s | 1.390 s | 0.9955 |
| 1,000 | 4,000 | 1.857 s | 0.554 s | 2.411 s | **0.9968** |

The best speed/recall point among these measurements is dataset- and target-
dependent. Candidate width 100 already reaches 99.41% R@50 in 0.43 s, while
width 1,000 approaches the exact FlatIP overlap at roughly 2.4 s. Search
thread testing at `ef_search=200`, candidate_k=50, gave 0.703/0.386/0.237/
0.196/0.183 s for 1/2/4/8/16 threads; all returned identical labels, so 16
threads was used for the wider candidate experiments.

The observed VRAM measurements during CAGRA load and conversion were:

| Measurement | Bytes | MiB |
|---|---:|---:|
| Baseline before CAGRA load | 158,466,048 | 151.1 |
| Peak total | 1,557,266,432 | 1,485.1 |
| Peak increase | 1,398,800,384 | **1,334.0** |
| After destroying CAGRA index | 158,466,048 | 151.1 |

This confirms the intended detach behavior experimentally: after the HNSW
serialization completes, destroying the CAGRA object returns the observed GPU
allocation to baseline. The HNSW search-only process then runs with no GPU
device and produces the same labels as the GPU-enabled process. The peak is a
sampled measurement at 10 ms intervals, so it should be treated as an observed
requirement rather than a formal allocator upper bound.

This converted CPU format is not the custom centroid-code HNSW index already
used by the `--index hnswlib` CPU path. It is a separate cuVS wrapper around
the CAGRA graph, and it is not readable by upstream Python hnswlib. The
production pipeline still uses CAGRA directly for `--hnswlib_gpu`; this
experiment establishes a possible persisted CPU fallback/serving format but
does not silently switch the production backend.

### FAISS HNSW over the same original-window int8 representation

The next experiment compares the cuVS path with a conventional FAISS HNSW
index. It uses the same data and the same candidate representation as the
CAGRA benchmark: the original float16 node embeddings are first assembled
into normalized float32 windows, then rounded with a global scale of 850.0 to
the signed-int8 range. The original float16-derived windows are retained for
the exact rerank and are never replaced by the int8 values for final scoring.
This is window quantization for candidate selection; it is not the centroid
code representation used by the compact custom-distance HNSWLIB path.

FAISS 1.10.0 provides `IndexHNSWSQ`, so the true compact test uses
`QT_8bit_direct_signed`. The Python API still accepts float32 input arrays;
the scalar-quantizer index is what stores one byte per coordinate internally.
For an apples-to-apples quality control, `IndexHNSWFlat` was also given the
same integer-valued coordinates. That control has essentially the same search
math but stores float32 values, so its larger index is expected. Both indexes
used inner product, `M=32`, `efConstruction=200`, and 16 OpenMP threads. The
exact rerank uses the original normalized windows and the reference is the
same FlatIP top-50 label set used by the CAGRA measurements.

The full build/search command was:

```bash
docker run --rm --ipc=host \
  -v "$PWD":/work -v /mnt/ssd_samsung:/mnt/ssd_samsung -w /work \
  community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e \
  python3 /work/bin/benchmark_faiss_hnsw_int8.py \
    --database-windows /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flat_windows_float16_w11_s1.npy \
    --queries /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/gpu_original_queries.float32.npy \
    --reference-labels /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flatip_reference_labels.int64.npy \
    --outdir /mnt/ssd_samsung/ginflow-hnsw-research/faiss-hnsw-int8-full \
    --variant sq-int8 --metric ip --int8-scale 850 \
    --ef-search-values 50,100,200,400,800 --candidate-k 50 --output-k 50 \
    --m 32 --ef-construction 200 --num-threads 16
```

The `flat-int8` control uses the same command with
`--variant flat-int8` and a different output directory. The L2 comparison
uses `--variant sq-int8 --metric l2`. Existing indexes can be reused for
candidate-width experiments with `--load-index`; this avoids rebuilding the
839,188-vector graph.

At candidate_k=50, the IP results were:

| Path | Index bytes | Peak RSS (MiB) | Build | Search + rerank | Final R@1 | Final R@5 | Final R@10 | Final R@50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Exact FlatIP reference | 4,726,306,861 | — | 4.59 s | 24.672 s | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| FAISS `IndexHNSWSQ`, int8, IP, ef=200 | 1,409,946,230 | 6,259 | 735.742 s | **0.381 s** | 0.9668 | 0.9789 | 0.9793 | 0.9762 |
| FAISS `IndexHNSWFlat`, int8-valued float32, IP, ef=200 | 4,954,676,306 | 9,569 | 523.783 s | **0.373 s** | 0.9688 | 0.9813 | 0.9818 | 0.9782 |
| Converted cuVS CPU HNSW, int8, ef=200 | 1,409,835,936 | — | 15.034 s conversion* | **0.233 s** | 0.9688 | 0.9871 | 0.9898 | 0.9895 |

The cuVS conversion row is not a fresh graph build: the CAGRA graph was built
on the GPU in 39.037 s and then serialized to the CPU HNSW format in 15.034 s.
Thus its build-plus-conversion wall time was about 54.1 s, whereas FAISS built
its HNSW graph from scratch on the CPU. The timings were collected on
different code paths and are useful engineering measurements, not a
controlled hardware-normalized benchmark. The cuVS CPU row used 8 search
threads; the FAISS rows used 16.

The important storage result is that FAISS `IndexHNSWSQ` and converted cuVS
HNSW are almost identical in size: 110,294 bytes apart for this dataset. The
float-backed control is 3.51 times larger than the true scalar-quantized
FAISS index and used about 3.23 GiB more peak resident memory, while its
recall differed by less than 0.2 percentage points at R@50. This confirms
that the direct-signed scalar-quantizer representation is the memory-saving
part; the graph topology is the dominant remaining cost.

The FAISS L2 variant was also measured with the same true int8 storage:

| Metric | Build | Index bytes | Peak RSS | Search + rerank at ef=200 | Final R@50 at ef=200 | Final R@50 at ef=800 |
|---|---:|---:|---:|---:|---:|---:|
| IP | 735.742 s | 1,409,946,230 | 6,259 MiB | 0.381 s | 0.9762 | 0.9798 |
| L2 | 964.652 s | 1,409,946,230 | 6,255 MiB | 0.441 s | 0.9664 | 0.9672 |

IP is the better FAISS candidate metric for this FlatIP reference. The L2
variant is not equivalent after int8 rounding because quantization changes the
per-window norms; the exact rerank cannot recover candidates that L2 omitted.

Widening the FAISS IP candidate pool improves final overlap, at the cost of
both HNSW traversal and original-vector reranking:

| Candidate pool | ef search | Search + rerank | Final R@50 |
|---:|---:|---:|---:|
| 50 | 200 | 0.381 s | 0.9762 |
| 100 | 800 | 1.222 s | 0.9893 |
| 200 | 800 | 1.369 s | 0.9916 |
| 500 | 2,000 | 4.562 s | 0.9936 |
| 1,000 | 4,000 | 16.635 s | **0.9937** |

The converted cuVS CPU path reached 0.9941 at candidate_k=100, 0.9955 at
candidate_k=500, and 0.9968 at candidate_k=1,000 in its own measurements.
Therefore, for this Rouskin workload, FAISS HNSW is a viable CPU int8
candidate index with approximately the same compact storage as the converted
cuVS format, but the cuVS-built graph delivered higher recall and faster
search at comparable candidate widths. The experiment does not change the
production index selector: FAISS HNSW remains a separate research baseline,
while the current pipeline's compact centroid-code HNSW and optional CAGRA
paths retain their existing contracts.

### Node-level product quantization with custom C++ hnswlib

This experiment quantizes the original 128-dimensional node embeddings before
forming candidate windows. It is different from both the centroid-code path
and the FAISS/CAGRA int8-window path:

```text
original float16 node vectors
        |
        +--> product quantizer over each 128-d node
                |
                +--> packed M-subquantizer code per node
                        |
                        +--> packed code window in custom C++ hnswlib
                        +--> exact rerank with original float16-derived windows
```

For `M` subquantizers and `nbits` bits per subquantizer, each subvector has
`128/M` dimensions and `2**nbits` centroids. The HNSW payload stores only
the packed subcodes. The custom distance looks up the dot product between
matching subquantizer centroids at each window position and sums those values.
The PQ codebook and lookup table are auxiliary artifacts; the original
float16 node embeddings are never replaced and remain the source for exact
reranking and downstream SW/alignment.

The full grid used 839,188 database windows, 512 query windows, `w=11`,
897,588 nodes, a 250,000-node training sample, 20 PQ iterations, `M=32`,
`efConstruction=200`, 16 threads, candidate pools 50/100, and `efSearch`
200/800. Final scores were always computed from the original normalized
windows. Peak RSS below is the observed child-process peak during search; the
index size includes the HNSW graph and packed window codes.

```bash
docker run --rm -e LD_LIBRARY_PATH=/opt/conda/lib --ipc=host \
  -v "$PWD":/work -v /mnt/ssd_samsung:/mnt/ssd_samsung \
  -v /home/nicolas/programs/GINFINITY:/home/nicolas/programs/GINFINITY:ro \
  -w /work \
  community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e \
  python3 /work/bin/benchmark_pq_hnswlib.py \
    --embeddings /home/nicolas/programs/GINFINITY/experiments/rouskin_sample_6k_quantization/embeddings.float32.npy \
    --structures /work/tests/data/rouskin_sample_6k.tsv \
    --query-selections /work/tests/data/queries_rouskin_6k.tsv \
    --database-windows /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flat_windows_float16_w11_s1.npy \
    --queries /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/gpu_original_queries.float32.npy \
    --reference-labels /mnt/ssd_samsung/ginflow-hnsw-research/benchmark-compact-rerank-thresholds/flatip_reference_labels.int64.npy \
    --outdir /mnt/ssd_samsung/ginflow-hnsw-research/pq-hnsw-rouskin \
    --results-name results.json \
    --executable /mnt/ssd_samsung/ginflow-hnsw-research/pq_hnswlib \
    --pq-m-values 4,8,16 --nbits-values 4,8 \
    --candidate-k-values 50,100 --ef-search-values 200,800 \
    --hnsw-m-values 32 --ef-construction-values 200 \
    --threads 16 --rerank-batch-size 32
```

At `candidate_k=100`, `efSearch=800`, the quantizer/packing trade-off was:

| PQ layout | Node code | Window code (`w=11`) | Index bytes | Search RSS | Candidate R@50 | Final R@50 | Search + rerank |
|---|---:|---:|---:|---:|---:|---:|---:|
| `M=4, 4-bit` | 2 B | 22 B | 250.3 MB | 348.6 MiB | 0.2844 | 0.3777 | 0.613 s |
| `M=4, 8-bit` | 4 B | 44 B | 268.7 MB | 367.1 MiB | 0.4318 | 0.5569 | 0.698 s |
| `M=8, 4-bit` | 4 B | 44 B | 268.7 MB | 366.2 MiB | 0.4681 | 0.6012 | 0.875 s |
| `M=8, 8-bit` | 8 B | 88 B | 305.7 MB | 403.5 MiB | 0.5823 | 0.7262 | 0.718 s |
| `M=16, 4-bit` | 8 B | 88 B | 305.7 MB | 401.5 MiB | 0.6364 | 0.7991 | 0.859 s |
| `M=16, 8-bit` | 16 B | 176 B | 379.5 MB | 476.1 MiB | 0.7424 | 0.8945 | 0.897 s |

The graph dominates storage: doubling the code payload does not double the
index size. More importantly, the candidate pool is a stronger recall control
than `efSearch` once the PQ layout is fixed. The 50/100-pool grid is therefore
only the compact first pass, not the final high-recall setting.

For `M=16`, 8-bit, HNSW `M=32`, the wider candidate sweep used
`efSearch=1,000/2,000/4,000` and candidate pools 50/100/200/500/1,000. The
table below selects `efSearch=4,000` so the pool sizes can be compared at the
same graph exploration setting:

| Candidate pool | Search + rerank | Search RSS | Candidate R@50 | Final R@50 |
|---:|---:|---:|---:|---:|
| 50 | 1.580 s | 477.5 MiB | 0.7472 | 0.7472 |
| 100 | 1.669 s | 477.6 MiB | 0.7475 | 0.9015 |
| 200 | 1.726 s | 478.1 MiB | 0.7477 | 0.9597 |
| 500 | 1.941 s | 479.9 MiB | 0.7489 | 0.9821 |
| 1,000 | 2.236 s | 482.9 MiB | 0.7491 | **0.9889** |

The candidate R@50 remains around 0.75 because it measures the approximate
PQ ranking before exact reranking. The final R@50 rises as the exact reranker
gets a wider pool. Missing hnswlib neighbors are represented as sentinels and
are ignored by the reranker; they do not become database ID zero.

The graph-degree sweep used `M=16/32/64`, `efConstruction=200`,
`candidate_k=1,000`, and `efSearch=4,000` with the same `M=16`, 8-bit PQ
layout. It measures whether more graph links can improve traversal recall:

| HNSW M | Index bytes | Build | Build RSS | Search + rerank | Search RSS | Final R@50 |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 272.3 MB | 185.9 s | 514 MiB | 2.390 s | 380.9 MiB | 0.9821 |
| 32 | 379.5 MB | cached; 261.9 s uncached | 616 MiB uncached | 2.295 s | 482.9 MiB | 0.9889 |
| 64 | 594.3 MB | 465.9 s | 821 MiB | 2.565 s | 687.8 MiB | **0.9941** |

`M=64` is the highest-recall PQ point measured here, but its extra 215 MB of
index storage and roughly 2.5x build time over `M=16` buy only 1.2 percentage
points of R@50 over `M=16`. `M=32` is the balanced choice; `M=64` is justified
when approximately 99.4% overlap is worth the memory and build cost. The
benchmark does not claim that these settings are universal outside the
Rouskin workload.

Compared with the other measured paths, the best node-PQ result is:

| Path | Index bytes | Search + exact rerank | Final R@50 |
|---|---:|---:|---:|
| Exact FlatIP | 4,726 MB | 24.67 s | 1.0000 |
| PQ `M=16`, 8-bit, HNSW `M=32`, pool 1,000 | 379.5 MB | 2.24–2.30 s | 0.9889 |
| PQ `M=16`, 8-bit, HNSW `M=64`, pool 1,000 | 594.3 MB | 2.56 s | **0.9941** |
| Converted cuVS CPU HNSW int8, pool 1,000 | 1,409.8 MB | 2.41 s | 0.9968 |
| GPU CAGRA int8, pool 50 | 1,396.4 MB | 0.0298 s GPU candidate search | 0.9900 |

The PQ path is substantially smaller than the int8 window indexes and is
competitive in CPU query time. Its remaining gap to exact FlatIP comes from
the candidate-generation representation; exact reranking cannot recover a
window that never enters the pool. All detailed machine-readable PQ results
are kept in the SSD cache under
`/mnt/ssd_samsung/ginflow-hnsw-research/pq-hnsw-rouskin/`, including
`results.json`, `results-wide.json`, and `results-hnsw-m.json`.

## Validation completed

The following checks have passed in the current working tree:

- targeted host unit tests: 7 tests passed across the GPU representation and
  compact C++ HNSW helpers, with optional legacy Python-hnswlib tests skipped
  when hnswlib is not installed;
- GPU representation tests cover automatic/non-clipping int8 scaling, separate
  window-manifest staging, manifest coordinates, and original-row handling;
- FAISS HNSW int8 benchmark helper tests cover scale/rounding, clipping
  rejection, original-window reranking, recall, and index reuse arguments;
- `tests/test_parameter_values.py` and JSON schema validation;
- `nf-test test tests/search.nf.test --tag hnswlib --profile +docker`;
- CPU and GPU stub tests for both `BUILD_HNSWLIB_INDEX` and `SEARCH_HNSWLIB`;
- `nf-test test tests/hnswlib_gpu.nf.test --profile +docker` for GPU-profile
  validation;
- pinned cuVS 24.10 real GPU build/search smoke with separate window and
  manifest staging, persisted-index reload, and exact original-embedding
  reranking;
- cuVS C API CAGRA-to-HNSW conversion smoke and full 839,188-window Rouskin
  conversion, including sampled VRAM, persisted HNSW reload, CPU-only search,
  and byte-for-byte comparison of GPU-enabled versus no-GPU labels;
- candidate-width and CPU-thread sweeps for the converted HNSW file, with
  FlatIP recall and exact reranking evaluated by `bin/benchmark_cagra_hnsw_cpu.py`;
- real smoke build/search with compact C++ HNSW, original float16 preservation,
  exact rerank, clustering, alignment, and report generation;
- real query-only search against the published compact database with persisted
  centroids, exact rerank, clustering, alignment, and report generation;
- real compact smoke search with `--plot_sw true`: two SW SVGs were produced
  from a database whose packed residue arrays remained float16;
- vectorized production reranking on the real compact smoke database, including
  short-result handling when a small smoke graph returns fewer candidates than
  requested;
- conversion of all 512 Rouskin selection rows into a sliced query structures
  table, followed by a real 512-query Nextflow run through search, clustering,
  alignment, and report generation;
- reduced-parameter real 1,200-record test-profile database build: 437,316
  windows, compact index, and original float16 embeddings retained; and
- full-parameter real 1,200-record test-profile database build with
  `k=2048`, `M=32`, `ef_construction=200`, and `ef_search=100`: 437,316
  windows, a 130,404,480-byte compact index, and original float16 embeddings
  retained; and
- the compact C++ standalone test, which verifies returned distances against
  an explicit Python implementation of the registered matrix sum and checks
  that the index payload is smaller than dense vector storage.

The complete multi-backend command in `AGENTS.md` passed all 33 tests in 393.5
seconds with `--profile +docker`. The HNSW changes therefore have regression
coverage alongside the existing FAISS/ScaNN/NGT/cuVS modules.

## Reproducibility checksums

```text
tests/data/rouskin_sample_6k.tsv
  7aebc1af275c33803ce23c98840e5f467c62bdf04566064b1e14a11be129aac9
tests/data/queries_rouskin_6k.tsv
  1bbf3c1ff814abe253b678688a310c32c62fa1d52e2fbc31cf1744899483301e
GINFINITY embeddings.float32.npy
  7cfdd2c0860fce95d3baf6909fbc9f94a91d2467c70abdf967a85206908783d7
```

The Rouskin codebooks came from
`/home/nicolas/programs/GINFINITY/experiments/rouskin_sample_6k_quantizer_comparison/`
and the benchmark JSON/CSV files under the SSD cache record the detailed
machine-readable values.

## Plan and TODO checklist

### Implementation

- [x] Add centroid fitting and persisted float16 centroids.
- [x] Add float32 `k x k` centroid similarity generation.
- [x] Add quantized node-code and quantized-window modules.
- [x] Add C++ hnswlib custom `DISTFUNC` over raw uint16 code windows.
- [x] Vendor/pin hnswlib 0.8.0 headers and license.
- [x] Add compact database build/search and query-only centroid application.
- [x] Keep original embeddings and records in every built database.
- [x] Add optional original-vector reranking after candidate selection.
- [x] Add `--index hnswlib`/`hnsw` dispatch, metadata, and unused-parameter
  warnings.
- [x] Publish quantization and code-window artifacts for inspection.
- [x] Add the optional cuVS CAGRA GPU companion with int8 candidate windows
  and exact original-embedding reranking.
- [x] Add a research-only node-level PQ encoder, packed-window builder, and
  custom C++ hnswlib distance over positional subquantizer lookup tables.
- [x] Preserve and use original float16-derived windows for exact PQ reranking,
  including short-result/sentinel handling.

### Tests

- [x] Unit-test centroid determinism, effective k, code ranges, window slicing,
  similarity sums, and original dtype preservation.
- [x] Unit-test C++ custom distances and exact rerank ordering.
- [x] Run `tests/test_pq_hnswlib.py`: packed-code round trips, node-major
  window layout, positional C++ distance sums, and missing-neighbor reranking.
- [x] Test hnswlib Python legacy compatibility in the benchmark environment.
- [x] Test HNSWLIB dispatch with nf-test.
- [x] Run real smoke build/search and real query-only rerank.
- [x] Run a reduced real test-profile build.
- [x] Run the complete documented nf-test command across all backends and
  profiles: 33 tests passed in 393.5 seconds with `--profile +docker`.

### Optimization and comparison

- [x] Compare exact FlatIP with compact code-only HNSW.
- [x] Sweep k from 1,024 through 16,384.
- [x] Sweep M=32 and M=64 and multiple ef-search values.
- [x] Measure candidate pools of 1,000, 2,000, 3,000, and 5,000 with original
  reranking.
- [x] Sweep seed thresholds from 0.70 through 0.95 and record hit counts,
  overlap, precision, and recall.
- [x] Measure window sizes 7, 11, and 15, recording the w=15 query-count caveat.
- [x] Establish a >97% R@50 profile on the 6k query set.
- [x] Compare exact/approximate seed counts, overlap, precision, and recall
  against the FlatIP baseline across thresholds 0.70 through 0.95.
- [x] Benchmark the GPU CAGRA companion against FlatIP and compact CPU HNSW on
  the Rouskin 6k window set.
- [x] Convert the full CAGRA graph to cuVS's CPU HNSW format through the C API
  and measure conversion time, serialized size, CPU search, and recall.
- [x] Verify the converted HNSW file in a no-GPU process and compare its labels
  with the GPU-enabled conversion process.
- [x] Sweep converted-HNSW candidate width, `ef_search`, and CPU thread count
  against the FlatIP reference.
- [x] Benchmark FAISS `IndexHNSWSQ` with direct-signed int8 windows against a
  float-backed `IndexHNSWFlat` control.
- [x] Compare FAISS IP and L2 candidate metrics, index size, build time, peak
  host memory, search time, candidate width, and exact-rerank recall.
- [x] Sweep node-PQ `M`/`nbits`, candidate pools 50 through 1,000, ef-search,
  and HNSW graph degree 16/32/64 on the full Rouskin query set.
- [x] Add CPU fallback and GPU stub/unit coverage for both HNSWLIB modules.
- [ ] Run a larger downstream SW/alignment diff against the FlatIP baseline at
  the production threshold.
- [x] Decide the production policy: retain the requested k=2048 default and
  document the explicit k=4096/candidate_k=5000 rerank profile for high recall.
- [ ] Measure total Nextflow wall time and peak memory, not only index-call
  timings.
- [ ] Repeat the GPU/FlatIP comparison through the complete Nextflow Rouskin
  path, including clustering, alignment, and report generation.

### Cleanup and maintenance

- [x] Run `repowise update` after the final meaningful edits.
- [ ] Keep only the benchmark scripts needed for reproducibility; move stable
  experiment code to a dedicated `benchmarks/` directory if desired.
- [ ] Remove failed/duplicate SSD caches after the result JSON/CSV files are
  archived. No repository files depend on those caches.
- [x] Recheck `git diff --check`, schema validation, Python compilation, and
  the complete regression suite before merging.

## Experiment log

| Date | Experiment | Result |
|---|---|---|
| 2026-08-19 | Python binding and upstream hnswlib review | Python exposes only built-in spaces; C++ exposes custom `SpaceInterface`/`DISTFUNC`. |
| 2026-08-19 | Dense centroid-expanded Python HNSW sweep | Best tested codebook-only result was R@50 0.4736 at k=4096; storage was about 4.96 GB. |
| 2026-08-19 | Compact C++ HNSW implementation | Raw code payload is `w*2` bytes; standalone score/index smoke passed. |
| 2026-08-19 | Compact code-only k/M/ef sweep | k=16,384/M=64 reached R@50 0.5440; increasing graph parameters alone was insufficient. |
| 2026-08-20 | k=4096, candidate_k=1000, original rerank | R@50 0.9041. |
| 2026-08-20 | k=4096, candidate_k=5000, original rerank | Fresh paired runs ranged from R@50 0.9724 to 0.9780; the threshold-artifact run was 0.9755. |
| 2026-08-20 | Candidate-pool sweep | R@50 0.9020/0.9438/0.9608/0.9724 for candidate_k 1,000/2,000/3,000/5,000. |
| 2026-08-20 | Threshold sweep | At threshold 0.8: 14.510 exact and approximate hits/query, micro-recall and precision 0.9732. |
| 2026-08-20 | Window-size sweep | R@50 0.9602 (w=7), 0.9755 (w=11), 0.9873 (w=15 with 484 queries). |
| 2026-08-20 | Vectorized rerank hardening | Grouped original-window scoring and short-result handling passed unit tests and a real Nextflow smoke run. |
| 2026-08-20 | Pipeline-native Rouskin queries | Converted 512 selection rows to sliced structures and completed a real query-only Nextflow run through reporting. |
| 2026-08-20 | Real compact smoke/query-only runs | C++ build, custom search, original-vector rerank, downstream alignment, and reports passed. |
| 2026-08-20 | GPU CAGRA companion benchmark | cuVS 24.10 int8 CAGRA search took 0.0298 s for 512 queries and reached R@50 0.9900 after exact reranking; serialized load was 0.813 s. |
| 2026-08-20 | GPU HNSWLIB module integration | Original window arrays/manifests, GPU CAGRA metadata, exact rerank, CPU fallback, and GPU stubs/unit tests passed. |
| 2026-08-20 | CAGRA-to-cuVS-CPU-HNSW conversion | Full int8 graph serialized in 15.034 s to a 1,409,835,936-byte HNSW file; observed peak VRAM increase was 1,334 MiB and returned to baseline after CAGRA destruction. |
| 2026-08-20 | Converted HNSW CPU-only search | No-GPU search produced identical labels; candidate_k=50 reached R@50 0.9895 in 0.233 s including exact rerank, and candidate_k=1,000 reached 0.9968 in 2.411 s. |
| 2026-08-20 | FAISS HNSWSQ true-int8 IP benchmark | `IndexHNSWSQ(QT_8bit_direct_signed)` built in 735.742 s, occupied 1,409,946,230 bytes, and reached R@50 0.9762 at ef=200/candidate_k=50; ef=800 reached 0.9798. |
| 2026-08-20 | FAISS float-backed int8 control | `IndexHNSWFlat` over the same integer-valued coordinates built in 523.783 s and occupied 4,954,676,306 bytes; recall was within 0.2 percentage points of true int8. |
| 2026-08-20 | FAISS HNSWSQ L2 and candidate-width sweep | L2 reached R@50 0.9664 at ef=200; IP reached 0.9893/0.9916/0.9936/0.9937 at candidate pools 100/200/500/1,000 with the listed wider ef values. |
| 2026-08-20 | Node-level PQ custom-distance sweep | `M=16`, 8-bit PQ with HNSW `M=32`, candidate_k=1,000, and exact rerank reached R@50 0.9889 in 2.24 s with a 379.5 MB index; the 50/100 first-pass grid ranged from 0.3777 to 0.8945 R@50. |
| 2026-08-20 | Node-level PQ HNSW-degree sweep | At candidate_k=1,000, HNSW `M=16/32/64` reached R@50 0.9821/0.9889/0.9941 with 272.3/379.5/594.3 MB indexes; `M=64` peaked at 821 MiB during its 465.9 s build. |
