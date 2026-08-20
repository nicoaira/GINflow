# Quantized-node HNSW research plan, implementation, and results

Status: implemented and validated on the smoke pipeline; the Rouskin 6k
benchmark exceeds 97% FlatIP recall with the recommended candidate-rerank
profile. The complete documented nf-test regression suite and the full default
k=2048 test-corpus build have passed. Optional larger downstream comparisons
and whole-pipeline resource measurements remain follow-up work.

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

## Validation completed

The following checks have passed in the current working tree:

- host unit tests: 10 tests passed, with the two legacy Python-hnswlib tests
  skipped when hnswlib is not installed;
- benchmark Conda environment unit tests: all 10 passed, including the legacy
  Python round-trip tests and the vendored C++ custom-distance tests;
- `tests/test_parameter_values.py` and JSON schema validation;
- `nf-test test tests/search.nf.test --tag hnswlib --profile +docker`;
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

### Tests

- [x] Unit-test centroid determinism, effective k, code ranges, window slicing,
  similarity sums, and original dtype preservation.
- [x] Unit-test C++ custom distances and exact rerank ordering.
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
- [ ] Run a larger downstream SW/alignment diff against the FlatIP baseline at
  the production threshold.
- [x] Decide the production policy: retain the requested k=2048 default and
  document the explicit k=4096/candidate_k=5000 rerank profile for high recall.
- [ ] Measure total Nextflow wall time and peak memory, not only index-call
  timings.

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
