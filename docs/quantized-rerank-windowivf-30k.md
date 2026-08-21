# Batched exact reranking, 30k benchmark, and WindowIVF research artifact

Status: implementation, smoke tests, the full 30k PQ-HNSW matrix, and the
first WindowIVF sweep are complete. Large generated files belong under
`/mnt/ssd_samsung/ginflow-hnsw-research`, not in git.

See the [headline results table](quantized-hnsw-research.md#headline-results)
for a consolidated comparison with the earlier 6k HNSW, CAGRA, FAISS, and
node-PQ experiments.

## Goal

Measure a two-stage search path on the larger Rouskin sample:

```text
node embeddings (original float16, 128 values/node)
  ├─► node PQ / centroid candidate representation
  │      └─► custom-distance HNSW or standalone WindowIVF
  └─► original normalized windows
         └─► exact batched reranking of only the candidate labels
```

The original node arrays remain authoritative. Quantization is used only to
choose candidates. Exact reranking reconstructs each selected window from the
original residue embeddings and computes the normal dot product of normalized
concatenated windows. SW, alignment, plots, and the rest of the pipeline keep
using the original embeddings.

The benchmark reports candidate recall and final recall against a blockwise
FlatIP-equivalent reference, with emphasis on R@100 and R@500. Candidate pools
of 1,000, 5,000, and 10,000 are tested even when the final pipeline seed count
is 50 or 500.

## Pipeline change: `RERANK_CANDIDATES`

`modules/rerank_candidates/main.nf` invokes `bin/rerank_candidates.py` after
`SEARCH_HNSWLIB` when `--exact_rerank true` (or the compatibility alias
`--hnswlib_rerank true`) is set. It accepts the candidate seed table, the
published database, and the query window artifacts.

The reranker has two bounded loops:

1. `--exact_rerank_batch_size` controls how many query windows are processed
   together.
2. `--exact_rerank_candidate_batch_size` controls how many candidate columns
   are materialized at once.

CPU mode uses multiple workers around NumPy kernels. CUDA mode uses CuPy for
the batched dot products and one GPU worker. Both modes return the same stable
score ordering. `rerank_metrics.json` records the backend, batching, worker
count, candidate width, and elapsed time.

The exact kernel is shared by the pipeline and research benchmarks in
`bin/rerank_core.py`; this avoids measuring a different implementation from
the one used in production.

`--hnswlib_candidate_k` (or `--hnswlib_gpu_candidate_k`) is the preselection
width. `--seed_k` is the number retained after exact scoring and thresholding.
For a recall experiment with R@500, use `--seed_k 500` or use the standalone
benchmark’s `--output-k 500`; a 50-row seed table cannot measure a 500-neighbour
final recall.

## Query convention

`tests/data/queries_rouskin_structures.tsv` is the shared, ordinary full-
molecule query input for both 6k and 30k pipeline runs. It has 512 molecule
rows and is directly ingestible by `--query`.

ANN recall needs a fixed window coordinate. The previous benchmark’s 512
`transcript_id`/`window_offset` pairs were recovered into:

```text
/mnt/ssd_samsung/ginflow-hnsw-research/rouskin_shared_queries/query_selections.tsv
```

The IDs were checked against the new shared structure file. This external
coordinate map does not replace the pipeline query file and is deliberately
kept out of the repository. If it is deleted, recreate it from immutable git
history before cleanup or make a new deterministic selection map and record
that choice in the benchmark metadata.

## Resumable artifact layout

```text
/mnt/ssd_samsung/ginflow-hnsw-research/
├── rouskin-30k-pipeline-sq/           # published DB: windows, embeddings, records
├── rouskin-shared-query-pipeline/     # published query windows
├── rouskin-30k-inputs/                # packed arrays and exact FlatIP-equivalent reference
├── rouskin-30k-overnight/             # PQ-HNSW/CAGRA logs and JSON results
├── tooling/pq_hnswlib                 # compiled custom HNSW driver
└── tooling/window_ivf                 # compiled standalone WindowIVF driver
```

The directory labels are shown without spaces; the actual names are exactly
the paths above. Each configuration writes its own codebook, node codes,
window-code cache, index, raw labels, timing file, and result row. Re-running
the command uses those files when their metadata signatures match.

## Prepare the 30k and shared-query artifacts

The preparation uses a low-memory FAISS SQ index as a convenient published
container for the original embeddings/windows. The packer then builds a
separate blockwise exact FlatIP-equivalent reference from those original
windows; FlatIP is not the production candidate selector. The blockwise
implementation avoids a second 23 GiB in-memory copy while preserving the
same inner-product ranking. It is intentionally recorded as
`FlatIP-equivalent` rather than pretending that a serialized FAISS
`IndexFlatIP` was built for this 30k run.

```bash
cd /home/nicolas/programs/ginflow

nextflow run main.nf -profile docker -resume -c conf/rouskin_30k_research.config \
  --input tests/data/rouskin_sample_30k.tsv \
  --outdir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-pipeline-sq \
  --index faiss --faiss_index sq --faiss_sq_type 8bit --shard_size 400 \
  --ginfinity_full_precision false --window_size 11 --window_stride 1 \
  --plot_backend none

nextflow run main.nf -profile docker -resume \
  --input tests/data/queries_rouskin_structures.tsv \
  --outdir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-shared-query-pipeline \
  --index faiss --faiss_index flatip --shard_size 64 \
  --ginfinity_full_precision false --window_size 11 --window_stride 1 \
  --plot_backend none
```

Then pack the published artifacts and build the exact reference:

```bash
docker run --rm \
  -v /home/nicolas/programs/ginflow:/work \
  -v /mnt/ssd_samsung:/mnt/ssd_samsung \
  -w /work \
  community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e \
  python3 bin/pack_research_inputs.py \
  --database /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-pipeline-sq/faiss \
  --structures tests/data/rouskin_sample_30k.tsv \
  --query-structures tests/data/queries_rouskin_structures.tsv \
  --query-selections /mnt/ssd_samsung/ginflow-hnsw-research/rouskin_shared_queries/query_selections.tsv \
  --query-windows-dir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-shared-query-pipeline/windows \
  --query-manifests-dir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-shared-query-pipeline/windows \
  --outdir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs \
  --original-dtype float16 --reference-k 500 --reference-batch-size 32768
```

Check `inputs.json` before starting a long sweep. It must report the actual
database window count, 512 queries, dimension 1408, and `original_dtype`.

## PQ-HNSW overnight matrix

Compile the custom driver once:

```bash
g++ -O3 -std=c++17 -fopenmp -static-libgcc -static-libstdc++ \
  -I vendor/hnswlib-0.8.0 vendor/hnswlib-0.8.0/pq_hnswlib.cpp \
  -Wl,-Bstatic -lgomp -Wl,-Bdynamic \
  -o /mnt/ssd_samsung/ginflow-hnsw-research/tooling/pq_hnswlib
```

The full matrix considers:

| Dimension | Values | Why |
|---|---:|---|
| node-PQ M | 8, 16 | 16 subquantizers use more codebook detail |
| node-PQ bits | 4, 8 | 0.5/1 byte per node for M=16 |
| HNSW M | 32, 64 | graph degree/recall trade-off |
| ef-search | 1,000; 5,000; 10,000 | search effort |
| candidate pool | 1,000; 5,000; 10,000 | exact rerank recall ceiling |
| final output | 500 | permits R@500 |

Launch it with:

```bash
python3 bin/run_rouskin_overnight.py \
  --embeddings /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/node_embeddings.float16.npy \
  --structures tests/data/rouskin_sample_30k.tsv \
  --query-selections /mnt/ssd_samsung/ginflow-hnsw-research/rouskin_shared_queries/query_selections.tsv \
  --database-windows /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/database_windows_float16_w11_s1_source_order.npy \
  --queries /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/query_windows.float32.npy \
  --reference-labels /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/reference_labels.npy \
  --pq-executable /mnt/ssd_samsung/ginflow-hnsw-research/tooling/pq_hnswlib \
  --outdir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-overnight \
  --pq-m-values 8,16 --nbits-values 4,8 \
  --hnsw-m-values 32,64 \
  --ef-search-values 1000,5000,10000 \
  --candidate-k-values 1000,5000,10000 \
  --ef-construction-values 40 --output-k 500 --rerank-workers 8 \
  --rerank-batch-size 8 --candidate-batch-size 1024 \
  --sample-size 500000 --pq-niter 25 --hnsw-threads 16 \
  --original-dtype float16
```

For a short validation run, reduce to `--pq-m-values 16 --nbits-values 4
--hnsw-m-values 32 --ef-search-values 5000 --candidate-k-values 1000`.
`results.json` contains `candidate_recall_at_100`,
`candidate_recall_at_500`, `final_recall_at_100`, and
`final_recall_at_500` for every row, plus build/search/rerank timing and peak
RSS.

## GPU CAGRA experiment

The existing CAGRA companion accepts `--dtype int8` and now supports the same
batched exact reranker. On an 8 GiB RTX 3070 laptop GPU, a 30k database with
about 4.1 million 1408-dimensional windows needs approximately 5.40 GiB
(5.79 decimal GB) for the raw int8 matrix alone (`4,115,576 * 1,408`
bytes), before graph workspace,
temporary buffers, and allocator overhead. It may not fit. The driver should
be run as a separate, resumable attempt and an OOM is an informative resource
result, not a reason to discard the CPU benchmark.

```bash
python3 bin/run_rouskin_overnight.py \
  ...same input arguments as above... \
  --cagra \
  --cagra-index /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-cagra-int8.index
```

The GPU branch tries candidate pools 1,000/5,000/10,000, reranks to 500 with
CuPy, and records candidate/final R@100/R@500. If 30k does not fit, keep the
6k CAGRA result as the GPU comparison and report the 30k memory failure.

The first 30k attempt used int8 scale 127, inner-product metric, and
`k=1000`. Upload reached about 5.7 GiB on the 8 GiB RTX 3070, then CAGRA
failed during graph construction with CUDA `out_of_memory`; no serialized
index was produced. This is an observed capacity result, not an inference
from the raw-matrix estimate.

## Standalone WindowIVF

`bin/window_ivf.cpp` is the first custom IVF implementation. It deliberately
does not depend on FAISS internals yet. It stores sampled coarse
representatives as packed node-PQ windows, list offsets, reordered packed
window codes, and original window labels.

Coarse routing, probing, and scanning all call the same positional lookup
function as the HNSW path. `--nprobe` controls the number of lists and the
output is still a candidate pool for exact reranking. This makes it a clear
correctness/performance baseline before adopting FAISS `InvertedLists` and
scanner classes.

Compile and smoke-test it with:

```bash
g++ -O3 -std=c++17 -fopenmp -static-libgcc -static-libstdc++ \
  bin/window_ivf.cpp -Wl,-Bstatic -lgomp -Wl,-Bdynamic \
  -o /mnt/ssd_samsung/ginflow-hnsw-research/tooling/window_ivf
python3 -m unittest -v tests/test_window_ivf.py
```

The first 30k WindowIVF experiment is intentionally modest because the
prototype’s coarse assignment is brute-force:

```bash
python3 bin/benchmark_window_ivf.py \
  --embeddings /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/node_embeddings.float16.npy \
  --structures tests/data/rouskin_sample_30k.tsv \
  --query-selections /mnt/ssd_samsung/ginflow-hnsw-research/rouskin_shared_queries/query_selections.tsv \
  --database-windows /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/database_windows_float16_w11_s1_source_order.npy \
  --queries /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/query_windows.float32.npy \
  --reference-labels /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-inputs/reference_labels.npy \
  --executable /mnt/ssd_samsung/ginflow-hnsw-research/tooling/window_ivf \
  --outdir /mnt/ssd_samsung/ginflow-hnsw-research/rouskin-30k-window-ivf \
  --pq-m-values 16 --nbits-values 4 \
  --nlist-values 16,32 --nprobe-values 1,4,16 \
  --candidate-k-values 1000,5000,10000 --output-k 500
```

Do not interpret this initial nlist sweep as the final IVF design. The
production follow-up is to use FAISS `InvertedLists` storage and scanner
lifecycle while replacing the distance computation with the positional lookup
function.

## Observed 30k results

The completed PQ-HNSW run contains 72 rows: 8 combinations of node-PQ/HNSW
parameters times 3 search efforts times 3 candidate pools. It used all 512
shared queries, 4,115,576 database windows, exact scoring from the original
float16 node embeddings, and `output_k=500`.

The strongest rows by candidate pool were:

| PQ M/bits | HNSW M | effective ef | pool | final R@100 | final R@500 | search + exact rerank | index |
|---|---:|---:|---:|---:|---:|---:|---:|
| 16/8 | 64 | 10,000 | 1,000 | 0.9674 | 0.9005 | 5.51 s | 2.91 GB |
| 16/4 | 64 | 10,000 | 5,000 | 0.9827 | 0.9599 | 9.41 s | 2.55 GB |
| 16/4 | 64 | 10,000 | 10,000 | 0.9855 | 0.9792 | 14.05 s | 2.55 GB |

The best overall final recall was the last row. The 16/4 node PQ also beat
16/8 in the highest-recall configuration here, so the extra 8-bit code
precision did not justify its larger code and index for this dataset. A
smaller M32 graph reached R@100=0.9538/R@500=0.9291 at a 5,000 pool in 8.19 s
with a 1.50 GB index, making it a reasonable lower-memory operating point.
The exact rerank portion of the 10,000-pool rows was about 9.4--10.3 s for
all 512 queries; the candidate search was about 4--5 s.

The exact reference search took 29.19 s for the 512 queries using the
blockwise `IndexFlatIP`-equivalent implementation. The reference labels are
therefore directly comparable in ranking semantics, but the benchmark does
not claim that a full serialized FAISS FlatIP index fit through the 30k build
memory peak.

The standalone WindowIVF prototype used 16/4 node PQ. Representative rows
were:

| nlist | nprobe | pool | final R@100 | final R@500 | search + exact rerank | index | build peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 4 | 10,000 | 0.9440 | 0.9273 | 42.77 s | 395 MB | 756 MiB |
| 16 | 16 | 10,000 | 0.9903 | 0.9857 | 133.71 s | 395 MB | 756 MiB |
| 32 | 4 | 10,000 | 0.9194 | 0.8944 | 26.16 s | 395 MB | 756 MiB |
| 32 | 16 | 10,000 | 0.9884 | 0.9821 | 72.99 s | 395 MB | 756 MiB |

WindowIVF uses much less index storage than PQ-HNSW because it stores packed
window codes and lists instead of a large HNSW graph. Its current scanner is
still a straightforward custom-distance scan, however, so latency is much
higher at the probe counts needed for high recall. This prototype is useful
as a correctness and storage baseline, not yet as the production index.

## TODO and acceptance criteria

- [x] Keep original 128-dimensional float16 node arrays in the database.
- [x] Separate exact reranking from candidate search as a named Nextflow
  process.
- [x] Support bounded CPU multi-worker and GPU batch reranking.
- [x] Add process, unit, custom C++ smoke, and HNSW regression tests.
- [x] Use the same full-molecule query table for 6k and 30k pipeline runs.
- [x] Complete and archive the 30k PQ-HNSW matrix with R@100/R@500.
- [x] Compare final rerank recall and timing against the exact FlatIP-equivalent reference.
- [x] Record whether int8 CAGRA fits in the available GPU memory.
- [x] If CAGRA fits, compare GPU search and GPU rerank timings; if not, retain
  the explicit OOM/peak-memory result.
- [x] Complete the first 30k WindowIVF nlist/nprobe experiment.
- [ ] Replace brute-force coarse assignment with a faster coarse-training
  strategy if WindowIVF is competitive.
- [ ] Prototype FAISS `InvertedLists`/scanner integration behind a separate
  extension target; do not couple it to the production pipeline until recall
  and serialization tests pass.
- [x] Add a full Nextflow integration test with `--exact_rerank true` and
  assert that `rerank_metrics.json` is published.
- [ ] Archive only JSON/CSV summaries and delete raw 30k code/index caches
  after results are backed up.

## Cleanup map

Safe-to-recreate artifacts are everything below the research directory except
the final JSON/CSV summaries and the external query coordinate map. Before
deleting them, preserve:

```text
inputs.json
*/results.json
*/flatip.json
*/run-complete.json
query_selections.tsv
```

The repository additions are the reranker, benchmark drivers, C++ prototype,
tests, and this document. No production result depends on a hidden cache.
