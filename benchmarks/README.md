# Reproducible vector-index benchmark cache

This directory separates the expensive, shared data preparation from individual
FAISS, ScaNN, NGT, and cuVS benchmark runs. Generated corpora stay outside the
Git worktree; the committed code records exactly how they were made.

## Cache once, then benchmark many backends

Set an external cache root with sufficient space. The overnight run uses the
large SSD rather than the repository filesystem:

```bash
export GINFLOW_BENCHMARK_CACHE=/mnt/ssd_samsung/ginflow-benchmarks
```

Materialize graph, residue-embedding, and window artifacts exactly once for
each requested source table. The command uses the pipeline's existing
`PREPARE_WINDOWS` workflow with a persistent Nextflow work directory and
`-resume`; it stops before any index is built.

```bash
python3 benchmarks/prepare_cache.py cache-windows \
  --dataset rouskin_6k \
  --input tests/data/rouskin_sample_6k.tsv

python3 benchmarks/prepare_cache.py cache-windows \
  --dataset rouskin_30k \
  --input tests/data/rouskin_sample_30k.tsv
```

The command refuses to reuse a cache whose source checksum, window parameters,
or deterministic preparation implementation fingerprint differs. That
fingerprint hashes the benchmark entry workflow, `PREPARE_WINDOWS`, the three
runtime module definitions and environments, `bin/slice_windows.py`, and the
applicable Nextflow configuration/profile files. It uses live file contents,
not only a Git revision, so uncommitted implementation changes cannot silently
reuse stale windows. Source paths and timestamps remain provenance only: the
same source bytes can be relocated without changing the request identity.

Each dataset cache is locked while Nextflow runs. Its
`window-cache-request.json` is written as `in_progress` before execution and
becomes `complete` only after matching graph, embedding, and window shard sets
are published. An interrupted same-identity command can resume with Nextflow
`-resume`; a mismatched identity fails without mixing artifacts. The command
never removes an existing cache.

After windows exist, create the backend-neutral inputs:

```bash
python3 benchmarks/prepare_cache.py flatten --dataset rouskin_6k
python3 benchmarks/prepare_cache.py select-queries --dataset rouskin_6k --query-count 512
python3 benchmarks/prepare_cache.py ground-truth --dataset rouskin_6k --engine auto
```

Repeat the same three commands for `rouskin_30k`. The default query count is
512 and is configurable. Use the same query count, seed, and cache for every
backend compared in one figure.

## Methodology and memory bounds

`flatten` makes a single `float32` NumPy `.npy` memmap (`flat/vectors.npy`),
not an in-memory concatenated database. It reads one compressed window shard at
a time and normalizes rows in `--row-chunk` blocks (default 4,096 rows). The
output normalization is performed even though GINflow's window manifests already
declare normalized vectors, so all backends receive exact unit-length vectors.

Query selection is deterministic and record-balanced:

1. Records are ranked by SHA-256 of their identifier and a public seed.
2. One window per record is selected before any record can receive a second
   query.
3. Its within-record offset is independently SHA-256-derived.

This avoids letting long transcripts dominate the query set while remaining
reproducible from the cache metadata.

Ground truth is exact cosine top-100 over the flattened normalized vectors. The
implementation holds only a configurable database block and a configurable
query block, plus `Q × 100` retained candidates. Ties are descending score then
ascending flat vector ID, including when an equal-score group crosses database
blocks. CPU is always available. With CuPy and a visible CUDA device,
`--engine auto` uses a GPU but still transfers only one database block at a
time; the GEMM stream is synchronized before each host result is merged. GPU
memory stays bounded by
`--database-chunk × dimension × 4` bytes plus cuBLAS workspace. The selected
engine, chunk sizes, elapsed time, and score/ID checksums are written to the
ground-truth manifest.

For the 30k cache, choose a conservative GPU block for an 8 GiB GPU, for
example `--database-chunk 16384 --query-chunk 16`; it is exact but may take
longer. Increase the CPU block only after accounting for the score matrix
(`query_chunk × database_chunk × 4` bytes) and the backend's own RAM use.

## Cache layout and provenance

For `<cache>/<dataset>/`, the important artifacts are:

```text
artifacts/                 # published graphs, embeddings, and windows
nextflow-work/             # persistent, resumable Nextflow work cache
window-cache-request.json  # source/config/workflow checksums
flat/vectors.npy           # normalized, mmap-able database matrix
flat/records.json          # record-to-vector ranges
flat/windows.tsv           # flat vector mapping
flat/flatten-manifest.json
queries/queries.npy
queries/queries.tsv
queries/query-selection.json
ground-truth/ground-truth.npz
ground-truth/ground-truth.json
provenance.json
results/                   # one result JSON per backend/configuration/repeat
```

`provenance.json` includes the source TSV checksum; GINFINITY model/checkpoint
evidence from the window manifests; flattened-vector payload checksum; query
selection checksum; exact-neighbour checksum; Git revision; Python/NumPy
versions; and host CPU/RAM/GPU evidence. The stable IDs named
`embedding_cache_id`, `query_selection_id`, and `ground_truth_cache_id` are
the compatibility keys for result aggregation. They derive only from semantic
inputs and artifact content. Timestamps, elapsed time, absolute cache paths,
and the exact CPU/GPU engine remain provenance fields outside those identities,
so the IDs are stable across repeated generation and cache relocation while
still detecting changed vectors, selected queries, or neighbour IDs.

## Backend result contract

Every timed repeat must be a JSON object satisfying
[`result-schema.json`](result-schema.json). Use the helpers in
`benchmark_utils.py` to calculate recall and create records rather than
reimplementing their fields. Successful records require `k=100`, cosine
metric, a measured build time/index byte count, warmup and timed query counts,
QPS/latency, recall@100, peak RSS/VRAM (when observable), and all three cache
provenance IDs.

Successful rows use `fixed-query-batches-v1`: every timed call has the same
configurable `query_batch_size`; `latency_ms.mean/p50/p95` are actual
milliseconds per query batch; and QPS is total timed queries divided by those
same batch-call durations. The batch size is stored in both `measurement` and
provenance and participates in result comparability, so do not compare points
produced with different batch sizes. Choose positive warm-up and timed query
counts that are exact multiples of the batch size. The default 512 selected
queries works cleanly with a batch size of 32.

Run at least three timed repeats after warm-up. Error and skipped configurations
use `latency_ms: null` and preserve an explanation in `error`; they belong in
the report but are not recall-vs-QPS plot points. All successful points are
retained in the full-range plot; the optional recall-threshold plot is only a
focused operating-range view.

## Fast checks

The fixture tests exercise flattening, deterministic record balancing, exact
bounded-block top-k, result validation, and a dry-run cache invocation without
requiring GINFINITY, Docker, or a GPU:

```bash
python3 -m unittest discover -s benchmarks/tests -v
```

Use `inspect` to retrieve all cache evidence without rereading the source data:

```bash
python3 benchmarks/prepare_cache.py inspect --dataset rouskin_30k
```

## Live progress dashboard

Benchmark runners write one JSON object per configuration repeat under
`<cache>/<dataset>/results/<backend>/`.  To follow partial runs, keep the
self-contained dashboard renderer running while the containers work:

```bash
python3 benchmarks/live_dashboard.py \
  --cache-dir /mnt/ssd_samsung/ginflow-benchmarks \
  --output /mnt/ssd_samsung/ginflow-benchmarks/dashboard/index.html \
  --refresh-seconds 15 --watch
```

Serve that directory from another terminal and open the printed local URL:

```bash
python3 -m http.server 8765 \
  --directory /mnt/ssd_samsung/ginflow-benchmarks/dashboard
```

Open <http://localhost:8765/>.  The page refreshes every 15 seconds and shows
complete and partial repeat groups, errors/skips, all successful recall-versus
QPS points, raw float32 database bytes, serialized index bytes, and separate
build/search RSS and VRAM measurements.  The chart legend identifies each
backend; recall@100 is the horizontal axis and queries/s is the logarithmic
vertical axis.  Hover or focus a point for its full configuration, parameter,
repeat, resource, hardware, and provenance metadata.  Enable **Highlight GPU
indexes** to keep GPU points vivid while fading CPU points.  The HTML is
replaced atomically, so a browser never sees a half-written page.  Enable
**Show Pareto curves** to draw the non-dominated recall/QPS frontier separately
for each backend.  While a runner container is active, the page also shows a
running-job card with backend/dataset, current process stage, elapsed time, CPU,
RSS, and an estimated terminal-configuration progress based on result JSONs.
Raw JSON and backend-specific Markdown reports remain available below the same
cache root.
