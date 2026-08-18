# Benchmark result validation and plotting

The backend runners emit one JSON object for every timed repeat (`schema_version: "ginflow-benchmark-v1"`). The canonical field definition is [`result-schema.json`](result-schema.json); these reporting tools consume individual `.json` files, JSONL/NDJSON, CSV, or TSV recursively. Non-result manifests such as the JSON schema are ignored deliberately.

Validate before plotting:

```bash
python3 benchmarks/validate_results.py \
  /mnt/ssd_samsung/ginflow-benchmarks/rouskin_sample_6k/results \
  --output /mnt/ssd_samsung/ginflow-benchmarks/rouskin_sample_6k/validation.json
```

Render the static SVGs and populated Markdown table only after validation succeeds:

```bash
python3 benchmarks/plot_results.py \
  /mnt/ssd_samsung/ginflow-benchmarks/rouskin_sample_6k/results \
  --output-dir /mnt/ssd_samsung/ginflow-benchmarks/rouskin_sample_6k/report
```

The output directory contains:

- `validation.json` — machine-readable validation status, errors, warnings, and per-configuration aggregates.
- `benchmark_results.md` — a compact result table with skipped/error configurations retained.
- `recall_qps_all_<dataset>_<hardware>_<scope>.svg` — full-range vector plot containing every successful aggregate in a comparable scope.
- `recall_qps_gt_<threshold>_<dataset>_<hardware>_<scope>.svg` — optional focused vector plot containing only points strictly above the requested recall threshold. The stable scope suffix prevents overwriting distinct dataset names that sanitize to the same filename.

## What is accepted

For a `status: "ok"` row, the validator requires the fixed benchmark definition (`metric: "cosine"`, `k: 100`), a positive dataset window count and dimension, a non-empty backend/configuration/run/repeat identity, positive warm-up and timed-query counts, positive finite QPS, and recall@100 in `[0, 1]`. It also requires the query and exact-ground-truth SHA-256 digests, index/build measurements, ordered per-batch latency percentiles, and provenance containing `git_commit`, `runner_version`, `hardware_id`, `embedding_cache_id`, `query_selection_id`, and `ground_truth_cache_id`. Peak RSS and VRAM are retained when observable and otherwise remain `null` rather than being fabricated.

Successful rows must declare the fixed measurement protocol: `measurement.protocol: "fixed-query-batches-v1"`, a positive query-batch size and timed-batch count, matching timed-query totals, per-batch latency units, and the QPS scope. This prevents a figure from silently mixing throughput measured with different query batching.

At least two valid timed repeats per configuration are required by default. Use `--min-repeats 1` only for a smoke check, never for a reported benchmark. Duplicate repeat indices and incompatible database/query/ground-truth, cache, runner, or query-batching identities inside a dataset/hardware scope are errors. `skipped` and `error` rows are retained for coverage reporting but never plotted; they must identify the backend, dataset, configuration, run, and reason.

## How the plots are kept honest

Each point is the median recall@100 and median QPS across its timed repeats. Its horizontal whisker is the measured QPS min–max; the table exposes both recall and QPS ranges. The table also reports the raw normalized float32 database payload (`dataset_window_count × dimension × 4`), serialized index bytes, **build peak** RSS/VRAM from `measurement.build_peak_*`, and **search peak** RSS/VRAM from the top-level `peak_*` fields for every configuration. Build and search peaks are not interchangeable. The current search peak is sampled after an in-process serialized reload and may include allocator pages retained from construction; it is a warm-process upper bound until a cold `SEARCH_*`-equivalent process is measured. The raw payload is exact for the pinned cache; measured peaks are scoped to the host/container/driver/thread/batch settings, and cuVS VRAM is device-wide `nvidia-smi` usage rather than process attribution. The plotter draws no connecting lines or fitted frontier, so it cannot imply performance at an unmeasured setting. It uses a logarithmic QPS axis and emits both a full-range plot and, when eligible points exist, a focused plot whose **median** recall@100 is strictly greater than `0.80`. `--recall-threshold` changes the focused eligibility rule and its y-axis/grid annotation. Results with a different hardware ID always receive a separate plot.

Run the synthetic checks without any index library installed:

```bash
python3 benchmarks/tests/test_results_tools.py -v
```

The fixture is intentionally synthetic and is only a contract test. It must never be copied into a benchmark report as a real measurement.
