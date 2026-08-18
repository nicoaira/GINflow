# GINflow vector-index benchmark results

Use this table only after `benchmarks/validate_results.py` accepts the raw per-repeat result files. The generated `benchmark_results.md` is the authoritative populated version; this is a compact reporting template, not a source of measurements.

| Dataset / hardware scope | Backend | Configuration | Run ID | Raw database | Recall@100, median (min–max) | QPS, median (min–max) | Timed repeats | Index size | Build time | Build peak RSS | Build peak VRAM | Search peak RSS | Search peak VRAM | Availability / caveat |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `<dataset_id> / <hardware_id>` | `<faiss|scann|ngt|cuvs>` | `<parameter_label>` | `<run_id>` | `<bytes>` | `<median (min–max)>` | `<median (min–max)>` | `<n>` | `<bytes>` | `<seconds>` | `<bytes>` | `<bytes or CPU>` | `<bytes>` | `<bytes or CPU>` | `<measured|skipped|error and reason>` |

All successful configurations belong in the table and full-range recall-vs-QPS plot. A companion focused plot may show only configurations with median recall@100 strictly greater than 0.80. Keep below-threshold and unavailable configurations in the table and explain their status. Do not compare rows unless the validator confirms the same dataset window count, dimension, metric, k, warm-up/timed query counts, fixed query-batch protocol, query-set hash, exact-ground-truth hash, cache provenance, runner revision, and hardware scope.
