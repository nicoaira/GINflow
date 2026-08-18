# Vector-index benchmark TODO

**Status:** Task 1 is complete and committed (`234f0c3`). Task 2 has its research guide, cache/measurement harness, reporting tools, and backend runners, but the real index runs and final report are still outstanding.

This file is the hand-off checklist for finishing Task 2. It records only work that is still required; generated caches and benchmark results stay outside git.

## Verified state

- [x] Dedicated build/search processes exist for FAISS, ScaNN, NGT, and cuVS; workflow dispatch, GPU validation, backend detection, and focused tests were completed in Task 1.
- [x] Hardware and tuning guidance is documented in [`docs/indexes.md`](indexes.md), with version scope and conservative memory/VRAM caveats.
- [x] The reproducible cache, exact-ground-truth, fixed-query-batch timing, result schema, validator, plotter, and Markdown reporting helpers are implemented under [`benchmarks/`](../benchmarks/).
- [x] The 6k cache is complete: 117 graph/embedding/window shards, 839,188 normalized vectors (4,726,306,816 bytes / 4.402 GiB raw float32 payload), 512 deterministic queries, and exact GPU top-100 ground truth. Its external cache request is `window-cache-request-2bf2dbe6ab4e0a05`; the ground-truth cache is `ground-truth-a42ca20a4ae4f91e`.
- [x] The 30k cache is complete end-to-end: request `window-cache-request-6fe176b0e111d942`, 4,115,576 normalized vectors (23,178,924,032 bytes / 21.587 GiB raw float32 payload), 512 deterministic queries (`query-selection-8b51ce6be33b7a8b`), and GPU exact top-100 ground truth (`ground-truth-0301b276bea26b64`).
- [x] Synthetic and pinned-container runner tests exist for all four backends. Real 6k results are being retained externally; no final cross-backend conclusion has been claimed yet.
- [x] Real 6k frontiers now exist for FAISS (69 successful rows, 23 configurations), ScaNN (54 successful rows, 18 configurations), and cuVS (75 successful rows, 25 configurations plus explicit error/skip rows); their external reports include raw payload, serialized index, build/search memory, and full/focused plots.
- [x] A live, self-contained dashboard is available at `/mnt/ssd_samsung/ginflow-benchmarks/dashboard/index.html`; `benchmarks/live_dashboard.py --watch` rewrites it atomically every 15 seconds and includes partial groups, errors/skips, all successful points, database/index sizes, and separate build/search memory columns.
- [ ] Real 6k NGT/QG/QBG is still running; the QBG build has already demonstrated roughly 23 GiB RSS on this host before search measurement.
- [x] 30k dry-run capacity plans are recorded: ScaNN and all NGT structures exceed the host safe RAM estimate, cuVS exceeds the 8 GiB GPU limit because the raw upload is 23,178,924,032 bytes, and only compressed FAISS IVFPQ settings pass the initial CPU gate.

## Remaining work

### 1. Close the runner audit before spending benchmark time

- [x] Fix and test the ScaNN and cuVS `MemAvailable` preflight parsing.
- [x] Enforce the documented minimum of three repeats for `rouskin_6k` and `rouskin_30k` in FAISS, ScaNN, and cuVS, including smoke/plan overrides.
- [x] Make the NGT production path benchmark the serialized-and-reloaded index (`serialize_index` then `load_index`), not the still-live build object.
- [x] Rerun the host runner suite: 46 tests passed, 7 expected pinned-environment tests skipped. Pinned FAISS, ScaNN, NGT, and cuVS smoke/integration checks had already passed in their backend images; the independent final audit agent was interrupted by workspace credit exhaustion, so preserve that limitation in the delivery notes.
- [ ] Tighten the cuVS GPU capacity preflight: on the 8 GiB card, two 6k CAGRA settings passed the estimate but produced CUDA `std::bad_alloc`; retain those error rows and calibrate a larger safety margin/workspace model before calling the gate predictive.

### 2. Finish the 30k benchmark cache

- [x] Flatten the completed 30k window shards into the normalized vector memmap and record count/window count.

  ```bash
  rtk python3 benchmarks/prepare_cache.py flatten \
    --cache-dir /mnt/ssd_samsung/ginflow-benchmarks --dataset rouskin_30k
  ```

- [x] Select the deterministic 512-record-balanced query set.

  ```bash
  rtk python3 benchmarks/prepare_cache.py select-queries \
    --cache-dir /mnt/ssd_samsung/ginflow-benchmarks \
    --dataset rouskin_30k --query-count 512
  ```

- [x] Compute exact cosine top-100 ground truth with bounded chunks using the GPU engine (35.5 s, database chunk 32,768, query chunk 16).

  ```bash
  rtk python3 benchmarks/prepare_cache.py ground-truth \
    --cache-dir /mnt/ssd_samsung/ginflow-benchmarks \
    --dataset rouskin_30k --k 100 --engine auto
  ```

- [x] Inspect both cache manifests and verify matching input/cache/query/ground-truth hashes before running any index.

### 3. Run the real index frontiers

For each dataset/backend/parameter configuration, retain at least three timed repeats, fixed query-batch size/count, warmup count, recall@100, QPS, latency percentiles, index-build metadata, and the exact cache/ground-truth provenance. Resource accounting is mandatory: record the raw database size (`window_count × dimension × 4` bytes for the normalized float32 matrix), serialized index size on disk, **build peak** host RSS/VRAM, **search peak** host RSS/VRAM, and the sampler/measurement scope. Build and search are separate resource questions: construction may need training/graph/quantization workspace and temporary copies, while search needs the loaded index, query batches, and search scratch. Never substitute one for the other. If a resource is not observable, write `null` with a reason; if a capacity gate refuses a configuration, retain the estimated requirement, available resource, safe limit, and explicit skip reason. Keep skipped and capacity/error rows in the result files rather than silently dropping them.

The final per-configuration table must therefore include at least: dataset ID, indexed-window count, vector dimension, raw database bytes (and GiB), backend/index type, every build/search parameter, serialized index bytes (and GiB), build time, build peak RSS/VRAM, search peak RSS/VRAM, repeat count, recall/QPS/latency statistics, and resource provenance. Do not replace measured peak memory with the raw matrix size: report both separately. The raw payload is exact for this float32 cache; measured peaks are precise for the pinned container, host, driver, thread/batch settings, and sampler, not universal guarantees for every deployment.

- [ ] Run [`benchmarks/run_faiss.py`](../benchmarks/run_faiss.py) on 6k, then the feasible 30k CPU/GPU configurations. Cover the FlatIP, IVFFlat/nprobe, IVFPQ, and HNSW/efSearch ladders; mark GPU configurations that exceed available VRAM as skipped.
- [ ] Run [`benchmarks/run_scann.py`](../benchmarks/run_scann.py) on 6k and 30k using production tree-AH/SOAR frontiers. Do not report brute-force or unsupported modes as production points.
- [ ] Run [`benchmarks/run_ngt.py`](../benchmarks/run_ngt.py) on 6k and 30k for production NGT/QG/QBG. Keep native epsilon/edge experiments visibly separate and explicitly experimental; apply the QBG workspace/RAM gate.
- [ ] Run [`benchmarks/run_cuvs.py`](../benchmarks/run_cuvs.py) on GPU for the feasible 6k CAGRA/IVF/IVF-PQ frontiers. For 30k, record capacity skips where the raw 1408-dimensional upload or index workspace cannot fit the 8 GiB GPU; do not substitute an unrecorded smaller dataset.
- [ ] Confirm every real result row validates against [`benchmarks/result-schema.json`](../benchmarks/result-schema.json), including query and ground-truth hashes and the strict metric/`k=100` contract.
- [ ] Audit every real result file for the resource fields above; verify that each configuration has a database-size value and either measured host/GPU memory or an explicit unobservable/capacity reason.
- [ ] Keep build and search memory columns distinct in the final report (`measurement.build_peak_*` versus top-level search `peak_*`); document that cuVS/NVIDIA VRAM is device-wide `nvidia-smi` usage and therefore an upper bound rather than process-attributed memory.
- [ ] Add/execute a cold serialized-load search measurement for each configuration. The current runners measure timed search after an in-process build/serialize/reload, so their search RSS can include allocator pages retained from construction; treat it as a warm-process upper bound until an isolated `SEARCH_*`-equivalent process records the cold search footprint.

### 4. Generate the benchmark report

- [ ] Aggregate the validated JSON/JSONL/CSV outputs with [`benchmarks/validate_results.py`](../benchmarks/validate_results.py) and the reporting helpers.
- [ ] Generate deterministic full-range recall-vs-QPS SVGs containing all successful points, plus focused SVGs containing only points with **recall > 0.80** (strict inequality). Render/inspect each SVG once for readable axes, legends, thresholds, and one-point scopes.
- [ ] Complete a user-facing report (for example `docs/index-benchmarks.md`) with:
  - dataset/cache IDs, query-selection IDs, ground-truth hashes, and software/container revisions;
  - host CPU/RAM/GPU/VRAM and the memory projections from `docs/indexes.md`;
  - the exact parameter grids and repeat/timing protocol;
  - recall/QPS/latency tables and plots for 6k and 30k;
  - a database-size and memory table for every backend/parameter configuration (raw vectors, serialized index, peak RSS, peak VRAM, and estimates/skips);
  - capacity skips, failed configurations, and the 8 GiB limitations;
  - tuning recommendations tied to hardware, database size, recall target, and latency target;
  - reproducible commands and a statement that build/RSS measurements are harness measurements, not universal production capacity guarantees.
- [ ] Ensure the report does not imply a backend is comparable when its cache, metric, query set, ground truth, hardware, or runner revision differs.

### 5. Final verification and delivery

- [ ] Re-run the focused Task 1 module/workflow tests after the module-resource/symlink cache fix, plus the complete `benchmarks/tests` suite.
- [ ] Run `rtk git diff --check`, schema/config validation, and `rtk repowise update`; review the resulting diff for accidental generated data or external-cache paths.
- [ ] Keep `/mnt/ssd_samsung/ginflow-benchmarks` and raw result directories untracked; commit only source, tests, documentation, report metadata, and intentionally selected plots.
- [ ] Commit Task 2 separately after the real benchmark report is complete, with a message describing the research, harness, measurements, and plots.

## Acceptance gate

Task 2 is complete only when both datasets have reproducible cache identities and exact top-100 ground truth, all four backend runners have validated real results or explicit capacity skips, the report retains all successful measurements and also provides the focused strict `recall > 0.80` view, the tuning/hardware guidance is tied to those measurements, and the final report plus tests are committed.
