# AGENTS.md

GINflow – Automated Agent & Contributor Guide

Purpose
GINflow is now a BLAST-style Nextflow pipeline for local similarity searches across RNA secondary structures. It consumes a transcript table, generates ginfinity node embeddings, derives sliding-window vectors, indexes the database windows with FAISS, seeds/query windows, clusters hits, runs banded affine-gap Smith–Waterman alignments, and emits a TSV of top-scoring high-scoring pairs (HSPs).

High-Level Workflow
1. **Metadata extraction** (`extract_meta_map`): capture ID, sequence, and structure columns for downstream reporting.
2. **Node embeddings** (`generate_node_embeddings`): call `ginfinity-generate-node-embeddings` (CPU/GPU) to obtain per-node vectors.
3. **Window vectors** (`generate_window_vectors`): concatenate `k` consecutive node embeddings (after batching/merging `split_size` chunks), normalise, and split into database vs query windows based on the ID list.
4. **FAISS index** (`build_faiss_index`): construct the configurable FAISS backend (Flat, IVF, IVFPQ, OPQ+IVFPQ, or HNSW) over database windows.
5. **Seeding** (`query_faiss_index`): search query windows, keep neighbours ≥ `seed_similarity_threshold`, and record diagonals.
6. **Diagonal clustering** (`cluster_seeds`): enforce the two-hit rule by keeping ≥ `cluster_min_seeds` within `cluster_span` nt while allowing small diagonal drift (`cluster_diagonal_tolerance`, `cluster_max_diagonal_span`).
7. **Alignment** (`align_candidates`): gather clustered regions, sample background μ/σ, run banded Smith–Waterman with affine gaps, emit DP traces, and summarise alignment statistics plus sequence/structure snippets.
8. **Outputs**: copy intermediate artefacts (embeddings, windows, seeds, clusters) and final `alignments.tsv`/`alignment_dp.jsonl` into `outdir/` alongside standard Nextflow reports.

Key Modules (modules/)
- `generate_node_embeddings`: ginfinity wrapper producing node-level embeddings.
- `generate_window_vectors`: sliding-window concatenation and metadata emission.
- `build_faiss_index`: new configurable FAISS index builder (Flat/IP/L2, IVF, IVFPQ, OPQ+IVFPQ, HNSW).
- `query_faiss_index`: batch query over FAISS, cosine filtering, and seed TSV generation.
- `cluster_seeds`: diagonal clustering with tolerance-controlled drift and span gating.
- `align_candidates`: background sampling, banded SW, alignment metrics & TSV export.
- Legacy plotting/aggregation/report modules remain in the repository but are no longer wired into `workflows/rna_similarity.nf` (kept for future reporting work).

Containers
containers/ holds Dockerfiles (one per functional group) and pre-built Singularity .sif images (containers/sifs/). Profiles docker or singularity activate the corresponding runtime. The sif_push.sh script can upload built images to a registry; ensure environment vars (if any) are set before use.

Configuration & Parameters (nextflow.config)
Important params (all overridable on the CLI):
- **Input / metadata**: `--input`, `--queries`, `--outdir`, `--header`, `--split_size`, `--id_column`, `--structure_column_name`, `--sequence_column`, `--keep_cols`.
- **Embeddings**: `--node_embeddings_tsv` (reuse cache), `--ginfinity_model_path`, `--num_workers`, `--inference_batch_size`, `--use_gpu` (ginfinity only).
- **Window vectors**: `--window_size`, `--window_stride`.
- **FAISS**: `--index_type`, `--faiss_metric`, `--faiss_nlist`, `--faiss_pq_m`, `--faiss_pq_bits`, `--faiss_opq_m`, `--faiss_hnsw_m`, `--faiss_hnsw_efc`, `--faiss_k`, `--faiss_nprobe`, `--faiss_use_gpu`.
- **Seeding & clustering**: `--seed_similarity_threshold`, `--cluster_span`, `--cluster_min_seeds`, `--cluster_diagonal_tolerance`, `--cluster_max_diagonal_span`.
- **Alignment**: `--alignment_gamma`, `--band_width`, `--band_buffer`, `--band_max_width`, `--alignment_padding`, `--gap_open`, `--gap_extend`, `--xdrop`, `--score_min`, `--score_max`, `--background_samples`, `--random_seed`, `--top_n`.

Profiles
- test: Small synthetic dataset (tests/data/) for fast validation.
- gpu: Enables GPU flags (combine with docker or singularity).
- docker / singularity / conda / mamba: Software stack selection.
- slurm: Submit processes to a SLURM cluster (labels control resources: gpu, medium_memory, high_memory, high_cpu, lightweight).
- local: Forces local executor.

Testing (Quick Start)
Smoke run (CPU-only) that exercises the new BLAST workflow end-to-end:

```
nextflow run main.nf -profile smoke,docker
```

Use `-profile test,<stack>` for a slightly larger dataset. Add `,gpu` to accelerate ginfinity embedding (faiss/align remain CPU-bound).

Artifacts
- `node_embeddings.tsv` – ginfinity outputs (copied to `outdir/`).
- `database_windows.{npy,tsv}` / `query_windows.{npy,tsv}` – sliding-window vectors + metadata.
- `faiss_index/faiss.index` + `faiss_mapping.tsv` + `faiss_index_stats.json` – index artefacts.
- `seeds.tsv`, `clusters.tsv`, `cluster_members.tsv` – seed/cluster diagnostics.
- `alignments.tsv` – final ranked HSP list (includes coordinates, scores, cosine averages, gaps, sequences/structures, and gapped alignment strings).
- `alignment_dp.jsonl` – per-alignment DP trace (JSONL) for downstream audit or visualisation.
- `alignment_pairs.txt` – BLAST-style text dump of gapped alignments.
- Standard Nextflow trace/report/timeline/dag under `outdir/reports/`.

Scoring Details
`align_candidates.py` samples random node pairs (`--background_samples`) to estimate μ₀/σ₀, computes z-scores for cosine similarities, scales by `--alignment_gamma`, clamps to [`--score_min`, `--score_max`], and feeds the resulting matrix into a banded affine-gap Smith–Waterman implementation. Gap costs are `--gap_open` / `--gap_extend`; the minimum band width is `--band_width`, which expands by the observed diagonal span plus `--band_buffer` (capped by `--band_max_width` if set); padding defaults to `--alignment_padding`, and the extension halts once the best score drops by more than `--xdrop`.

Agent / Automation Guidance
When building autonomous agents (e.g. LLM-based assistants) to extend or maintain GINflow:
1. Discovery: Parse nextflow.config to extract params and profiles; map processes in main.nf / workflows/*.nf to module folders.
2. Safety: Validate presence and schema of input TSV (id_column, required columns) before starting compute-heavy stages.
3. Incremental Runs: Detect existing outputs (embeddings.tsv, faiss_index/) and skip regenerate unless upstream params changed.
4. GPU Selection: If params.use_gpu true and gpu_type changed, force rebuild or re-pull GPU-specific containers if they differ.
5. Resource Tuning: For SLURM, adjust labels rather than editing core process definitions (e.g. add withLabel selectors for new memory tiers).
6. Batching: Tune `--split_size` to control rows per embedding job; smaller batches increase concurrency but generate more intermediate files (header rows are mandatory for batching).
7. Caching: Leverage Nextflow resume (-resume) but also surface a summary of cache hit/miss per process in agent logs.
8. Validation Hooks: After `query_faiss_index` ensure each query emits ≤ `faiss_k` seeds and report when none survive the similarity threshold; after `align_candidates` flag NaNs/Infs and ensure `alignments.tsv` is non-empty before completion.
9. Extensibility: Additional post-alignment scoring or visualisation steps should slot in after `align_candidates`, consuming its TSV (or the upstream seeds/clusters files) without mutating existing artefacts.
10. Reporting: Keep HTML generation side-effect free (read-only on inputs, pure write to reports/). Agents adding new plots should append, not overwrite existing ones.
11. Provenance: Embed params.json snapshot (automatic: consider writing a small module exporting params to outdir/params_snapshot.json for full reproducibility).

Suggested Agent Tasks
- add_param_snapshot: Module to dump effective params.
- verify_inputs: Lightweight preflight (column presence, row counts, uniqueness of id_column).
- extend_faiss_metrics: Compute additional statistics (mean intra-cluster distance, etc.).
- enrich_report: Add comparative plots across runs (requires run registry file).

Quality & CI Suggestions
- Add a minimal GitHub Actions workflow invoking the test profile on push.
- Lint: nextflow config validate (nf-core schema if adopted) + python -m ruff (if Python scripts expanded).
- Determinism: Use hist_seed and fix any other RNG seeds in Python helpers.

Contributing
1. Open an issue describing change rationale (performance, feature, bug).
2. Add or update test data only if necessary; keep it small.
3. Ensure test profile passes locally.
4. Document new params in nextflow.config and README plus (if agent-focused) here.

Troubleshooting (Common)
- OOM in embedding step: Reduce inference_batch_size or apply high_memory label.
- Empty query results: Confirm queries file IDs exist in embeddings TSV.
- GPU not used: Confirm -profile gpu plus docker.runOptions includes --gpus all (or singularity --nv) and host has drivers.

License
See LICENSE.

Contact
Open issues / pull requests for improvements.

## CI & Automated Testing (Agents)

Purpose: Provide rapid feedback on branch changes, especially for agent-created branches (codex/*).

Key Workflows:
- GINflow CI Test (ci.yml): Runs on push / PR to main, master, and codex/** branches. Supports manual dispatch.

Branch-based Modes:
- codex/* branches default to `smoke` mode (very fast, minimal resources) with profiles: smoke,docker,(gpu if available)
- Other branches default to `test` mode (moderate coverage) with profiles: test,docker,(gpu if available)
- `full` mode may be selected manually via workflow_dispatch input when deeper validation needed.

Dispatch Inputs:
- mode: auto | smoke | test | full (auto selects smoke for codex/* else test)
- use_gpu: auto | yes | no (auto adds gpu profile if runtime present)
- extra_profiles: comma-separated list appended (e.g. `singularity` or `monitor`)
- resume: true|false to enable `-resume` cache reuse

Smoke vs Test vs Full:
- smoke: Minimal k (100), top_n (3), reduced workers/batch size, no plotting of distributions, skips SVG drawings, fastest (<2 min)
- test: Medium size (default small test dataset), some plots, moderate runtime
- full: Uses only explicitly requested profiles (start from docker) and full default params

Artifacts:
- ci_<mode> artifact containing: test_results/ or smoke_results/, .nextflow.log, ci_summary/summary.json
- summary.json fields: branch, sha, mode, profiles, gpu_runtime, requested_gpu, resume, success

Agent Workflow Loop (Recommended):
1. Create or update branch `codex/<task-name>` with changes.
2. Push commits; CI auto-runs in smoke mode.
3. Download summary.json & logs via GitHub API. Parse success.
4. If failure: parse .nextflow.log for the first process failure (search for 'ERROR ~' or process name) and patch code.
5. Re-push changes; repeat until success.
6. Optionally dispatch `mode=test` for deeper coverage before opening PR.
7. For performance-sensitive modifications (embedding logic), manually trigger `mode=full`.

Failure Parsing Heuristics:
- Container/image pull errors: retry once; if persistent, switch to conda profile (`extra_profiles=conda`).
- Memory errors (Exit status 137 / OOM): reduce inference_batch_size or add label adjustments (future enhancement: dynamic memory scaling).
- GPU missing: rerun with `use_gpu=no` if gpu profile added but runtime absent.

Security & Safety:
- Avoid adding secrets directly into workflow edits; use repository secrets for new credentials.
- Keep changes to nextflow.config additive; do not delete existing profiles (agents should append new profile blocks).

Extending CI (Agent Tasks):
- Add new quick validation modules: create `profile smoke` variations or stub-run mode (future optional).
- Append new artifacts (e.g., param snapshot) by writing JSON to `ci_summary/` and updating upload step.

Status Badge (add to README):
`![CI](https://github.com/nicoaira/GINflow/actions/workflows/ci.yml/badge.svg)`

Triggering Manually via API (for agents):
`POST /repos/nicoaira/GINflow/actions/workflows/ci.yml/dispatches` JSON body: `{ "ref": "codex/my-branch", "inputs": { "mode": "smoke" } }` with a token having `workflow` scope.

Cache Considerations:
- Resume flag reuses previous work/ hashes; only use when underlying input & params stable.
- Agents should omit resume after structural pipeline edits to force full rebuild.

End of CI Agent Section.
