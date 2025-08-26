# AGENTS.md

GINflow – Automated Agent & Contributor Guide

Purpose
GINflow is a Nextflow pipeline for computing similarity between RNA transcript secondary structures. It ingests a table of transcript metadata / structures, produces fixed-length embeddings, builds a FAISS index, performs nearest-neighbour queries, aggregates pairwise scores, and generates plots plus HTML reports.

High-Level Workflow
1. Input parsing: Load tabular transcript data (TSV/CSV) with an ID column (default transcript_id) and a structure column.
2. Window generation (generate_windows): Optionally splits long structures into fixed-size windows (param split_size) and masks low-confidence positions (mask_threshold).
3. Embedding inference (generate_embeddings): Runs a model to produce numeric embeddings (GPU optional via -profile gpu). Batch size controlled by inference_batch_size and num_workers.
4. FAISS index build (build_faiss_index): Creates a vector index plus mapping files.
5. Query phase (query_faiss_index): Looks up nearest neighbours (faiss_k) for provided queries (queries file) or all items.
6. Distance sorting & filtering (sort_distances, filter_top_contigs): Orders neighbour lists and trims to top_n.
7. Scoring (aggregate_score, merge_query_results): Computes per-pair scores and aggregated metrics using alpha*, beta*, gamma, percentile parameters.
8. Visualization (plot_distances, plot_score, draw_contig_svgs, draw_unagg_svgs): Generates distributions, per-contig SVGs (configurable with draw_contig_svgs / run_unaggregated_report), and summary figures.
9. Reporting (generate_aggregated_report, generate_unaggregated_report): Produces HTML reports (pairs tables, metrics, plots) in outdir/reports.

Key Modules (modules/)
- generate_windows: Sequence/structure segmentation & masking.
- generate_embeddings: Vector representation generation (CPU/GPU).
- build_faiss_index: FAISS index + mapping TSV.
- query_faiss_index: K-NN retrieval.
- sort_distances / merge_query_results: Post-processing of raw distances.
- filter_top_contigs: Selects top_n candidates per query.
- aggregate_score: Applies scoring formula and parameterized weighting.
- plot_distances / plot_score: Histogram and score distribution plotting.
- draw_contig_svgs / draw_unagg_svgs: Structure SVG rendering.
- generate_aggregated_report / generate_unaggregated_report: HTML report writers.

Containers
containers/ holds Dockerfiles (one per functional group) and pre-built Singularity .sif images (containers/sifs/). Profiles docker or singularity activate the corresponding runtime. The sif_push.sh script can upload built images to a registry; ensure environment vars (if any) are set before use.

Configuration & Parameters (nextflow.config)
Adjust params { ... } for defaults; override on the CLI or with profiles:
- Input & Output: --input <file>, --outdir <dir>, --header, --id_column, --keep_cols.
- Structure handling: --subgraphs, --L, --keep_paired_neighbors, --structure_column_name, --mask_threshold.
- Embeddings: --num_workers, --inference_batch_size, --use_gpu (or -profile gpu), --gpu_type (t4|a100).
- FAISS / Query: --faiss_k, --top_n, --queries, --embeddings_tsv, --faiss_index.
- Aggregation: --alpha1 --alpha2 --beta1 --beta2 --gamma --percentile.
- Plot/Report toggles: --plot_distances_distribution, --plot_score_distribution, --draw_contig_svgs, --run_aggregated_report, --run_unaggregated_report.

Profiles
- test: Small synthetic dataset (tests/data/) for fast validation.
- gpu: Enables GPU flags (combine with docker or singularity).
- docker / singularity / conda / mamba: Software stack selection.
- slurm: Submit processes to a SLURM cluster (labels control resources: gpu, medium_memory, high_memory, high_cpu, lightweight).
- local: Forces local executor.

Testing (Quick Start)
Run the bundled test workflow (produces test_results/):

nextflow run main.nf -profile test,docker,gpu

(If no GPU available, omit ,gpu.)

Artifacts
- Embeddings: <outdir>/embeddings.tsv
- FAISS index: <outdir>/faiss_index/*.index plus mapping TSV.
- Distances / Scores: pairs_scores*.tsv, distances.sorted.tsv
- Reports & Plots: <outdir>/reports/, plots/, drawings/

Scoring Parameters
The score combines structure-aware weighting (alpha1, alpha2) and distance shaping (beta1, beta2, gamma) and is further summarized via percentile selection. Tune these cautiously; re-run only downstream scoring & reporting if embeddings and raw distances unchanged (consider adding process-level publishDir reuse to avoid recompute).

Agent / Automation Guidance
When building autonomous agents (e.g. LLM-based assistants) to extend or maintain GINflow:
1. Discovery: Parse nextflow.config to extract params and profiles; map processes in main.nf / workflows/*.nf to module folders.
2. Safety: Validate presence and schema of input TSV (id_column, required columns) before starting compute-heavy stages.
3. Incremental Runs: Detect existing outputs (embeddings.tsv, faiss_index/) and skip regenerate unless upstream params changed.
4. GPU Selection: If params.use_gpu true and gpu_type changed, force rebuild or re-pull GPU-specific containers if they differ.
5. Resource Tuning: For SLURM, adjust labels rather than editing core process definitions (e.g. add withLabel selectors for new memory tiers).
6. Caching: Leverage Nextflow resume (-resume) but also surface a summary of cache hit/miss per process in agent logs.
7. Validation Hooks: After query_faiss_index ensure each query has <= faiss_k neighbors; after aggregate_score ensure no NaN; halt with explicit diagnostics.
8. Extensibility: New similarity metrics should slot in after query_faiss_index and before aggregate_score via an inserted module consuming distances TSV.
9. Reporting: Keep HTML generation side-effect free (read-only on inputs, pure write to reports/). Agents adding new plots should append, not overwrite existing ones.
10. Provenance: Embed params.json snapshot (automatic: consider writing a small module exporting params to outdir/params_snapshot.json for full reproducibility).

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
