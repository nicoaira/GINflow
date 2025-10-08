# Alignment Phase (align_candidates)

This document explains how GINflow’s alignment phase works and describes every parameter that influences it. The alignment phase consumes clustered seed hits and produces a ranked list of local high‑scoring pairs (HSPs) using a banded Smith–Waterman algorithm with affine gaps and X‑drop termination. Scores are derived from background‑adjusted cosine similarities between node embeddings.

- Inputs: clustered seeds, node embeddings, and transcript metadata
- Core steps: background estimation → region selection → banded affine‑gap SW with X‑drop → result summarisation
- Outputs: `alignments.tsv`, `alignment_dp.jsonl` (optional), `alignment_pairs.txt` (optional), `alignment_plots/` (optional)

## Inputs

The Nextflow module `modules/align_candidates/main.nf` invokes `bin/align_candidates.py` with:

- `--cluster-members` (TSV): per‑seed rows with at least `cluster_id`, `query_transcript`, `target_transcript`, `query_window_start`, `query_window_end`, `target_window_start`, `target_window_end`, `diagonal`, and seed similarity columns.
- `--cluster-summaries` (TSV): one row per cluster including `cluster_id`, `seed_count`, and `max_similarity`.
- `--embeddings` (TSV): node‑level embeddings with columns `transcript_id`, `node_index`, and `embedding_vector` (comma‑separated floats). Vectors are assumed to be L2‑normalised; dot products equal cosine similarity.
- `--meta` (TSV): transcript metadata containing sequence/structure columns; used to slice and emit context strings.

Column names for IDs, node indices, embeddings, sequence, and structure are configurable (see Parameters).

**Seeding (query_faiss_index)**
- Purpose: search query window vectors against a FAISS index of database windows and emit BLAST-style ungapped seeds.
- Inputs: `faiss.index`, `database_windows.npy/tsv`, `query_windows.npy/tsv` (from `generate_window_vectors`).
- Similarity: converts FAISS distances to cosine per metric: `ip` uses inner product directly; `l2` maps via `cosine = 1 - (d2/2)` assuming unit-normalized vectors.
- Filtering: keep neighbours with `similarity >= --similarity-threshold` and drop self-hits where `query_transcript == target_transcript`.
- Diagonal: `diagonal = target_window_start - query_window_start` using 0-based, half-open window starts; 0 means aligned starts, positive means target start is to the right of query, negative the opposite.
- Output: `seeds.tsv` with columns `query_transcript`, `target_transcript`, `query_window_start`, `query_window_end`, `target_window_start`, `target_window_end`, `query_window_index`, `target_window_index`, `similarity`, `rank`, `diagonal`. Also writes `query_stats.json` with counts (`query_windows`, `seeds_kept`, `top_k`, `similarity_threshold`).
- Parameters (Nextflow → CLI): `--faiss_k` → `--top-k`, `--seed_similarity_threshold` → `--similarity-threshold`, `--faiss_metric` → `--metric` (`ip|l2`), `--faiss_nprobe` → `--nprobe` (IVF only), `--faiss_use_gpu` → `--use-gpu`.

Notes
- Window metadata fields required: `transcript_id` (or configured `--id-column`), `window_index`, `window_start`, `window_end`.
- Ensure window vectors are approximately unit-length when using `ip` to interpret scores as cosine.
- Each query window returns up to `--top-k` neighbours before filtering; the final `seeds.tsv` may be empty if all hits fall below threshold.

**Seed Clustering (cluster_seeds)**
- Purpose: enforce a two-hit rule and group seeds into diagonal-consistent clusters for downstream banded alignment.
- Grouping: process seeds per `(query_transcript, target_transcript)` pair ordered by `query_window_start`.
- Acceptance criteria for adding a seed to the current active cluster:
  - Proximity: `query_window_start` and `target_window_start` gaps to the previous seed are both `<= --cluster-span` nucleotides.
  - Diagonal tolerance: the seed’s `diagonal` lies within `[active_diagonal_min - --diagonal-tolerance, active_diagonal_max + --diagonal-tolerance]`.
  - Maximum diagonal span: if `--max-diagonal-span > 0`, the span `candidate_max - candidate_min` after adding the seed must remain `<= --max-diagonal-span`.
- Cluster finalization: when a seed fails the criteria, close the active cluster and emit it only if it contains `>= --min-seeds` members; then start a new cluster with the current seed. After iterating, emit the last cluster if it passes the size threshold.
- Cluster summary fields (`clusters.tsv`): `cluster_id`, `query_transcript`, `target_transcript`, `diagonal` (median of member diagonals), `seed_count`, `query_start`, `query_end`, `target_start`, `target_end`, `max_similarity`, `diagonal_min`, `diagonal_max`, `diagonal_span`.
- Membership mapping (`cluster_members.tsv`): one row per seed with `cluster_id` plus all seed fields from `seeds.tsv`.
- Statistics: optional `cluster_stats.json` with total `clusters` and `seeds` retained.
- Parameters (Nextflow → CLI):
  - `--cluster_span` → `--cluster-span` (default 80): max gap between neighbouring seeds in both query and target coordinates.
  - `--cluster_min_seeds` → `--min-seeds` (default 2): minimum hits required to keep a cluster.
  - `--cluster_diagonal_tolerance` → `--diagonal-tolerance` (default 12): allowed diagonal drift around the running min/max when extending a cluster.
  - `--cluster_max_diagonal_span` → `--max-diagonal-span` (default 96; 0 disables): hard cap on total diagonal range within a cluster.

Tuning Tips (Seeding & Clustering)
- Few or no seeds: lower `--seed_similarity_threshold`, increase `--faiss_k`, or choose a more exhaustive index/metric (e.g., `flat_ip`, higher `--faiss_nprobe`). Verify that queries exist in the embeddings and that vectors are normalized.
- Over-fragmented clusters: increase `--cluster-span` and/or `--diagonal-tolerance` modestly; consider reducing `--max-diagonal-span` only if diagonal drift is implausible biologically.
- Over-merged clusters: decrease `--cluster-span` or `--diagonal-tolerance`; enforce a tighter `--max-diagonal-span`.
- Seed density: larger `--window_size` or `--window_stride` in `generate_window_vectors` reduce/increase seed density and affect clustering granularity.

## Algorithm

1) Background μ/σ estimation
- Purpose: derive a null model for cosine similarities so that per‑position scores are comparable across pairs.
- Method: uniformly sample `--background-samples` random node pairs (with `--random-seed`) across all available transcripts and compute dot products. Estimate μ₀ and σ₀ from these samples; if σ₀ < 1e‑6, set σ₀ = 1.0 to avoid division by zero.

2) Candidate region selection
- For each cluster, compute the median diagonal and its span from member seeds (`diagonal_min`, `diagonal_max`, `diagonal_span`).
- Expand seed bounds with `--padding` to obtain `[query_start, query_end)` and `[target_start, target_end)` slices from the node embeddings.
- Derive `diag_offset = median_diagonal + query_start − target_start` to align the band’s centre with the cluster’s implied offset between query and target.

3) Adaptive banded Smith–Waterman with affine gaps and X‑drop
- Band width: start from `--band-width` and adapt to `max(band_width, diagonal_span + band_buffer)`. The adaptive value is rounded up to an even number; if `--band-max-width > 0`, cap at that value.
- Scoring: at DP cell (i, j), compute dot = `query[i] ⋅ target[j]`, convert to a background‑adjusted score:
  score = clamp( `--gamma` × (dot − μ₀)/σ₀, `--score-min`, `--score-max` )
- Affine gaps: three‑state SW (match `M`, gap‑in‑query `X`, gap‑in‑target `Y`) with costs `--gap-open` and `--gap-extend`.
- X‑drop termination: after finishing a row, if `best_score − row_best > --xdrop`, stop extension early.
- Traceback: follow pointers from the best cell/state to recover the local alignment path.

4) Result summarisation and outputs
- Coordinates: convert local slice coordinates back to global `[query_start, query_end)`, `[target_start, target_end)`.
- Metrics: length, `avg_cosine` across matched cells, `gap_opens`, `gap_bases`, `seed_count`, `max_seed_similarity`, `diagonal_*`, and `band_width`.
- Context strings: if both structure strings are present, emit `aligned_query_structure`, `aligned_target_structure`, and an optional `alignment_mask` ("|" at matches); else fall back to sequence strings if available.
- Outputs:
  - `alignments.tsv`: top results globally (sorted by `score`, limited by `--top-n`).
  - `alignment_dp.jsonl` (optional): DP trace for each reported alignment (compact per‑step records).
  - `alignment_pairs.txt` (optional): BLAST‑style gapped alignment text for quick inspection.
  - `alignment_plots/` (optional): PNG heatmaps of banded scoring matrices when `--plot-scoring-matrices` is enabled.

## Parameters and Effects

Names map from Nextflow params to `align_candidates.py` CLI flags.

- `--alignment_gamma` → `--gamma` (float, default 1.5)
  Scales z‑scores before clamping. Higher values emphasise similarity differences; too large may saturate at `--score-max`.

- `--score_min` / `--score_max` → `--score-min` / `--score-max` (floats, defaults −4, 8)
  Clamp per‑cell scores to stabilise the DP against extreme z‑scores or outliers. Tightening reduces dynamic range; widening can increase sensitivity but risk runaway scores.

- `--band_width` → `--band-width` (int, default 96)
  Minimum DP band width in nucleotides. Wider bands tolerate larger indel/diagonal drift at the cost of compute. Lower bound of the adaptive band.

- `--band_buffer` → `--band-buffer` (int, default 32)
  Additional slack added to the observed diagonal span when adapting the band. Increase to accommodate noisy clusters; decrease for speed when clusters are tight.

- `--band_max_width` → `--band-max-width` (int, default 0 = unlimited)
  Upper bound on the band after adaptation. Use to cap worst‑case compute for long regions.

- `--xdrop` (float, default 50)
  X‑drop termination threshold. Larger values allow longer extensions (potentially slower); smaller values terminate earlier (faster, potentially missing weak tails).

- `--gap_open` / `--gap_extend` → `--gap-open` / `--gap-extend` (floats, defaults 12, 2)
  Affine gap penalties. Increase to discourage gaps; decrease to allow more gapped alignments. The ratio controls typical gap lengths.

- `--alignment_padding` → `--padding` (int, default 32)
  Nucleotides added to both sides of the seed‑derived windows before alignment. Increase if alignments truncate near edges; decrease to speed up.

- `--background_samples` (int, default 10000)
  Number of random node pairs to estimate background μ₀/σ₀. More samples stabilise estimates at a cost of a short upfront compute step.

- `--random_seed` (int, default 42)
  RNG seed for reproducible background sampling.

- `--top_n` → `--top-n` (int, default 50)
  Maximum number of alignments reported globally (after sorting by `score`).

- `--plot_scoring_matrices` (boolean, default false)
  When enabled, renders PNG heatmaps of the banded Smith–Waterman scoring matrices for the reported alignments into `alignment_plots/`, showing the full query/target windows with colour restricted to the evaluated band.

- Column/field selectors (strings)
  - `--id_column` (default `transcript_id`): transcript identifier column in embeddings/meta tables.
  - `--node_index_column` (default `node_index`): node order within a transcript in the embeddings table.
  - `--embedding_column` (default `embedding_vector`): comma‑separated float vector per node.
  - `--sequence_column` / `--structure_column`: column names in the meta table used to slice and emit context strings. Auto‑detected if omitted (`*sequence*` / `*struct*`).

## Output Columns (alignments.tsv)

Key fields written by `align_candidates.py` (when available):
- `query_id`, `target_id`, `cluster_id`, `score`, `alignment_length`
- `query_start`, `query_end`, `target_start`, `target_end`
- `avg_cosine`, `gap_opens`, `gap_bases`, `seed_count`, `max_seed_similarity`
- `diagonal_median`, `diagonal_min`, `diagonal_max`, `diagonal_span`, `band_width`
- Context strings (conditional):
  - `query_sequence`, `target_sequence`, `query_sequence_region`, `target_sequence_region`
  - `query_structure`, `target_structure`, `query_structure_region`, `target_structure_region`
  - `aligned_query_sequence`, `aligned_target_sequence`
  - `aligned_query_structure`, `aligned_target_structure`
  - `alignment_mask` ("|" at aligned positions when strings present)
  - `alignment_content` ("structure" or "sequence")

If no valid alignments survive scoring (`path` empty or `score <= 0`), an empty TSV is emitted; optional trace/text files are created but may be empty.

## Tuning Tips

- Truncated alignments: increase `--alignment_padding`, `--band_width`, and/or `--band_buffer`; relax `--band_max_width` if capping.
- Excessive runtime: reduce `--band_max_width`, `--band_buffer`, and `--alignment_padding`; increase `--xdrop` only cautiously (it increases runtime if set too high).
- Over‑gapped results: raise `--gap_open`/`--gap_extend`.
- Weak sensitivity: increase `--alignment_gamma` moderately; widen `--score_min`/`--score_max`; consider slightly larger band settings.
- Noisy μ/σ: increase `--background_samples`.

## Reproducibility

- Results depend on seed clustering and the random background sampling; fix `--random_seed` and Nextflow’s `-resume` for deterministic reruns on identical inputs/params.

## Notes and Edge Cases

- Adaptive band parity: the diagonal‑span‑derived target is rounded up to an even number to keep band symmetry around the central diagonal.
- Numerical stability: if σ₀ is extremely small, it is set to 1.0 to avoid exploding z‑scores.
- Content selection: structure alignments are preferred when both structure columns are present; otherwise, sequence strings are used when available.
