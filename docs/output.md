# nicoaira/ginflow: Output

Published under `--outdir` (or `-output-dir`) by the entry-workflow `output` block:

```
outdir/
├── graphs/
│   └── <shard>/
│       ├── <shard>.safetensors
│       └── <shard>.json
├── embeddings/
│   └── <shard>/
│       ├── <shard>.npz
│       └── <shard>.manifest.json
├── windows/
│   └── <shard>/
│       ├── <shard>.windows.npz
│       └── <shard>.windows.manifest.json
├── faiss/
│   ├── index.faiss
│   ├── windows.tsv
│   └── meta.json
├── seeds.tsv
├── clusters.tsv
├── cluster_members.tsv
├── alignments.tsv
├── alignments.txt
├── report.html
├── plots/
│   ├── rnartistcore/
│   └── r4rna/
├── samples.csv
└── pipeline_info/
```

`faiss/` is published on build runs and now also contains residue `embeddings.npz` plus `records.tsv` so a later `--query --database` run can align. `seeds.tsv`, `clusters.tsv`, and `alignments.tsv` are published on query runs.

## samples.csv

Index of every published shard. Columns: `id`, `graphs`, `metadata`, `embeddings`, `manifest`, `windows`, `windows_manifest`.

## graphs/

One GINFINITY graph shard per input chunk (`--shard_size` records), in a subdirectory named after the shard id.

- `*.safetensors` — node/edge tensors
- `*.json` — identifiers, sequences, structures, graph spec, and a content checksum

## embeddings/

Per-nucleotide embeddings from `ginfinity embed-graphs`, one subdirectory per shard.

- `*.npz` — one array per `transcript_id`, shape `(L, 128)`
- `*.manifest.json` — package/model versions, device, and record shapes

## windows/

Sliding windows of the node embeddings (`--window_size`, default 11, stride 1).

- `*.windows.npz` — one `(n_windows, 1408)` array per `transcript_id` (concatenated, L2-normalized)
- `*.windows.manifest.json` — window size, shapes, and sequences shorter than `w`

## faiss/

Reusable exact inner-product index (`IndexFlatIP`) of every database window.

- `index.faiss` — FAISS index
- `windows.tsv` — `faiss_id`, `transcript_id`, `start`, `end`
- `meta.json` — window geometry, model fingerprint, and counts
- `embeddings.npz` — per-nucleotide embeddings, one array per `transcript_id`
- `records.tsv` — `transcript_id`, `sequence`, `secondary_structure`
- `evd.json` — Karlin–Altschul λ, K, database residue count, and the reverse-sequence null fit

`--database` for a later search run should point at this directory.

## seeds.tsv

Query-mode hits above `--seed_min_similarity`. Columns: `query_id`, `query_start`, `query_end`, `target_id`, `target_start`, `target_end`, `score`, `rank`. Self-hits are kept.

## clusters.tsv

Diagonal HSP clusters. A cluster starts at the highest-scoring unused seed and grows by the best remaining seed whose query/target box is within `--cluster_span` and whose diagonal stays inside `--cluster_diagonal_tolerance` / `--cluster_max_diagonal_span`. Singletons are dropped.

## alignments.tsv / alignments.txt

GINFINITY-SW local alignments of each cluster crop, ranked by ascending database E-value. Coordinates are 0-based half-open on the original molecules. Extra columns: `bit_score`, `evalue` (K m N e^{−λS}), `evalue_pair` (same formula with the target length instead of the full database), plus the full `query_sequence` / `query_structure` / `target_sequence` / `target_structure`. `alignments.txt` is the six-line RNA rendering.

## report.html

Self-contained search report written on every query run. Open it in a browser: hits are grouped by query, filterable by E-value, and each hit shows the aligned span on both molecules. Structure plots are inlined when `--plot_backend` was set.

## plots/

Published when `--plot_backend` is `rnartistcore`, `r4rna`, or `both`. Each alignment can produce a query SVG and a target SVG, up to `--plot_max`. The aligned span is coloured with `--plot_highlight_colour`; unaligned residues stay gray.

- `plots/rnartistcore/*.svg` — RNArtistCore 2Ds
- `plots/r4rna/*.svg` — R4RNA arc diagrams
