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
│   ├── index.faiss   # FAISS types
│   ├── cuvs/         # stock CAGRA or IVF-Flat
│   ├── hnsw/         # CPU HNSW converted from CAGRA when --cagra_to_hnsw
│   ├── cagra.index   # PQ-CAGRA graph
│   ├── index.bin     # hnswlib PQ graph
│   ├── quantization/ # node SQ/PQ/OPQ artifacts
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
│   ├── r4rna/
│   └── sw/
├── samples.csv
└── pipeline_info/
```

`faiss/` is published on build runs and now also contains residue `embeddings.npz` plus `records.tsv` so a later `--query --database` run can align. `seeds.tsv`, `clusters.tsv`, and `alignments.tsv` are published on query runs.

When the input table has `start` / `end` windows, each window is stored under `{transcript_id}:{start}-{end}` and treated as its own molecule from embeddings onward. See [Sliced graphs](#sliced-graphs).

## samples.csv

Index of every published shard. Columns: `id`, `graphs`, `metadata`, `embeddings`, `manifest`, `windows`, `windows_manifest`.

## graphs/

One GINFINITY graph shard per input chunk (`--shard_size` records), in a subdirectory named after the shard id.

- `*.safetensors` — node/edge tensors (sliced shards also store `residue_index` and `node_roles`)
- `*.json` — identifiers, full source sequences/structures, graph spec, and a content checksum. A sliced row is stored under `{transcript_id}:{start}-{end}`.

## embeddings/

Per-nucleotide embeddings from `ginfinity embed-graphs`, one subdirectory per shard.

- `*.npz` — one array per identifier, shape `(L, 128)`, stored as float16 by default. `--ginfinity_full_precision` changes GINFINITY's model inference to float32; it does not change the default embedding storage dtype. For a slice, `L` is the core window (`end - start`), not the source molecule.
- `*.manifest.json` — package/model versions, device, and record shapes (`core_length`, and `start`/`end` when sliced)

## windows/

Sliding windows of the node embeddings (`--window_size`, default 11, stride 1).

- `*.windows.npz` — one `(n_windows, 1408)` array per `transcript_id` (concatenated, L2-normalized)
- `*.windows.manifest.json` — window size, shapes, and sequences shorter than `w`

## faiss/

Reusable window index of every database window. Default is exact inner product (`IndexFlatIP` / `--faiss_index flatip`). `meta.json` records `backend`, `index_type`, `--quantize`, and graph/IVF parameters. See [indexes.md](indexes.md).

- `index.faiss` — FAISS index (always stored as a CPU index, even after a GPU build)
- `cuvs/` — serialized GPU cuVS CAGRA or IVF-Flat
- `hnsw/` — CPU HNSW converted from uncompressed/SQ CAGRA (`--cagra_to_hnsw`)
- `cagra.index` — custom PQ-CAGRA graph (`--index cagra` with `--quantize pq|opq`)
- `index.bin` — hnswlib custom-distance PQ graph
- `quantization/` — node codebook, SDC table, OPQ rotation, SQ scales; not a substitute for `embeddings.npz`
- `windows.tsv` — `faiss_id`, `transcript_id`, `start`, `end`
- `meta.json` — window geometry, model fingerprint, counts, and index settings
- `embeddings.npz` — per-nucleotide embeddings, one array per identifier
- `records.tsv` — `transcript_id`, `sequence`, `secondary_structure`. For a slice the identifier is `{id}:{start}-{end}` and the sequence/structure are the core window (pairs that cross the cut are written as unpaired), so later search/alignment treat that window as its own molecule.
- `evd.json` — Karlin–Altschul λ, K, database residue count, and the reverse-sequence null fit

When `--index hnswlib` is selected, the workflow also publishes `quantization/`
and `windows_quantized/` at the output root. These contain the CPU centroid-code
candidate representation and are not used by SW. GPU HNSWLIB additionally
stores its int8-scaled original-window CAGRA data inside the serialized index;
the original float16 embeddings above remain authoritative for exact scoring,
SW, alignment, and plots.

`--database` for a later search run should point at this directory.

## seeds.tsv

Query-mode hits above `--seed_min_similarity`. Columns: `query_id`, `query_start`, `query_end`, `target_id`, `target_start`, `target_end`, `score`, `rank`. Self-hits are kept. Without exact reranking, HNSWLIB writes the raw position-wise centroid-score sum and applies the equivalent threshold multiplied by `--window_size`; with `--exact_rerank true`, it writes the normal full-window cosine score from the original float16 embeddings. The same run publishes `rerank_metrics.json` with batch/device/worker settings and elapsed time.

## clusters.tsv

Diagonal HSP clusters. A cluster starts at the highest-scoring unused seed and grows by the best remaining seed whose query/target box is within `--cluster_span` and whose diagonal stays inside `--cluster_diagonal_tolerance` / `--cluster_max_diagonal_span`. Singletons are dropped.

## alignments.tsv / alignments.txt

GINFINITY-SW local HSPs are collapsed by query-target pair, ranked by ascending database E-value. Coordinates are 0-based half-open on the independent subject (the full molecule, or the core window of a slice). Each row reports `total_score` (the sum of disjoint HSP scores), `max_score` (the strongest HSP), `alignment_count`, `evalue` (K m N e^{−λS_total}), and `evalue_pair` (same formula with the target length instead of the full database). The legacy `score` column aliases `total_score`. HSP provenance is retained in `cluster_ids`, `hsp_scores`, and `hsp_spans`; the row also includes the subject `query_sequence` / `query_structure` / `target_sequence` / `target_structure`, and gapped `query_aligned` / `target_aligned` strings for the first HSP. `alignments.txt` contains one pair header followed by the six-line rendering for each HSP.

## report.html

Self-contained search report written on every query run. Open it in a browser: hits are grouped by query, filterable by E-value, and paginated (10 / 25 / 50 / 100 / 150 per page, default 10). The default theme is light; `--report_theme dark` writes a gray-900 page, and a toggle in the masthead switches themes when viewing. The summary table and each hit card show database E, pair E, bits, total score, maximum HSP score, and HSP count, with the aligned span on both molecules. A floating **Back to top** control appears after scrolling and returns to the report header. RNArtistCore and SW plots use a two-column Query | Target (or cosine | SW scores) panel. R4RNA is a single full-width alignment arc plot. The masthead uses `docs/images/ginflow_logo.svg`; the favicon is `docs/images/ginflow_icon.svg`.

## plots/

Published when `--plot_backend` is `rnartistcore`, `r4rna`, or `both`, and/or when `--plot_sw` is set. Each query gets its own draw task. RNArtistCore writes a query SVG and a target SVG per HSP row. R4RNA writes one alignment-coordinate SVG per HSP row (query arcs above, target arcs below, identity ribbon between the backbones). SW plots write a cosine SVG and a score SVG (traceback on the scores). `--plot_max_pairs` (default 25) counts unique query-target pairs per query, not HSP rows or SVG files; every HSP in each selected pair is rendered, so the SVG count can be larger. The aligned span, shared pairs, matching bases, and SW path use `--plot_highlight_colour`.

- `plots/rnartistcore/*.svg` — RNArtistCore 2Ds
- `plots/r4rna/*.svg` — R4RNA alignment arc diagrams
- `plots/sw/*.svg` — crop cosine and substitution-score heatmaps

## Sliced graphs

If `structures.tsv` included `start` / `end` columns, every published artifact uses the **expanded** subjects, not the original TSV rows.

| Artifact | What you see for a window `[start, end)` |
|---|---|
| Identifier | `{transcript_id}:{start}-{end}` (full molecules keep the original id) |
| `graphs/*.json` | Full source sequence/structure; the graph node set is core + optional pairing context |
| `embeddings/*.npz` | One `(end - start, 128)` array — core nucleotides only |
| `windows/*` | Sliding 11-nt seeds over that core embedding |
| `faiss/records.tsv` | Core sequence and a pair-closed core structure (crossing pairs become `.`) |
| `faiss/windows.tsv`, `seeds.tsv`, `alignments.tsv` | Coordinates are 0-based on the core: `0` is `sequence[start]`, not the 5′ end of the source RNA |

Two overlapping windows on the same source row (`34,40` / `90,95`) appear as two ids, two embeddings, and two query groups in `report.html`. How to write the table: [usage.md — Sliced graphs](usage.md#sliced-graphs).

## See also

- [Clustering and alignment](clustering-and-alignment.md)
- [E-values](statistics.md)
- [Plotting and report](plotting.md)
- [Parameters](parameters.md)
