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
├── samples.csv
└── pipeline_info/
```

`faiss/` is published on build runs. `seeds.tsv` is published on query runs.

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

`--database` for a later search run should point at this directory.

## seeds.tsv

Query-mode hits above `--seed_min_similarity`. Columns: `query_id`, `query_start`, `query_end`, `target_id`, `target_start`, `target_end`, `score`, `rank`. Self-hits are kept.
