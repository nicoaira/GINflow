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
├── samples.csv
└── pipeline_info/
```

## samples.csv

Index of every published shard. Columns: `id`, `graphs`, `metadata`, `embeddings`, `manifest`. Paths point at the files under `graphs/` and `embeddings/`.

## graphs/

One GINFINITY graph shard per input chunk (`--shard_size` records), in a subdirectory named after the shard id.

- `*.safetensors` — node/edge tensors
- `*.json` — identifiers, sequences, structures, graph spec, and a content checksum

## embeddings/

Per-nucleotide embeddings from `ginfinity embed-graphs`, one subdirectory per shard.

- `*.npz` — one array per `transcript_id`
- `*.manifest.json` — package/model versions, device, and record shapes
