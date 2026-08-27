# FAQ and troubleshooting

## Launch errors

### `Specify either --input (optionally with --query) or --query --database, not all three.`

Mode is inferred from flags. Build, or search an existing database, or
build-and-search. You cannot pass `--database` together with `--input`.

### `--query requires --database … or --input`

A search needs something to search against. Pass `--database path/to/index`
or also pass `--input` so this run builds the index first.

### `--database must be a directory …`

Point `--database` at a previous run’s `index/` directory, not the
parent `--outdir` and not `index.faiss` itself. That directory must
contain `windows.tsv`, `meta.json`, residue `embeddings.npz` /
`records.tsv`, and a vector index (`index.faiss`, `cuvs/`,
`cagra.index`, or `index.bin`).

### `--index X does not match the --database (Y)`

The database’s `meta.json` (or on-disk artifact) already selected a
library. Omit `--index` on query-only runs, or pass the matching value.

### `--faiss_gpu requires -profile gpu`

GPU FAISS needs the `faiss-gpu` image and the NVIDIA runtime. Add
`gpu` to `-profile`. `--faiss_gpu` only applies to `--index faiss`
with `--faiss_index flatip` or `flatl2`.

### `--index cagra/ivf GPU build/search requires -profile gpu`

CAGRA and cuVS IVF are GPU-built. To search a CAGRA graph without a
GPU, build with `--cagra_to_hnsw true`, then query with
`--search_device cpu`.

### `--quantize pq cannot be used with --index faiss`

PQ/OPQ need a custom-distance graph: `--index cagra` (recommended) or
`--index hnswlib` (CPU-only build). FAISS and cuVS IVF only accept
`--quantize none` or `sq`.

### `--index hnswlib is the custom-distance PQ/OPQ CPU path`

Uncompressed or SQ windows belong on `--index faiss` or `--index cagra`,
not hnswlib.

### `--candidate_k must be >= --seed_k`

The ANN pool has to be at least as large as the post-rerank seed list.

### Unused index parameter warning

A flag you set does not apply to this `--index` / `--faiss_index`
combination (for example `--faiss_nlist` with `flatip`). It is ignored.
See [Window indexes](indexes.md).

## Runtime

### `ginfinity build-graphs` dies on unmatched brackets

The dot-bracket string must be balanced and the same length as the
sequence. Test tables are pair-closed. For sliced inputs, GINflow
rewrites crossing pairs as unpaired (`.`) on the published subject;
the *source* structure still has to parse.

### CUDA was requested but is unavailable

`-profile gpu` must be set so `EMBED_RNA_GRAPHS` gets
`environment.gpu.yml` and the CUDA Wave image. The published
ginfinity-only Wave tag is CPU PyTorch. Also pass `--allow-nondeterministic-cuda`
(the gpu profile does this).

### FAISS GPU image will not start

GPU FAISS is `pytorch::faiss-gpu=1.10.0` with a **CUDA 12.1** runtime,
for host drivers that report 12.1 or 12.2 (for example 535.x). A
conda-forge CUDA 12.9 FAISS build will not start on those drivers.

### Out of memory on index build

- Default FlatIP stores ~5.6 GB of float32 payload per million windows,
  plus graph/index overhead and original `embeddings.npz`.
- Switch to `--quantize sq` or `--quantize opq --index cagra`.
- See the payload table in [Window indexes](indexes.md).
- `conf/pq_cagra_48gb.config` raises quantizer/CAGRA process memory to
  48 GB (`-c`, not a params file).

### Self-hits dominate the report

Expected when the query is also in the database. Filter
`alignments.tsv` on `query_id != target_id`, or use a query table that
is not a subset of `--input`.

### No seeds / empty `alignments.tsv`

- Lower `--seed_min_similarity` (default 0.8 is strict).
- Raise `--candidate_k` / `--seed_k` on large databases.
- Confirm molecules are ≥ `--window_size` (default 11).
- Confirm `--quantize` / `--index` match how the database was built.

### Plots are missing from `report.html`

`--plot_backend` defaults to `none` and `--plot_sw` defaults to
`false`. Plots are also capped by `--plot_max_pairs` (25 unique pairs
**per query**).

### Nextflow version too old

The pipeline requires Nextflow ≥ 25.10.4 (workflow outputs). Do not
enable `nextflow.preview.output`.

## Concepts

### Is this sequence BLAST?

No. Seeds come from embedding-window cosine, and SW uses embedding
cosine, not a nucleotide matrix. E-values are calibrated on that
score. Sequence identity in `alignments.tsv` (`base_matches`) is a
descriptive statistic, not the alignment objective.

### Windows vs `start`/`end` slices

- `start` / `end` on the TSV change the **graph** GINFINITY builds
  (motif / domain subjects). Coordinates on outputs are relative to
  the core.
- `--window_size` / `--window_stride` cut **already embedded**
  nucleotides into 11-nt FAISS seeds.

Full guide: [Sliced graphs](usage.md#sliced-graphs).

### What is the `index/` directory?

It holds whatever `--index` you chose (FAISS, cuVS, PQ-CAGRA, or hnswlib)
plus residue embeddings and EVD. The generic name reflects that the pipeline
now supports several index libraries.

### Can I add sequences to an existing index?

Not in place. Rebuild with a concatenated structures table, or search
several databases separately and merge `alignments.tsv` yourself.

### Docker vs conda

Docker/Singularity is the usual path. Conda is fully supported and
already lists `nicolas.aira` for GINFINITY, GINFINITY-SW, RNArtistCore,
R4RNA, and PQ-CAGRA. GPU CAGRA/IVF/FAISS still need `-profile gpu` so
the NVIDIA runtime is attached.
