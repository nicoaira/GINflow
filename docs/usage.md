# nicoaira/ginflow: Usage

## Input

`--input` is one structures table. Default format is TSV with these required columns:

| Column | Meaning |
|---|---|
| `transcript_id` | Unique molecule id |
| `sequence` | RNA sequence (`A/C/G/U`; `T` is accepted and converted) |
| `secondary_structure` | Dot-bracket structure, same length as the sequence |

Extra columns (for example `rfam_family`) are allowed and ignored.

Override names with `--id-column`, `--sequence-column`, `--structure-column`, and `--delimiter`.

Large tables are split into shards of `--shard_size` records (default `50`) before graph construction.

## Run modes

The mode is inferred from the flags. Do not pass `--input`, `--query`, and `--database` together.

| Flags | Mode |
|---|---|
| `--input` | Build a FAISS database from the structures table |
| `--query` + `--database` | Embed the query table and search an existing `faiss/` directory |
| `--input` + `--query` | Build the database, search it, and publish it for later runs |

`--database` must be a previous run's `faiss/` directory (`index.faiss`, `windows.tsv`, `meta.json`). If `--input` and `--query` are the same file, embeddings are computed once.

Window search parameters:

| Flag | Default | Meaning |
|---|---|---|
| `--window_size` | `11` | Nucleotides per window |
| `--window_stride` | `1` | Step between window starts |
| `--seed_k` | `50` | Neighbours kept per query window before the threshold |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity |

Each window is the concatenation of `w` per-nucleotide 128-d vectors (1408-d), L2-normalized so the FAISS inner product is cosine similarity.

## Run

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
    --outdir results
```

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --query queries.tsv \
    --database results/faiss \
    --outdir search_results
```

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
    --query queries.tsv \
    --outdir results
```

CPU is the default embed device. For GPU:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --input structures.tsv \
    --outdir results
```

`-profile gpu` sets `accelerator = 1` on `EMBED_RNA_GRAPHS`. That switches the process to `environment.gpu.yml` and the GPU Wave image (`ginfinity` + `pytorch-gpu` + `cuda-version`), and passes `--device cuda --allow-nondeterministic-cuda`. Graph construction stays on the CPU image.

## Outputs

See [output.md](output.md).
