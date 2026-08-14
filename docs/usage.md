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

## Run

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
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
