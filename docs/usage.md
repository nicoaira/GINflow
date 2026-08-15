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

Large tables are split into shards of `--shard_size` records (default `50`) before graph construction. Query search uses the same size unless you set `--search_shard_size` (records per `SEARCH_FAISS` task).

## Run modes

The mode is inferred from the flags. Do not pass `--input`, `--query`, and `--database` together.

| Flags | Mode |
|---|---|
| `--input` | Build a FAISS database from the structures table |
| `--query` + `--database` | Embed the query table and search an existing `faiss/` directory |
| `--input` + `--query` | Build the database, search it, and publish it for later runs |

`--database` must be a previous run's `faiss/` directory (`index.faiss`, `windows.tsv`, `meta.json`, `embeddings.npz`, `records.tsv`). If `--input` and `--query` are the same file, embeddings are computed once.

Window search parameters:

| Flag | Default | Meaning |
|---|---|---|
| `--window_size` | `11` | Nucleotides per window |
| `--window_stride` | `1` | Step between window starts |
| `--seed_k` | `50` | Neighbours kept per query window before the threshold |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity |
| `--search_shard_size` | `--shard_size` | Query records per FAISS search task |

Each window is the concatenation of `w` per-nucleotide 128-d vectors (1408-d), L2-normalized so the FAISS inner product is cosine similarity.

Seeds are then clustered along nearby diagonals and each cluster is aligned with [GINFINITY-SW](https://github.com/nicoaira/GINFINITY-SW). Alignment runs on a padded crop of the cluster (`--align_pad`, default 32), not the full molecules.

| Flag | Default | Meaning |
|---|---|---|
| `--cluster_span` | `80` | Max gap from a seed to the current cluster box |
| `--cluster_min_seeds` | `2` | Drop singleton clusters |
| `--cluster_diagonal_tolerance` | `12` | Allowed diagonal drift when adding a seed |
| `--cluster_max_diagonal_span` | `96` | Max diagonal breadth of one cluster |
| `--cluster_max_seed_rank` | `10` | Ignore FAISS hits worse than this rank |
| `--align_pad` | `32` | Extra nucleotides on each side of the cluster crop |
| `--evd_samples` | `1000` | Reverse-sequence null alignments used to fit λ and K |
| `--evd_max_length` | `400` | Max null-sequence length during EVD calibration |

Every search writes `report.html` — a standalone page of ranked hits, aligned spans, and (if requested) structure and SW-matrix plots. The default theme is light (`--report_theme light`); pass `--report_theme dark` for a gray-900 report. The page also has a theme toggle. Hits are paginated (10 per page by default; 25, 50, 100, or 150). RNArtistCore and SW plots sit in a two-column Query | Target (or cosine | SW scores) grid. R4RNA is one full-width alignment arc plot per pair.

Optional structure plots (`--plot_backend rnartistcore`, `r4rna`, or `both`) draw the query and target 2Ds. RNArtistCore colours the aligned span (`--plot_highlight_colour`); the rest of the molecule is gray. R4RNA draws both structures on alignment coordinates (query arcs up, target arcs down) with a per-column base-identity ribbon. Default is `none`. `--plot_sw` adds the crop cosine matrix and the substitution-score matrix with the Smith–Waterman traceback. Each query runs its own draw task. `--plot_max_pairs` (default 25) is per query and counts alignment pairs. Plots are also inlined in the report. Each draw process uses `task.cpus` workers (6 with the default `process_medium` label).

| Flag | Default | Meaning |
|---|---|---|
| `--plot_backend` | `none` | `rnartistcore`, `r4rna`, `both`, or `none` |
| `--plot_sw` | `false` | Draw crop cosine + SW score matrices (traceback on the score plot) |
| `--plot_max_pairs` | `25` | Max alignment pairs plotted per query |
| `--plot_highlight_colour` | `#24B064` | Colour of the aligned span, shared pairs, matching bases, and SW traceback (nf-core green) |
| `--report_theme` | `light` | `light` or `dark` colour theme for `report.html` |

Alignments are ranked by ascending database E-value. E = K m N exp(−λS), where m is the query length and N is the number of residues in the searchable database. λ and K are fit at database-build time from Smith–Waterman scores of reversed real embeddings (preserves local embedding correlation, destroys homology). The legacy `K = exp(−λμ)` conversion is not used; K comes from a length-aware Gumbel MLE so that μ = ln(Kmn)/λ.

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

A ready-made pair lives in the repo: `-profile test` builds the 1200-sequence database, and `tests/data/example_queries.tsv` is four queries under 200 nt from RF00001, RF00003, RF01725, and RF01852. See the README for the exact commands.

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
