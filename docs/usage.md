# nicoaira/ginflow: Usage

## Input

`--input` and `--query` are structures tables. Default format is TSV with these **required** columns:

| Column | Meaning |
|---|---|
| `transcript_id` | Unique molecule id |
| `sequence` | RNA sequence (`A/C/G/U`; `T` is accepted and converted) |
| `secondary_structure` | Dot-bracket structure, same length as the sequence |

A table with only those three columns is encoded as full molecules. Extra columns (for example `rfam_family`) are allowed and ignored.

Optional `start` / `end` columns select **sliced graphs** — see [Sliced graphs](#sliced-graphs).

Override column names with `--id_column`, `--sequence_column`, `--structure_column`, and `--delimiter`.

Large tables are split into shards of `--shard_size` **input rows** (default `50`) before graph construction. A row that lists several windows expands to several graphs inside that shard. Query search uses the same size unless you set `--search_shard_size` (input rows per `SEARCH_FAISS` task).

## Sliced graphs

Use this when you want an embedding (and later a search/alignment subject) for a motif, domain, or other interval of a longer RNA, without throwing away pairing that sits just outside that interval.

This is **not** the same as `--window_size` / `--window_stride`. Those cut the *already computed* per-nucleotide embeddings into 11-nt FAISS seeds. `start` / `end` change the graph GINFINITY builds **before** embedding.

### How to write the table

Add two optional columns named `start` and `end` (or rename them with `--start_column` / `--end_column`). Both must be present if either is. Coordinates are **0-based half-open**, the same convention as a Python slice: the core is `sequence[start:end]`. Position `0` is the first nucleotide. The window must be non-empty and inside the sequence.

| What you want | `start` | `end` | What GINflow builds |
|---|---|---|---|
| Full molecule | empty | empty | One subject, original `transcript_id` |
| One interval | `50` | `113` | One subject, id `{transcript_id}:50-113` |
| Several intervals | `34,40` | `90,95` | Two subjects: `{id}:34-90` and `{id}:40-95` |

You can mix those cases in the same file. Omit the columns entirely if every row is a full molecule.

```tsv
transcript_id	sequence	secondary_structure	start	end
full-mol	GGGAAACCCUUUUGGG	......(((....)))		
one-slice	GGGAAACCCUUUUGGG	......(((....)))	9	16
two-slices	…long RNA…	…dot-bracket…	34,40	90,95
```

- `full-mol` — no window; treated as before.
- `one-slice` — core is `sequence[9:16]` (`UUUUGGG`, 7 nt). Published as `one-slice:9-16`.
- `two-slices` — two overlapping windows `[34, 90)` and `[40, 95)`. Each is its own database record **and** its own query. Whitespace around commas is ignored. The two lists must have the same length.

A worked mixed file is in the repo: [`tests/data/sliced_structures.tsv`](../tests/data/sliced_structures.tsv) (two full molecules, one annotated Rfam core, and one row with overlapping `34,40` / `90,95` windows).

### What happens at runtime

1. GINFINITY builds a graph for each window. By default it also keeps **crossing-pair context**: if a core nucleotide is paired with a base outside `[start, end)`, that partner (and a small neighbourhood) stays in the graph so message passing can use it.
2. After embedding, **only the core nucleotides** are kept. The embedding has shape `(end - start, 128)`.
3. From then on the pipeline treats that window as an independent molecule: FAISS records, seeds, clusters, alignments, and `report.html` all use the slice id and coordinates **relative to the core** (0 is the first nucleotide of the window, not of the source RNA).
4. The published sequence/structure for a slice is the core substring. Pairs that crossed the cut are written as unpaired (`.`) so the subject stays a legal, balanced structure.

### Flags

| Flag | Default | Meaning |
|---|---|---|
| `--start_column` | `start` | Name of the optional window-start column |
| `--end_column` | `end` | Name of the optional window-end column |
| `--no_slices` | `false` | Ignore `start`/`end` even if the file has those columns; encode full molecules |
| `--keep_paired_neighbours` | `true` | Keep crossing-pair partners outside the core in the graph |
| `--context_hops` | `4` | How far to expand around those partners |

`--context_hops 1` keeps only the partner. Higher values also walk backbone, pairing, and skip-2 edges around that partner. Raising the hop count makes the graph larger but does **not** change the embedding length: you still get one vector per core nucleotide.

Turn context off with `--keep_paired_neighbours false` (the graph is then just the core window).

`--input` and `--query` both accept sliced tables. You can search a sliced query against a full-molecule database, or the reverse: each slice is just another subject with its own id.

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
    --query queries.tsv \
    --keep_paired_neighbours true \
    --context_hops 4 \
    --outdir results
```

## Run modes

The mode is inferred from the flags. Do not pass `--input`, `--query`, and `--database` together.

| Flags | Mode |
|---|---|
| `--input` | Build a window database from the structures table |
| `--query` + `--database` | Embed the query table and search an existing `faiss/` directory |
| `--input` + `--query` | Build the database, search it, and publish it for later runs |

`--database` must be a previous run's `faiss/` directory (`windows.tsv`, `meta.json`, `embeddings.npz`, `records.tsv`, and a vector index such as `index.faiss`, `scann/`, `ngt/`, or `cuvs/`). cuVS databases require `-profile gpu` for the later search as well. If `--input` and `--query` are the same file, embeddings are computed once.

Window search parameters:

| Flag | Default | Meaning |
|---|---|---|
| `--window_size` | `11` | Nucleotides per window |
| `--window_stride` | `1` | Step between window starts |
| `--seed_k` | `50` | Neighbours kept per query window before the threshold |
| `--seed_min_similarity` | `0.8` | Minimum cosine similarity |
| `--search_shard_size` | `--shard_size` | Query records per search task |
| `--index` | `faiss` | Index **library**: `faiss`, `scann`, `ngt`, or `cuvs` |
| `--faiss_index` | `FlatIP` | FAISS **structure** (only with `--index faiss`) |
| `--ngt_index` | `NGT` | NGT **structure** (`NGT`, `QG`, or `QBG`) when `--index ngt` |
| `--cuvs_index` | `CAGRA` | cuVS **structure** (`CAGRA`, `IVF`, or `IVF-PQ`) when `--index cuvs`; requires `-profile gpu` |

Each window is the concatenation of `w` per-nucleotide 128-d vectors (1408-d), L2-normalized so an inner product is cosine similarity. These sliding windows are built **after** embedding; they are not the optional `start` / `end` columns on the structures table (see [Sliced graphs](#sliced-graphs)).

**Index library vs FAISS type, every library-specific flag, GPU, and unused-parameter warnings:** [Window indexes](indexes.md). Passing `--scann_*` with `--index faiss`, or `--faiss_nlist` with `FlatIP`, prints a launch warning and ignores the flag. `--faiss_index ScaNN` is an error; use `--index scann`.

Seeds are then clustered along nearby diagonals and each cluster is aligned with [GINFINITY-SW](https://github.com/nicoaira/GINFINITY-SW). Alignment runs on a padded crop of the cluster (`--align_pad`, default 32), not the full molecules.

| Flag | Default | Meaning |
|---|---|---|
| `--cluster_span` | `80` | Max gap from a seed to the current cluster box |
| `--cluster_min_seeds` | `2` | Drop singleton clusters |
| `--cluster_diagonal_tolerance` | `12` | Allowed diagonal drift when adding a seed |
| `--cluster_max_diagonal_span` | `96` | Max diagonal breadth of one cluster |
| `--cluster_max_seed_rank` | `10` | Ignore seed hits worse than this rank |
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

A ready-made pair lives in the repo: `-profile test` builds the 1200-sequence database, and `tests/data/example_queries.tsv` is four queries under 200 nt from RF00001, RF00003, RF01725, and RF01852. See the README for the exact commands. `tests/data/sliced_structures.tsv` is a small mixed table (full molecules + single slice + overlapping slices) you can pass as `--input` and/or `--query`.

CPU is the default embed device. For GPU:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --input structures.tsv \
    --outdir results
```

`-profile gpu` sets `accelerator = 1` on `EMBED_RNA_GRAPHS`. That switches the process to `environment.gpu.yml` and the GPU Wave image (`ginfinity` + `pytorch-gpu` + `cuda-version`), and passes `--device cuda --allow-nondeterministic-cuda`. Graph construction stays on the CPU image.

FAISS GPU is **opt-in** and separate from embedding:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --index faiss \
    --faiss_gpu \
    --faiss_index IVFFlat \
    --input structures.tsv \
    --outdir results
```

`--faiss_gpu` without `-profile gpu` is a pipeline error. With both, `BUILD_FAISS_INDEX` and `SEARCH_FAISS` switch to `environment.gpu.yml` (`pytorch::faiss-gpu=1.10.0`, CUDA 12.1 runtime) and the FAISS GPU Wave image. That runtime matches host drivers that report CUDA 12.1 or 12.2 (for example NVIDIA driver 535.x). A newer conda-forge CUDA 12.9 build will not start on those drivers.

ScaNN is CPU-only and uses its own image (`environment.scann.yml`, `scann==1.4.2`):

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --index scann \
    --input structures.tsv \
    --query queries.tsv \
    --outdir results
```

cuVS is GPU-only and uses the RAPIDS/CuPy image for both index construction and search. Choose `CAGRA`, `IVF` (IVF-Flat), or `IVF-PQ` with `--cuvs_index`; all cuVS runs require `-profile gpu`:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --index cuvs --cuvs_index CAGRA \
    --input structures.tsv --query queries.tsv \
    --outdir results
```

## Outputs

See [output.md](output.md).
