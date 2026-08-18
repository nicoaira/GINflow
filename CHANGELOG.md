# nicoaira/ginflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]


### `Added`

- ScaNN seed search (`--index scann`): [Google ScaNN](https://github.com/google-research/google-research/tree/master/scann) (`scann==1.4.2`) as a CPU-only alternative to FAISS. Auto-selects brute-force (<20k windows), AH+reorder (<100k), or tree+AH+reorder. ScaNN knobs are `--scann_leaves`, `--scann_leaves_to_search`, `--scann_reorder`, `--scann_ah_dim`, `--scann_anisotropic`, `--scann_soar`. Artifacts go to `faiss/scann/` instead of `index.faiss`. `--faiss_gpu` with ScaNN is an error.
- Launch warning when a library-specific flag does not apply to `--index` / `--faiss_index` (for example `--scann_reorder` with FAISS, or `--faiss_nlist` with `flatip`). Guide: [docs/indexes.md](docs/indexes.md).
- Additional FAISS index types (`--faiss_index`): `flatip` (default), `flatl2`, `hnsw`, `ivfflat`, `lsh`, `sq`, `pq`, `ivfsq`, `ivfpq`, `ivfpqr`, matching [Faiss indexes](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes). IVF/PQ/HNSW/LSH/SQ knobs are `--faiss_nlist`, `--faiss_nprobe`, `--faiss_pq_m`, `--faiss_pq_nbits`, `--faiss_pq_m_refine`, `--faiss_hnsw_m`, `--faiss_hnsw_ef_construction`, `--faiss_hnsw_ef_search`, `--faiss_lsh_nbits`, `--faiss_sq_type`.
- Optional FAISS GPU (`--faiss_gpu`) for `flatip`, `flatl2`, `ivfflat`, `ivfpq`, and `ivfsq`. Requires `-profile gpu` so `BUILD_FAISS_INDEX` / `SEARCH_FAISS` get `accelerator = 1` and the `faiss-gpu` image; otherwise the pipeline errors. GPU-incompatible types also error.
- Optional `start`/`end` columns on the structures table build GINFINITY sliced graphs (one independent subject/query per window, including several comma-separated windows on the same row). Defaults: `--keep_paired_neighbours` with `--context_hops 4`. Mixed examples are in `tests/data/sliced_structures.tsv`. Per-query alignment and plot filenames use `baseName` so two slices of the same accession do not collide.
- GINflow logo in `docs/images/ginflow_logo.svg` (README header and search report masthead) and `docs/images/ginflow_icon.svg` (report favicon).
- Metro-map pipeline schematic (`docs/images/ginflow_metro.svg`) generated with [nf-metro](https://github.com/seqeralabs/nf-metro) from `docs/images/ginflow_metro.mmd`.
- FAISS seed search: `GENERATE_WINDOWS` slices concatenated `w=11` windows, `BUILD_FAISS_INDEX` writes a reusable `IndexFlatIP` database, and `SEARCH_FAISS` returns seeds above `--seed_min_similarity`. Modes are inferred from `--input` / `--query` / `--database`.
- Seed clustering (`CLUSTER_SEEDS`) and GINFINITY-SW alignment (`ALIGN_CLUSTERS`) on a padded crop of each cluster. The FAISS directory now also packs residue embeddings and sequences so query-only runs can align.
- Database E-values: `ESTIMATE_EVD` fits Karlin–Altschul λ, K from reverse-sequence GINFINITY-SW scores; `ALIGN_CLUSTERS` ranks hits by ascending E = K m N exp(−λS).
- Optional alignment plots: `DRAW_RNARTISTCORE` and `DRAW_R4RNA` (`--plot_backend`). Conda packages `nicolas.aira::rnartistcore=0.4.6` (OpenJDK 17–21) and `nicolas.aira::r-r4rna=2.0.9`, with Wave-frozen Docker and Singularity images. Each process draws molecules in parallel with `task.cpus` (6 on `process_medium`).
- Optional SW-matrix plots: `DRAW_SW` (`--plot_sw`) writes crop cosine and substitution-score SVGs with the Smith–Waterman traceback on the score plot. Uses the existing GINFINITY-SW container. Inlined in `report.html`.
- Search HTML report: `WRITE_REPORT` writes a self-contained `report.html` on every query run (ranked hits, span rails, inlined plots).
- Diverse Rfam test tables: 10-sequence `-profile smoke_test` (`tests/data/smoke_test_structures.tsv`) and 1200-sequence `-profile test` (`tests/data/test_structures.tsv`), rebuilt from 12 Rfam families using each record's full sequence and structure, with 5-mer Jaccard < 0.4.
- GPU conda env `modules/embed_rna_graphs/environment.gpu.yml` and `tests/nextflow.gpu.config`, matching nf-core ribodetector's `task.accelerator` switch.
- `scripts/bump_ginfinity_containers.py` rebuilds the CPU and GPU Wave images for a new `ginfinity` release and pins the URLs in the module `environment*.yml` files and `main.nf` processes.
- GINFINITY embeddings are stored as float16; CPU model inference defaults to full precision, while the GPU profile defaults to float16. `--ginfinity_full_precision` enables full-precision inference explicitly.

### `Changed`

- Index construction and search now use a dedicated process pair for each backend: FAISS, ScaNN, NGT, and cuVS. Each pair owns only its backend environment, container, parameters, and accelerator requirements.
- Index selection is `--index faiss|scann|ngt|cuvs` plus the matching backend-specific index option. `--faiss_index scann` is rejected; use `--index scann`. ScaNN tree size is `--scann_leaves` / `--scann_leaves_to_search`, not `--faiss_nlist` / `--faiss_nprobe`.
- String-valued parameter choices are lowercase, including all `--faiss_index`, `--ngt_index`, `--cuvs_index`, embedding, plotting, and report-theme values.
- R4RNA plots are now one alignment-coordinate SVG per pair (query arcs up, target arcs flipped down, shared x, identity ribbon). The report shows that figure full-width instead of separate query and target diagrams. `alignments.tsv` keeps gapped `query_aligned` / `target_aligned` strings.
- Default git branch is `main` (`manifest.defaultBranch` and the schema `$id`).
- `--plot_max` is now `--plot_max_pairs` (default 25): max alignment pairs **per query**. Each pair plots both partners. `DRAW_R4RNA`, `DRAW_SW`, and `DRAW_RNARTISTCORE` run once per query.
- `ALIGN_CLUSTERS` and plot processes run per query (`SPLIT_CLUSTERS` then `MERGE_ALIGNMENTS` for the shared `alignments.tsv` / `report.html`). Index search stays batched by query shard; `--search_shard_size` sets query records per search task (default: `--shard_size`).
- `report.html` plot panel is a two-column Query | Target grid (one row per RNArtistCore, R4RNA, and alignment plot). Results are paginated: 10 per page by default, or 25 / 50 / 100 / 150.
- README, R4RNA / RNArtistCore plots, and `report.html` use the nf-core palette (green `#24B064`, yellow `#ECDC86`, brown `#3F2B29`, dark green `#396E35`, Bootstrap grays). `--plot_highlight_colour` now defaults to `#24B064`.
- `report.html` is a light theme by default. `--report_theme dark` writes a gray-900 report; a masthead toggle switches themes in the browser.

### `Fixed`

- Test TSVs now emit pair-closed, balanced `.()` full molecules so `ginfinity build-graphs` no longer dies on unmatched brackets.
- Sliced FAISS records and alignments now drop pairs that cross the window, so GINFINITY-SW formatting and structure plots receive a balanced subject.
- `scripts/bump_ginfinity_containers.py` no longer treats a published Wave image as a hard failure when the service reports `succeeded: false`, and it retries SIF HEAD 403s (Python-urllib User-Agent) with a ranged GET.
- `-profile gpu` now passes `--allow-nondeterministic-cuda` to `ginfinity embed-graphs`, which GINFINITY requires on CUDA.
- `-profile gpu` now follows the nf-core ribodetector pattern: `task.accelerator` switches `EMBED_RNA_GRAPHS` to `environment.gpu.yml` (`ginfinity` + `pytorch-gpu=2.6.0` + `cuda-version=12.6`) and the CUDA Wave image. The published ginfinity-only Wave tag is CPU PyTorch, which caused `CUDA was requested but is unavailable`.

### `Dependencies`

- Graph and embed modules use `nicolas.aira::ginfinity=1.2.1` (sliced graphs; float16 embeddings).
- FAISS modules use conda-forge `faiss-cpu=1.10.0` with MKL (`python=3.12`, `numpy=2.2.6`). GPU FAISS uses `pytorch::faiss-gpu=1.10.0` (CUDA 12.1 runtime) so it runs on host drivers that report CUDA 12.1/12.2 (e.g. 535.x). conda-forge `faiss-gpu` 1.10 is CUDA 12.9-only.
- Alignment uses `nicolas.aira::ginfinity-sw=1.0.1`.
- Plotting uses `nicolas.aira::rnartistcore=0.4.6` and `nicolas.aira::r-r4rna=2.0.9`.

### `Deprecated`
