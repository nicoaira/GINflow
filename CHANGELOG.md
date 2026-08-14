# nicoaira/ginflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]


### `Added`

- FAISS seed search: `GENERATE_WINDOWS` slices concatenated `w=11` windows, `BUILD_FAISS_INDEX` writes a reusable `IndexFlatIP` database, and `SEARCH_FAISS` returns seeds above `--seed_min_similarity`. Modes are inferred from `--input` / `--query` / `--database`.
- Seed clustering (`CLUSTER_SEEDS`) and GINFINITY-SW alignment (`ALIGN_CLUSTERS`) on a padded crop of each cluster. The FAISS directory now also packs residue embeddings and sequences so query-only runs can align.
- Database E-values: `ESTIMATE_EVD` fits Karlin–Altschul λ, K from reverse-sequence GINFINITY-SW scores; `ALIGN_CLUSTERS` ranks hits by ascending E = K m N exp(−λS).
- Optional alignment plots: `DRAW_RNARTISTCORE` and `DRAW_R4RNA` (`--plot_backend`). Conda packages `nicolas.aira::rnartistcore=0.4.6` (OpenJDK 17–21) and `nicolas.aira::r-r4rna=2.0.9`, with Wave-frozen Docker and Singularity images. Each process draws molecules in parallel with `task.cpus` (6 on `process_medium`).
- Search HTML report: `WRITE_REPORT` writes a self-contained `report.html` on every query run (ranked hits, span rails, inlined plots).
- Diverse Rfam test tables: 10-sequence `-profile smoke_test` (`tests/data/smoke_test_structures.tsv`) and 1200-sequence `-profile test` (`tests/data/test_structures.tsv`), rebuilt from 12 Rfam families using each record's full sequence and structure, with 5-mer Jaccard < 0.4.
- GPU conda env `modules/embed_rna_graphs/environment.gpu.yml` and `tests/nextflow.gpu.config`, matching nf-core ribodetector's `task.accelerator` switch.
- `scripts/bump_ginfinity_containers.py` rebuilds the CPU and GPU Wave images for a new `ginfinity` release and pins the URLs in the module `environment*.yml` files and `main.nf` processes.

### `Changed`

- Default git branch is `main` (`manifest.defaultBranch` and the schema `$id`).

### `Fixed`

- Test TSVs now emit pair-closed, balanced `.()` full molecules so `ginfinity build-graphs` no longer dies on unmatched brackets.
- `-profile gpu` now passes `--allow-nondeterministic-cuda` to `ginfinity embed-graphs`, which GINFINITY requires on CUDA.
- `-profile gpu` now follows the nf-core ribodetector pattern: `task.accelerator` switches `EMBED_RNA_GRAPHS` to `environment.gpu.yml` (`ginfinity` + `pytorch-gpu=2.6.0` + `cuda-version=12.6`) and the CUDA Wave image. The published ginfinity-only Wave tag is CPU PyTorch, which caused `CUDA was requested but is unavailable`.

### `Dependencies`

- FAISS modules use conda-forge `faiss-cpu=1.10.0` with MKL (`python=3.12`, `numpy=2.2.6`).
- Alignment uses `nicolas.aira::ginfinity-sw=1.0.1`.
- Plotting uses `nicolas.aira::rnartistcore=0.4.6` and `nicolas.aira::r-r4rna=2.0.9`.

### `Deprecated`
