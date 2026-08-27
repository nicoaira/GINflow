# Getting started

GINflow is a Nextflow pipeline. You do not install the Python tools
yourself: `-profile docker` (or `singularity` / `conda`) pulls pinned
containers.

## Requirements

- [Nextflow](https://www.nextflow.io/) **≥ 25.10.4**
  (`manifest.nextflowVersion`). The pipeline uses workflow outputs, not
  the old `nextflow.preview.output` flag.
- A container runtime:
  - **Docker** — default recommendation
  - **Singularity / Apptainer** — typical on HPC
  - **Conda / Mamba** — no containers; uses the `nicolas.aira` channel
    plus conda-forge / bioconda / rapidsai / nvidia
- Optional: an NVIDIA GPU and `-profile gpu` for CUDA embeddings,
  FAISS Flat on GPU, CAGRA, and cuVS IVF

If you are new to Nextflow, start with the
[nf-core environment setup](https://nf-co.re/docs/get_started/environment_setup/overview).

## Install Nextflow

```bash
# One common path; see nextflow.io for current installers
curl -s https://get.nextflow.io | bash
./nextflow -version
```

Confirm the version is 25.10.4 or newer.

## First run: smoke test

From a clone of this repository, 10 diverse full Rfam molecules
(`tests/data/smoke_test_structures.tsv`):

```bash
nextflow run . \
    -profile docker,smoke_test \
    --outdir smoke_out
```

That **builds** a tiny database. It does not search. Open
`smoke_out/index/` when it finishes: you should see `index.faiss`,
`windows.tsv`, `meta.json`, `embeddings.npz`, `records.tsv`, and
`evd.json`.

## First search: test database + short queries

1200 diverse full molecules (`-profile test`) and four queries under
200 nt from RF00001, RF00003, RF01725, and RF01852
(`tests/data/example_queries.tsv`). Those four records are also in the
test table, so the top hit for each query should be itself.

```bash
# Build
nextflow run . \
    -profile docker,test \
    --outdir testdb

# Search
nextflow run . \
    -profile docker \
    --query tests/data/example_queries.tsv \
    --database testdb/index \
    --outdir test_search
```

Or both in one run:

```bash
nextflow run . \
    -profile docker,test \
    --query tests/data/example_queries.tsv \
    --outdir test_search
```

Then:

1. Open `test_search/report.html` in a browser.
2. Ranked pairs are in `test_search/alignments.tsv`.
3. Nextflow’s own execution report is under `test_search/pipeline_info/`.

To try `start` / `end` windows, use
`tests/data/sliced_structures.tsv` as `--input` and/or `--query`.

## Run from GitHub without cloning

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
    --outdir results
```

Nextflow clones the pipeline into `~/.nextflow/assets/`. Pin a revision
with `-r <tag-or-commit>` once releases exist. The current version is
`1.0.0dev`.

## Your own data

Prepare a TSV with at least `transcript_id`, `sequence`, and
`secondary_structure`. Extra columns are allowed. See [Input](usage.md#input).

```bash
# Build a reusable FAISS database
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
    --outdir results

# Search it later
nextflow run nicoaira/ginflow \
    -profile docker \
    --query queries.tsv \
    --database results/index \
    --outdir search_results
```

Parameters belong on the CLI or in a Nextflow `-params-file`. Do **not**
put pipeline parameters in a `-c` config file; see the
[nf-core note on parameter files](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Profiles you will actually use

| Profile | Purpose |
|---|---|
| `docker` | Pull Wave-frozen images |
| `singularity` / `apptainer` | Same images as SIFs; `autoMounts` on |
| `conda` / `mamba` | Conda envs from each module’s `environment.yml` |
| `gpu` | CUDA embeddings; required for CAGRA/IVF build and `--faiss_gpu` |
| `smoke_test` | 10-sequence input, small resources |
| `test` | 1200-sequence input |
| `wave` | Build/pull via Seqera Wave |

Combine with a comma: `-profile docker,gpu` or `-profile singularity,test`.

CPU embeddings are the default. GPU:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --input structures.tsv \
    --outdir results
```

More: [Profiles and hardware](profiles.md).

## What “done” looks like

| You passed | You should see |
|---|---|
| `--input` only | `index/`; add `--save_graphs`, `--save_embeddings`, `--save_windows`, or `--save_quantized_windows` for intermediate shards |
| `--query` | plus `seeds/seeds.tsv`, `seeds/clusters.tsv`, `alignments.tsv`, `report.html` |
| `--plot_backend` / `--plot_sw` | `plots/` SVGs, inlined in the report |

[Output](output.md) describes every file.

## Resume and logs

```bash
nextflow run . -profile docker,test --outdir testdb -resume
```

`-resume` reuses cached tasks when inputs and scripts are unchanged.
Execution timeline, report, trace, and DAG always go to
`<outdir>/pipeline_info/`.
