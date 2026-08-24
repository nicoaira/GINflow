# Profiles and hardware

GINflow is configured like an nf-core pipeline (`conf/base.config` plus
named profiles in `nextflow.config`) but is **not** published on
nf-core.

## Container / environment profiles

Pick **one** runtime:

| Profile | Effect |
|---|---|
| `docker` | Docker; `runOptions = -u $(id -u):$(id -g)` |
| `singularity` | Singularity with `autoMounts` |
| `apptainer` | Apptainer with `autoMounts` |
| `podman` / `shifter` / `charliecloud` | Alternative engines |
| `conda` | Per-process `environment.yml`; channels `conda-forge`, `bioconda`, `nicolas.aira`, `rapidsai`, `nvidia` |
| `mamba` | Same as conda with `conda.useMamba` |
| `wave` | Seqera Wave freeze (`conda,container`) + OCI auto-pull for Singularity/Apptainer |
| `arm64` | `process.arch = arm64` and Wave (needed for images) |
| `emulate_amd64` | Docker `--platform=linux/amd64` |

`-profile docker` and `-profile conda` need **no extra local installs**
for index backends: no `pip`, no local CUDA toolkit, no source build.

Custom PQ-CAGRA code is on Anaconda as:

| Package | Channel | Used by |
|---|---|---|
| `pq-cagra-adc` | `nicolas.aira` | GPU PQ-CAGRA build and GPU ADC search |
| `pq-cagra-adc-cpu` | `nicolas.aira` | CPU ADC walker (`--search_device cpu` / `--cagra_to_hnsw`) |

Recipe notes: [`packaging/conda/README.md`](../packaging/conda/README.md).

## GPU

`-profile gpu` does three things:

1. `params.embed_device = 'cuda'`
2. `params.ginfinity_full_precision = false` (float16 inference)
3. NVIDIA runtime flags (`--gpus all` on Docker, `--nv` on
   Singularity/Apptainer) and `accelerator = 1` on `EMBED_RNA_GRAPHS`

That switches `EMBED_RNA_GRAPHS` to `environment.gpu.yml` and the GPU
Wave image (`ginfinity` + `pytorch-gpu=2.6.0` + `cuda-version=12.6`).
Graph construction stays on the CPU image.

FAISS GPU is **opt-in and separate**:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --index faiss \
    --faiss_gpu \
    --faiss_index flatip \
    --input structures.tsv \
    --outdir results
```

`--faiss_gpu` without `-profile gpu` is a pipeline error. GPU FAISS is
exact FlatIP / FlatL2 only (CUDA 12.1 runtime, host drivers that report
12.1 or 12.2, for example NVIDIA 535.x). FAISS IVF and HNSW are
CPU-only; GPU IVF is `--index ivf`; GPU graph is `--index cagra`.

CAGRA and cuVS IVF **build** always need a GPU and `-profile gpu`.
Search a CAGRA graph on CPU after `--cagra_to_hnsw true`:

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --index cagra \
    --cagra_to_hnsw true \
    --input structures.tsv \
    --query queries.tsv \
    --outdir results
```

Later query-only runs can use `--search_device cpu` without a GPU.
PQ/OPQ CAGRA still needs a GPU to **build**; CPU-only PQ **builds** use
`--index hnswlib`.

`--exact_rerank_device cuda` also requires `-profile gpu`.

Who can build and who can search: [Window indexes](indexes.md#devices-who-can-build-and-who-can-search).

## Test profiles

| Profile | Input | Notes |
|---|---|---|
| `smoke_test` | `tests/data/smoke_test_structures.tsv` (10 molecules) | `--shard_size 5`, `--evd_samples 200`, 4 CPU / 15 GB / 1 h cap |
| `test` | `tests/data/test_structures.tsv` (1200 molecules) | `--shard_size 400` |

Both tables are diverse full Rfam molecules (12 families, 5-mer Jaccard
< 0.4). Rebuild with `python3 scripts/build_test_structures.py` if the
source `molecules.jsonl` files change.

`conf/pq_cagra_48gb.config` raises memory for `FIT_NODE_QUANTIZER` and
`BUILD_PQ_CAGRA_INDEX` to 48 GB. Include it with `-c` (that file is
process resources, not pipeline parameters).

## Resource labels

`conf/base.config` (multiplied by `task.attempt`, one retry on
130–145 / 104 / 175–177):

| Label | CPUs | Memory | Time | Typical processes |
|---|---:|---:|---:|---|
| `process_single` | 1 | 6 GB | 4 h | `BUILD_RNA_GRAPHS` |
| `process_low` | 2 | 12 GB | 4 h | light Python |
| `process_medium` | 6 | 36 GB | 8 h | embed, plot, SW |
| `process_high` | 12 | 72 GB | 16 h | heavy index work |
| `process_high_memory` | (base) | 200 GB | — | `FIT_NODE_QUANTIZER` |
| `process_long` | — | — | 20 h | long walltime |

Override on HPC with a site config (`-c`) that only sets `process`
directives, executor, and queues — not `params`.

## Executors

The default executor is `local`. On a cluster, add a config, for
example:

```groovy
process {
    executor = 'slurm'
    queue    = 'cpu'
}
```

Keep GPU processes on a GPU partition; they set `accelerator = 1` when
the gpu profile / `--faiss_gpu` / CAGRA search requires it.

## Software pins (current)

| Component | Pin |
|---|---|
| GINFINITY | `nicolas.aira::ginfinity=1.2.1` |
| GINFINITY-SW | `nicolas.aira::ginfinity-sw=1.2.0` |
| FAISS CPU | conda-forge `faiss-cpu=1.10.0` (Python 3.12, NumPy 2.2.6, MKL) |
| FAISS GPU | `pytorch::faiss-gpu=1.10.0` (CUDA 12.1) |
| RNArtistCore | `nicolas.aira::rnartistcore=0.4.6` |
| R4RNA | `nicolas.aira::r-r4rna=2.0.9` |
| Nextflow | ≥ 25.10.4 |
| nf-schema | 2.7.2 |

After publishing a new `ginfinity` version:

```bash
python3 scripts/bump_ginfinity_containers.py --version X.Y.Z
```

That rebuilds CPU/GPU Wave images and pins URLs in the embed/graph
module `environment*.yml` files and process containers.
