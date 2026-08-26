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
| `fusion` | Fusion S3 filesystem with Wave task provisioning; requires `TOWER_ACCESS_TOKEN` |
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

`ESTIMATE_EVD` and `ALIGN_CLUSTERS` are sized dynamically from their
database/query payloads, alignment crop/cell limits, and worker count. Their
selectors live in [`conf/modules/estimate_evd.resources.config`](../conf/modules/estimate_evd.resources.config)
and [`conf/modules/align_clusters.resources.config`](../conf/modules/align_clusters.resources.config),
which are included by `nextflow.config`. They are not `process_high` tasks,
so profile configs should not replace their memory with a fixed selector.

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

## External nf-core profiles and AWS Batch

The pipeline loads custom profiles from [`nf-core/configs`](https://github.com/nf-core/configs)
using the same `custom_config_version` / `custom_config_base` mechanism as
nf-core pipelines. The default is the `master` branch. Set
`NXF_OFFLINE=true` to skip the remote include; local custom-config paths still
work in offline mode.

This makes the shared `awsbatch` profile available, together with other
institution profiles. The repository-local `conf/awsbatch.config` is loaded
after the shared profile and supplies the ginflow queue routing:

- ordinary CPU processes use `ginflow-cpu-queue-virginia`;
- CPU FAISS uses `ginflow-cpu-ondemand-queue-virginia`;
- GPU-only CAGRA/CuVS builds use the GPU Spot tier selected by the resource
  policy once automatic tier routing is enabled;
- GPU embedding, GPU FAISS, and GPU search can target the corresponding GPU
  Spot tier until the GPU On-Demand quota is approved.

The AWS Batch GPU queues are homogeneous by GPU-memory tier:

| Tier | Spot queue | On-Demand queue | Spot instances | On-Demand instance |
|---|---|---|---|---|
| 16 GB-class T4 | `ginflow-gpu-t4-16gb-spot-queue` | `ginflow-gpu-t4-16gb-ondemand-queue` | `g4dn.xlarge`–`g4dn.8xlarge` | `g4dn.2xlarge` |
| 24 GB-class L4 | `ginflow-gpu-g6-24gb-spot-queue` | `ginflow-gpu-g6-24gb-ondemand-queue` | `g6.2xlarge`–`g6.8xlarge` | `g6.2xlarge` |
| 48 GB-class L40S | `ginflow-gpu-g6e-48gb-spot-queue` | `ginflow-gpu-g6e-48gb-ondemand-queue` | `g6e.2xlarge`–`g6e.8xlarge` | `g6e.2xlarge` |
| 96 GB-class RTX PRO Server 6000 | `ginflow-gpu-g7e-96gb-spot-queue` | `ginflow-gpu-g7e-96gb-ondemand-queue` | `g7e.2xlarge`–`g7e.8xlarge` | `g7e.2xlarge` |

All On-Demand environments have `maxvCpus=8` and contain only an 8-vCPU
instance type. AWS Batch rejected the RTX PRO 4500 `g7` family in this
region, so the intermediate tier uses the supported L4 `g6` family instead.

The names ending in `-env-virginia` are compute environments, not queue
names. The tier-specific queue parameters are available for the forthcoming
memory-aware selector; the legacy `aws_gpu_queue` and
`aws_gpu_ondemand_queue` aliases currently point to the T4 queues. Override
any of these parameters when using different queues.

The measured resource baseline, dynamic sizing formulas, and the pending AWS
run checklist are in [the AWS Batch resource report](awsbatch-resource-optimization.md)
and [todo list](awsbatch-todo.md).

```bash
nextflow config . -profile awsbatch
```

AWS Batch needs an S3 work directory. The profile defaults to `us-east-1` and
does not mount a host AWS CLI path. Override the queue parameters on the
command line when needed:

```bash
nextflow run . \
    -profile awsbatch,smoke_test \
    -w s3://<bucket>/ginflow/work/cpu-${RUN_ID} \
    --aws_cpu_queue <cpu-queue> \
    --awsregion us-east-1 \
    --input tests/data/smoke_test_structures.tsv \
    --query tests/data/smoke_test_structures.tsv \
    --outdir s3://<bucket>/ginflow/results/cpu-${RUN_ID}
```

For a GPU smoke test, add `gpu --faiss_gpu`; the profile routes GPU work to
the GPU Spot queue while the On-Demand quota is pending:

```bash
nextflow run . \
    -profile awsbatch,gpu,smoke_test \
    -w s3://<bucket>/ginflow/work/gpu-${RUN_ID} \
    --aws_gpu_queue ginflow-gpu-t4-16gb-spot-queue \
    --awsregion us-east-1 \
    --faiss_gpu \
    --input tests/data/smoke_test_structures.tsv \
    --query tests/data/smoke_test_structures.tsv \
    --outdir s3://<bucket>/ginflow/results/gpu-${RUN_ID}
```

For non-Fusion staging, the task runtime must provide the AWS CLI because
`aws.batch.cliPath` is empty. Set `--awscli` if the AWS CLI is installed at a
known host-AMI path, or use the Fusion-based AWS setup below.

### Fusion/Wave S3 staging

The opt-in `fusion` profile replaces AWS CLI file staging with the Fusion
filesystem and enables Wave for the task runtime. Combine it with the imported
`awsbatch` profile, in this order:

```bash
export TOWER_ACCESS_TOKEN=<seqera-platform-token>

nextflow run . \
    -profile awsbatch,fusion,smoke_test \
    -w s3://<bucket>/ginflow/work/fusion-${RUN_ID} \
    --aws_cpu_queue <cpu-queue> \
    --awsregion us-east-1 \
    --input tests/data/smoke_test_structures.tsv \
    --query tests/data/smoke_test_structures.tsv \\
    --outdir s3://<bucket>/ginflow/results/fusion-${RUN_ID}
```

The token is read from the launch environment and is not stored in the
repository. Fusion still uses the AWS Batch instance/task role for the S3
bucket; grant it `s3:ListBucket` and object read/write/delete permissions for
the work bucket. The hosted Wave service is free, and Seqera Cloud currently
includes a monthly Fusion free tier of up to 100 TB of throughput; Fusion
stops accepting new runs after the quota is reached. AWS Batch and S3 charges
remain separate from Seqera licensing. See the [Fusion licensing](https://docs.seqera.io/fusion/licensing)
and [Wave](https://docs.seqera.io/wave) documentation for current limits.

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
