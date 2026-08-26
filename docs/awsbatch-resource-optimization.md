# AWS Batch resource optimization

Status: locally validated; the GPU Spot smoke run completed successfully. The
30k AWS benchmark is running with the same Fusion/Wave setup. Its first
PQ-CAGRA host was terminated by Spot and AWS Batch requeued the task; the
replacement attempt is currently running.

Date of this report: 2026-08-26.

## Executive summary

The repository now has an imported, repository-local `awsbatch` profile. It
routes GPU work to the Spot queue while the GPU On-Demand quota is pending,
uses the nf-core process-label convention, records requested and observed
resources in the trace, and scales the large index/quantization requests from
the actual number of nodes or windows.

The original S3 prerequisite was resolved by reusing the existing approved
bucket `s3://ginflow-nextflow-work-us-east-1`. The restricted pipeline user
successfully completed a write/read/delete preflight under a temporary
prefix. The existing bucket-access policy now includes the least-privilege
`s3:GetObjectTagging`, `s3:PutObjectTagging`, and
`s3:AbortMultipartUpload` actions, because Nextflow work files are tagged
`nextflow.io/temporary=true` and S3 publication uses `CopyObject`.
The GPU smoke run completed through AWS Batch and Fusion/Wave:

```text
run: ginflow-gpu-smoke-20260826T071853Z
queue: ginflow-gpu-queue-virginia (Spot)
watch: https://cloud.seqera.io/user/aira-nicolas94/watch/3fo7CCGg7ZqGoP
state: SUCCEEDED; all outputs and pipeline_info artifacts published
```

The 30k AWS run remains in progress. Do not treat the existing local 30k
reference as an AWS result.

## AWS account and capacity checks

The AWS CLI connection is valid:

```text
account: 111823360707
caller:  arn:aws:iam::111823360707:user/nextflow-batch-user
region:  us-east-1
```

The relevant Batch resources were `ENABLED` and `VALID` when checked:

| Workload | Queue | Compute environment | Capacity | Allocation | Max vCPUs |
|---|---|---|---|---|---:|
| CPU default | `ginflow-cpu-queue-virginia` | `ginflow-cpu-compute-env-virginia` | Spot | `SPOT_CAPACITY_OPTIMIZED` | 256 |
| CPU FAISS | `ginflow-cpu-ondemand-queue-virginia` | `ginflow-cpu-ondemand-env-virginia` | On-Demand | `BEST_FIT_PROGRESSIVE` | 256 |
| GPU default | `ginflow-gpu-queue-virginia` | `ginflow-gpu-compute-env-virginia` | Spot | `SPOT_CAPACITY_OPTIMIZED` | 32 |
| GPU future fallback | `ginflow-gpu-ondemand-queue-virginia` | `ginflow-gpu-ondemand-env-virginia` | On-Demand | `BEST_FIT_PROGRESSIVE` | 32 |

The GPU environments offer one GPU on each of:

| Instance | vCPUs | RAM | GPU |
|---|---:|---:|---:|
| `g4dn.xlarge` | 4 | 16 GiB | 1 |
| `g4dn.2xlarge` | 8 | 32 GiB | 1 |
| `g4dn.4xlarge` | 16 | 64 GiB | 1 |

The current GPU Spot and On-Demand compute-environment instance-type lists
contain only these three types; they do not yet admit `g4dn.8xlarge`. The EC2
API reports that `g4dn.8xlarge` would provide 32 vCPUs and 128 GiB host RAM,
but it still has one 16-GiB T4. Adding it would solve host-RAM/quota
constraints, not the 150k default PQ-CAGRA vRAM constraint.

### Homogeneous GPU queues created on 2026-08-26

The following new queues are `ENABLED` and `VALID`. Each queue maps to exactly
one GPU-memory tier. Spot environments have `maxvCpus=32`; every On-Demand
environment has `maxvCpus=8` and only an 8-vCPU `.2xlarge` instance type.

| GPU tier | Spot queue | On-Demand queue | Instance types |
|---|---|---|---|
| T4, 16 GB-class | `ginflow-gpu-t4-16gb-spot-queue` | `ginflow-gpu-t4-16gb-ondemand-queue` | Spot: `g4dn.xlarge`, `.2xlarge`, `.4xlarge`, `.8xlarge`; On-Demand: `.2xlarge` |
| L4, 24 GB-class | `ginflow-gpu-g6-24gb-spot-queue` | `ginflow-gpu-g6-24gb-ondemand-queue` | Spot: `g6.2xlarge`, `.4xlarge`, `.8xlarge`; On-Demand: `.2xlarge` |
| L40S, 48 GB-class | `ginflow-gpu-g6e-48gb-spot-queue` | `ginflow-gpu-g6e-48gb-ondemand-queue` | Spot: `g6e.2xlarge`, `.4xlarge`, `.8xlarge`; On-Demand: `.2xlarge` |
| RTX PRO Server 6000, 96 GB-class | `ginflow-gpu-g7e-96gb-spot-queue` | `ginflow-gpu-g7e-96gb-ondemand-queue` | Spot: `g7e.2xlarge`, `.4xlarge`, `.8xlarge`; On-Demand: `.2xlarge` |

AWS Batch rejected the proposed RTX PRO 4500 `g7` tier as an unsupported
instance-type set for this account/region. The supported L4 `g6` tier is the
available intermediate replacement. The current PQ-CAGRA extension still
needs a CUDA-architecture rebuild before using L4 or L40S hosts.

### GPU instance alternatives

The AWS EC2 API reports these relevant `us-east-1` offerings. GPU memory is
shown in binary GiB from the API's MiB value; AWS product pages commonly label
the same devices as 24, 48, or 96 GB classes.

| Instance | vCPUs | Host RAM | GPU | GPU memory | 150k PQ-CAGRA fit |
|---|---:|---:|---|---:|---|
| `g5.2xlarge` | 8 | 32 GiB | A10G | 22.4 GiB | GPU too small |
| `g6.2xlarge` | 8 | 32 GiB | L4 | 22.4 GiB | GPU too small |
| `g7.2xlarge` | 8 | 32 GiB | RTX PRO 4500 | 32 GiB | GPU too small |
| `g6e.2xlarge` | 8 | 64 GiB | L40S | 44.7 GiB | Below the 49-GiB request |
| `g7e.2xlarge` | 8 | 64 GiB | RTX PRO Server 6000 | 96 GiB | Fits, but likely an expensive overprovision |
| `g6e.8xlarge` | 32 | 256 GiB | L40S | 44.7 GiB | Below the 49-GiB request |
| `g7.8xlarge` | 32 | 128 GiB | RTX PRO 4500 | 32 GiB | GPU too small |
| `g7e.8xlarge` | 32 | 256 GiB | RTX PRO Server 6000 | 96 GiB | Fits with substantial headroom |

The estimator requests 32 GiB of host memory and 49 GiB of device memory for
the current 150k input and default PQ-CAGRA parameters. The raw modeled CUDA
allocations are about 39.0 GiB; the request adds a 25% safety margin. This
routes the build to the 96-GB-class queue. Changing the queue alone is not
sufficient: its compute environment must offer a compatible instance type.

The `pq-cagra-adc` 0.1.1 binary is compiled for CUDA architectures
`70-real`, `75-real`, `80-real`, `86-real`, `89-real`, and `89-virtual`.
This covers T4, RTX 3070, A10G, L4, and L40S with native code; the SM89 PTX
also provides a forward-JIT path for newer compatible drivers. Blackwell
(`g7e`, SM120) still needs a separate CUDA 12.8+ build and runtime validation.
See the official [AWS accelerated-computing
specifications](https://docs.aws.amazon.com/ec2/latest/instancetypes/ac.html)
and [AWS Batch GPU job documentation](https://docs.aws.amazon.com/batch/latest/userguide/gpu-jobs.html).

The Spot price-history query returned the following observations for
`us-east-1a` on 2026-08-25: `$0.3871/h` (`g4dn.xlarge`), `$0.4299/h`
(`g4dn.2xlarge`), and `$0.5416/h` (`g4dn.4xlarge`). These are historical
Spot observations, not a billing quote. The account cannot call
`pricing:GetProducts`, so exact On-Demand comparisons are intentionally not
claimed here.

## Baseline benchmark evidence

The existing local benchmark artifacts are under
`/mnt/ssd_samsung/ginflow-benchmarks/pq-cagra-r4rna-6k`. The command recovered
from shell history was:

```bash
nextflow run main.nf \
  -profile docker,gpu -resume \
  -c conf/pq_cagra_48gb.config \
  -w /mnt/ssd_samsung/ginflow-pq-cagra-wt/6k \
  --input tests/data/rfam_pdb_benchmark/rouskin_sample_6k.cleaned.tsv \
  --query tests/data/rfam_pdb_benchmark/queries.tsv \
  --outdir /mnt/ssd_samsung/ginflow-benchmarks/pq-cagra-r4rna-6k \
  --index cagra --quantize pq --pq_m 16 --pq_nbits 4 \
  --plot_backend r4rna
```

The input contains 3,900 records, 757,397 residues, and 718,397 theoretical
default sliding windows. The completed database metadata reports 710,119
windows and 754,659 nodes after short records are excluded.

The final trace was `pipeline_info/execution_trace_2026-08-22_06-34-06.txt`.
The preparation tasks came from the two earlier traces in the same directory;
cached rows are not treated as measurements.

| Process | Completed tasks | Median duration | Peak RSS | Median CPU |
|---|---:|---:|---:|---:|
| `BUILD_RNA_GRAPHS` | 92 measured across prep traces | 3.8–5.6 s | 203 MB max | 108–144% |
| `EMBED_RNA_GRAPHS` | 91 | 3.7 s | 891 MB p95 | 143% |
| `GENERATE_WINDOWS` | 91 | 2.2 s | 194 MB max | 101% |
| `FIT_NODE_QUANTIZER` | 1 | 17.6 s | 1.20 GB | 638% |
| `GENERATE_QUANTIZED_WINDOWS` | 1 | 3.2 s | 47 MB | 99% |
| `BUILD_PQ_CAGRA_INDEX` | 1 | 20 m 49 s | 713 MB | 101% |
| `SEARCH_PQ_CAGRA` | 1 | 6.8 s | 918 MB | 146% |
| `RERANK_CANDIDATES` | 1 | 41.5 s | 1.20 GB | 157% |
| `ESTIMATE_EVD` | 1 | 8.1 s | 322 MB | 795% |
| `CLUSTER_SEEDS` | 1 | 0.82 s | 14 MB | 98% |
| `ALIGN_CLUSTERS` | 20 | 3.6 s | 316 MB p95 | 288% |
| `MERGE_ALIGNMENTS` | 1 | 0.60 s | 14 MB | 93% |
| `SPLIT_ALIGNMENTS` | 1 | 0.40 s | 14 MB | 84% |
| `DRAW_R4RNA` | 20 | 0.52 s | 14 MB | 99% |
| `WRITE_REPORT` | 1 | 0.71 s | 15 MB | 96% |

The old traces did not record requested `cpus`, `memory`, `time`, or queue.
The effective generic nf-core-style requests were therefore inferred from
the prior configuration: 6 GB for `process_single`, 12 GB for `process_low`,
36 GB for `process_medium`, and 200 GB for `process_high_memory`. Those
defaults were much larger than the measured RSS for most tasks.

### AWS GPU Spot smoke result

The clean retry `ginflow-gpu-smoke-20260826T071853Z` completed successfully
with 18 tasks and published its trace as
`pipeline_info/execution_trace_2026-08-26_09-19-09.txt`. All Batch jobs had
one successful attempt and exit code zero. The GPU work was submitted to
`ginflow-gpu-queue-virginia` with `GPU=1`; the GPU host selected by Spot was a
`g4dn.4xlarge` in `us-east-1a`.

| Process | Tasks | Request | Peak RSS | CPU | Batch queue | Runtime |
|---|---:|---|---:|---:|---|---:|
| `BUILD_RNA_GRAPHS` | 2 | 1 CPU / 2 GB | 189 MB | 111–112% | CPU Spot | 1m55–1m56 |
| `EMBED_RNA_GRAPHS` | 2 | 2 CPU / 4 GB / GPU | 650–856 MB | 80–97% | GPU Spot | 4m46–6m36 |
| `GENERATE_WINDOWS` | 2 | 1 CPU / 2 GB | 1.4–1.5 MB | 87–95% | CPU Spot | 1m28–3m07 |
| `FIT_NODE_QUANTIZER` | 1 | 4 CPU / 4 GB | 1.5 MB | 206% | CPU Spot | 1m37 |
| `GENERATE_QUANTIZED_WINDOWS` | 1 | 1 CPU / 2 GB | 1.5 MB | 38% | CPU Spot | 18.4 s |
| `BUILD_PQ_CAGRA_INDEX` | 1 | 4 CPU / 4 GB / GPU | 147 MB | 107% | GPU Spot | 6m05 |
| `SEARCH_PQ_CAGRA` | 2 | 2 CPU / 4 GB / GPU | 206–214 MB | 209–212% | GPU Spot | 17.6–27.8 s |
| `RERANK_CANDIDATES` | 2 | 2 CPU / 4 GB | 265–273 MB | 157–159% | CPU Spot | 1m47–1m48 |
| `ESTIMATE_EVD_BUILD` | 1 | 4 CPU / 8 GB | 158 MB | 114% | CPU Spot | 1m56 |
| `CLUSTER_SEEDS` | 1 | 2 CPU / 4 GB | 1.5 MB | 82% | CPU Spot | 15.4 s |
| `ALIGN_CLUSTERS` | 1 | 4 CPU / 8 GB | 169 MB | 87% | CPU Spot | 17.3 s |
| `MERGE_ALIGNMENTS` | 1 | 1 CPU / 2 GB | 1.4 MB | 41% | CPU Spot | 14.7 s |
| `WRITE_REPORT` | 1 | 1 CPU / 2 GB | 1.5 MB | 31% | CPU Spot | 17.9 s |

The trace `realtime` was short, while `duration` included AWS Batch startup,
image pull, and task execution. For example, the 30k-scale inference must
not use smoke startup duration as a CPU or memory estimate. The smoke confirms
that the frozen GPU images initialize and that both CAGRA build and search
complete on AWS; the task logs identify the produced index as `PQ_CAGRA`.

### Live AWS 30k status (intermediate)

The requested uncleaned input run is:

```text
run: ginflow-pq-cagra-r4rna-30k-20260826T074256Z
queue: ginflow-gpu-queue-virginia (Spot)
watch: https://cloud.seqera.io/user/aira-nicolas94/watch/4vD7Xs6uBnXOa8
output: s3://ginflow-nextflow-work-us-east-1/ginflow/results/pq-cagra-r4rna-30k-20260826T074256Z
```

The first `BUILD_PQ_CAGRA_INDEX` attempt requested 4 vCPUs, 10 GiB, and one
GPU. AWS reported `Host EC2 (instance i-0859c4764d1218c47) terminated.` and
returned the job to `RUNNABLE`; it was not an OOM or application exit. With
`aws.batch.maxSpotAttempts = 5`, Batch restarted the same job on Spot host
`i-0a368949dd07b4a2f`, a `g4dn.2xlarge` in `us-east-1b`, at 09:13 UTC. The
replacement was observed `RUNNING` with the same 4-vCPU/10-GiB/GPU request.
CloudWatch reported approximately 12.9% average EC2 CPU while the index was
building, about one core on that 8-vCPU host; this is preliminary evidence for
reducing the PQ-CAGRA CPU request after the final trace confirms it.

The first `FIT_NODE_QUANTIZER` attempt in this run used the old remote-path
fallback and was killed at 4 GiB (exit 137). Its automatic retry used 8 GiB
and succeeded. The workflow now passes the node count as a `val` input. The
allocation-derived PQ formula requests 7 GiB for this uncleaned 30k input
before retry scaling.

### Existing local 30k reference

The target directory already contains a completed local reference run:
`/mnt/ssd_samsung/ginflow-benchmarks/pq-cagra-r4rna-30k`. It used
`rouskin_sample_30k.cleaned.tsv` (not the requested uncleaned
`tests/data/rouskin_sample_30k.tsv`), the old `conf/pq_cagra_48gb.config`,
and the local executor. It is useful empirical evidence, but it is not an
AWS Batch result and is not an exact reproduction of the requested input.

Its metadata reports 3,097,344 windows and 3,309,204 nodes. The trace gives
these additional measurements:

| Process | Completed tasks | Median duration | Peak RSS | Median CPU |
|---|---:|---:|---:|---:|
| `BUILD_RNA_GRAPHS` (database) | 424 | 5.3 s | 203 MB p95 | 117% |
| `EMBED_RNA_GRAPHS` (database) | 424 | 3.1 s | 681 MB p95 | 149% |
| `GENERATE_WINDOWS` (database) | 424 | 1.9 s | 103 MB p95 | 100% |
| `FIT_NODE_QUANTIZER` | 1 | 28.4 s | 4.8 GB | 447% |
| `GENERATE_QUANTIZED_WINDOWS` | 1 | 9.8 s | 48 MB | 100% |
| `BUILD_PQ_CAGRA_INDEX` | 1 | 1 h 31 m 43 s | 2.9 GB | 101% |
| `SEARCH_PQ_CAGRA` | 1 | 10.3 s | 2.3 GB | 133% |
| `RERANK_CANDIDATES` | 1 | 52.2 s | 2.3 GB | 145% |
| `ESTIMATE_EVD` | 1 | 14.1 s | 888 MB | 257% |
| `ALIGN_CLUSTERS` | 20 | 11.5 s | 954 MB p95 | 155% |

The requested uncleaned TSV has 4,407,536 residues and 4,115,576 theoretical
windows, so the resource formulas intentionally size the AWS PQ-CAGRA
request above this cleaned-input reference. The AWS run must still measure
GPU VRAM separately from host RSS.

## Resource policy now implemented

The standard labels follow the nf-core convention: `process_single`,
`process_low`, `process_medium`, `process_high`, `process_high_memory`, and
`process_gpu`. The task `tag` remains the per-task identifier used in logs;
AWS `resourceLabels` identify the pipeline and GPU resource class for cost
tracking. The convention follows the labels used by the
[nf-core/rnaseq base configuration](https://github.com/nf-core/rnaseq/blob/master/conf/base.config).

The AWS profile starts with these smaller defaults:

| Label | CPUs | Memory | Time |
|---|---:|---:|---:|
| `process_single` | 1 | 2 GB | 4 h |
| `process_low` | 2 | 4 GB | 4 h |
| `process_medium` | 4 | 8 GB | 8 h |
| `process_high` | 8 | 16 GB | 16 h |
| `process_high_memory` fallback | 4 | 16 GB | 24 h |

Named overrides further reduce small tasks and give the index builders
input-scaled memory. GPU-only CAGRA/CuVS builders have both the standard
`process_gpu` label and the GPU Spot queue. `EMBED_RNA_GRAPHS`, GPU FAISS, GPU
search, and CUDA reranking are selected by their device parameters and are
also routed to the GPU Spot queue.

## Dynamic sizing formulas

The window count is now emitted by `GENERATE_QUANTIZED_WINDOWS` as an explicit
value and passed to PQ-CAGRA/HNSW builds. FAISS and cuVS builders sum
`n_windows` from their window manifests. The node quantizer sums record
lengths from embedding manifests. Every dynamic request is multiplied by
`task.attempt`, preserving the retry-growth behavior used by nf-core.

PQ-CAGRA sizing is implemented in
`conf/modules/build_pq_cagra_index.resources.config` and
`conf/modules/search_pq_cagra.resources.config`. Let `N` be database windows,
`R` residue nodes, `Q` query windows, `W` window size, `D` embedding dimension,
`M` PQ subquantizers, `B` PQ bits, `I` intermediate graph degree, and `G` final
graph degree. The build estimates are:

```text
host bytes = 0.5 GiB + 2*N*W*M + 4*N*G + 128*N + 4*R*D
GPU bytes  = 0.5 GiB + N*(W*M + 12*I + 4*G + 4) + 4*M*(2^B)^2
request    = ceil(1.25 * bytes / GiB), with a 2-GiB host minimum
```

The search estimator includes the serialized index, Python target/query maps,
three overlapping float32 query copies, result tuples, and CPU per-worker
LUT/hash/traversal workspaces. Its GPU estimate includes the resident index,
query vectors, result buffers, hash tables, and padded traversal workspace.
Both add the same 0.5-GiB base and 25% margin. Search values below assume one
1,000-window query shard and the default parameters.

| Dataset | Nodes | Windows | Build host | Build GPU | Search host | Search GPU |
|---|---:|---:|---:|---:|---:|---:|
| Existing 6k metadata | 0.755 M | 0.710 M | 2 GB | 3 GiB | 2 GB | 2 GiB |
| Existing cleaned 30k metadata | 3.309 M | 3.097 M | 6 GB | 8 GiB | 3 GB | 3 GiB |
| `rouskin_sample_30k.tsv` | 4.408 M | 4.116 M | 7 GB | 11 GiB | 4 GB | 4 GiB |
| `rouskin_sample_150k.tsv` | 22.467 M | 20.967 M | 32 GB | 49 GiB | 16 GB | 12 GiB |

The following processes now also have allocation-derived module configs:

- `FIT_NODE_QUANTIZER` models the overlapping normalized node matrices,
  sample, k-means distance workspaces, assignments, and codes. SQ, PQ, and
  OPQ have separate peaks.
- `GENERATE_QUANTIZED_WINDOWS` uses the largest input shard and record. PQ/OPQ
  includes the temporary bit expansion used by `pack_pq_codes`; SQ includes
  float32 reconstruction and the largest accumulated window shard.
- `BUILD_CUVS_INDEX` models host and device memory separately. CAGRA includes
  the two float32 window matrices, NN-descent graph/distance buffers, the
  optimized graph, and the final device dataset copy. IVF includes the input,
  populated lists, and the default 50% k-means train set.
- `SEARCH_CUVS` includes the resident dataset/index, target map, largest query
  record, accumulated hits, and CAGRA/IVF traversal workspaces.
- `RERANK_CANDIDATES` has distinct CPU and CUDA paths. Host RAM includes the
  database target map, residue embeddings, all query/candidate metadata, and
  one candidate materialization per active CPU worker. CUDA VRAM contains one
  bounded query/candidate batch rather than the complete database.
- `CLUSTER_SEEDS` scales from the exact merged seed-row count.
- `EMBED_RNA_GRAPHS` uses each graph shard's node count, edge count, tensor
  size, the configured microbatch bounds, precision, and CPU/GPU execution
  path.

Every estimator adds a fixed runtime allowance, 25% headroom, and
`task.attempt` retry growth for host RAM. CUDA estimates are exposed as
`task.ext.gpu_memory_gb`; they do not replace Nextflow's host `memory`
directive.

For the cleaned local 30k metadata (3.309 M nodes, 3.097 M windows, 9,087
query windows, default CAGRA and rerank settings), the new estimates are:

| Process | Host RAM | GPU VRAM |
|---|---:|---:|
| PQ quantizer | 5 GiB | — |
| cuVS CAGRA build, uncompressed | 51 GiB | 43 GiB |
| cuVS CAGRA search, uncompressed | 23 GiB | 22 GiB |
| Exact rerank, 2 CPU workers | 4 GiB | — |
| Exact rerank, CUDA | 4 GiB | 2 GiB |
| Cluster 12,181 seeds | 1 GiB | — |

These cuVS values are intentionally much larger than PQ-CAGRA because cuVS
indexes the 1,408-dimensional float32 windows. The 43-GiB CAGRA build estimate
routes to the 48-GB-class queue; a 16-GiB T4 cannot build that uncompressed
30k index.

The 30k PQ-CAGRA build routes to the 16-GB-class T4 queue. The 150k build
routes to the 96-GB-class queue because the 48-GB-class L40S exposes 44.7 GiB,
below the 49-GiB request. Search routing is independent: CPU search stays on
the CPU queue, while GPU search uses the smallest GPU-memory tier that meets
the query-shard-dependent estimate.

Consequently, the approved 8-vCPU On-Demand capacity does not make the 150k
default PQ-CAGRA run viable on `g4dn.2xlarge`: its 16-GiB T4 is below the
49-GiB device-memory request. The run requires a compatible 96-GB-class GPU,
or a memory-reduced/chunked CAGRA implementation with
explicitly revalidated accuracy and performance.

For the default PQ-CAGRA build (`window_size=11`, `pq_m=16`, `pq_nbits=4`,
`intermediate_graph_degree=128`, `graph_degree=64`, and NN-Descent enabled for
large inputs), the dominant temporary device buffers are the uint8 codes,
three intermediate-degree uint32 arrays, the graph-degree uint32 graph, and
reverse-edge counts:

```text
N * (window_size*pq_m + 3*intermediate_graph_degree*4
    + graph_degree*4 + 4) bytes
```

This term is approximately 1.30 GiB for 0.710 M windows, 7.55 GiB for the
requested 4.116 M-window 30k input, and 38.51 GiB for 20.967 M windows at
150k. The estimator adds the similarity table, a 0.5-GiB base allowance, and
25% headroom before selecting a GPU tier. These are allocation-derived
estimates, not measured `nvidia-smi` peaks.

The completed 30k AWS trace should be used to validate the modeled peak and
adjust the documented safety margin if allocator or library overhead differs.

The current profile also applies the local CPU evidence before the AWS run:
embedding and reranking are capped at 2 CPUs, while their host memory is now
input-scaled, and EVD estimation uses at most 4 CPUs. GPU index builders retain
4-CPU requests because their host-side staging/graph work has not yet been
measured on Batch. These CPU requests remain provisional until the AWS trace
reports queue wait, utilization, and GPU memory.

The first 30k AWS quantizer attempt exposed why remote manifests must not be
parsed from a process-selector closure: Fusion supplied a symlink metadata
file, the count became zero, and the 4 GB fallback was OOM-killed (exit 137).
The workflow now parses manifest counts before task submission and passes
`n_nodes`/`n_windows` as `val` inputs. The retry at 8 GB completed, and the
new PQ allocation model resolves this uncleaned 30k input to 7 GB before the
attempt multiplier.

## Validation performed

- `nextflow config . -profile awsbatch` resolves `us-east-1`, the CPU default
  queue, empty `aws.batch.cliPath`, `maxSpotAttempts=5`, and `/tmp` volumes.
- `-profile awsbatch,gpu` routes PQ-CAGRA build and GPU search to the smallest
  homogeneous GPU-memory queue that satisfies `task.ext.gpu_memory_gb`; CPU
  search stays on the CPU Spot queue.
- Standalone Nextflow tasks resolved the allocation-derived values at the
  requested 30k and 150k counts: build 7/11 GiB and 32/49 GiB host/GPU;
  1,000-window search shards 4/4 GiB and 16/12 GiB host/GPU. The same tasks
  verified 16-GB- and 96-GB-class build routing and CPU/GPU search routing.
- A full local PQ-CAGRA CPU-search stub run completed through report writing,
  confirming that the added scalar sizing inputs are wired across the named
  workflow and remain compatible with process stubs.
- A full local CPU HNSW PQ run completed with the AWS profile and local
  executor. Its trace measured 4 GB for node-PQ fitting (13.6 MB RSS) and
  6 GB for HNSW build (375 MB RSS).
- A full local CPU FAISS smoke build completed with the AWS profile and local
  executor. FAISS requested 8 GB and used 14.7 MB RSS for the small input.
- The first live 30k node-quantizer attempt requested 4 GB because the
  process-selector closure could not parse a Fusion manifest path; AWS Batch
  reported `OutOfMemoryError`/exit 137. Its automatic retry requested 8 GB
  and completed. The count is now passed as a workflow `val`, so this fallback
  is no longer used for normal AWS runs.
- After the remote-manifest fix, the six search workflow tests passed, along
  with the FAISS, cuVS, PQ-CAGRA, HNSWLIB, and quantizer module tests.
- The local 30k reference was used to reduce GPU embedding from 4 to 2 CPUs,
  reranking from 4 to 2 CPUs, and EVD estimation from the default 8 to a
  4-CPU cap.
- The AWS CLI preflight against `ginflow-nextflow-work-us-east-1` succeeded
  for bucket location and temporary object write/read/delete.
- The launch user initially failed S3 `CopyObject` because the work objects
  carry `nextflow.io/temporary=true`; adding `GetObjectTagging`,
  `PutObjectTagging`, and `AbortMultipartUpload` to the existing bucket
  policy fixed the publisher preflight without exceeding the IAM inline-policy
  aggregate limit.
- The clean GPU Spot smoke run succeeded with GPU embedding, PQ-CAGRA build,
  PQ-CAGRA search, downstream CPU work, and all S3 outputs published.
- Batch selected `g4dn.4xlarge` Spot hosts for the one-GPU tasks; no Spot
  interruption or retry occurred in the successful run.
- Targeted nf-test modules passed: `BUILD_PQ_CAGRA_INDEX`,
  `BUILD_HNSWLIB_INDEX`, and `GENERATE_QUANTIZED_WINDOWS`.
- A real local CUDA smoke attempt reached the GPU image but failed with
  `CUDA was requested but is unavailable`; the host has an RTX 3070, but the
  local Docker runtime did not expose it. AWS GPU validation remains required.

## AWS launch prerequisites

Use Fusion for the AWS run so task staging does not depend on a host AWS CLI:

1. Provide a bucket, for example `s3://ACCOUNT-APPROVED-BUCKET/ginflow/`.
2. Grant the launch user bucket-location, bucket-list, and object
   read/write/delete permissions for the selected prefix.
3. Grant the Batch `ecsInstanceRole` equivalent object permissions. Fusion and
   ordinary Batch task staging use the Batch instance role inside the task
   environment.
4. Keep `TOWER_ACCESS_TOKEN` in `.env`; do not commit it. Load it only in the
   launch shell.

The current profile deliberately leaves `aws.batch.cliPath` empty. With
Fusion, `aws.batch.cliPath` is disabled by the `fusion` profile. Without
Fusion, the task runtime must contain an AWS CLI or `--awscli` must point to a
valid host-AMI path.

## Planned run commands

Replace `BUCKET` with the approved bucket and keep the generated run IDs:

```bash
set -a
. ./.env
set +a

RUN_ID=$(date -u +%Y%m%dT%H%M%SZ)
nextflow run . \
  -profile awsbatch,fusion,gpu,smoke_test \
  -w s3://BUCKET/ginflow/work/gpu-smoke-${RUN_ID} \
  --input tests/data/smoke_test_structures.tsv \
  --query tests/data/smoke_test_structures.tsv \
  --index cagra --quantize pq --pq_m 16 --pq_nbits 4 \
  --outdir s3://BUCKET/ginflow/results/gpu-smoke-${RUN_ID} \
  -with-report -with-trace -with-timeline -with-dag \
  -name ginflow-gpu-smoke-${RUN_ID}
```

Then run the 30k reproduction with the same query set used by the prior 6k
benchmark:

```bash
RUN_ID=$(date -u +%Y%m%dT%H%M%SZ)
nextflow run . \
  -profile awsbatch,fusion,gpu \
  -w s3://BUCKET/ginflow/work/pq-cagra-r4rna-30k-${RUN_ID} \
  --input tests/data/rouskin_sample_30k.tsv \
  --query tests/data/rfam_pdb_benchmark/queries.tsv \
  --index cagra --quantize pq --pq_m 16 --pq_nbits 4 \
  --plot_backend r4rna \
  --outdir s3://BUCKET/ginflow/results/pq-cagra-r4rna-30k-${RUN_ID} \
  -with-report -with-trace -with-timeline -with-dag \
  -name ginflow-pq-cagra-r4rna-30k-${RUN_ID}
```

After completion, sync the result and `pipeline_info/` to
`/mnt/ssd_samsung/ginflow-benchmarks/pq-cagra-r4rna-30k` for comparison with
the existing 6k artifacts.

## Measurement method for the final report

For each AWS run, retain:

- Nextflow trace: requested `cpus`, `memory`, `time`, `queue`, `tag`,
  `peak_rss`, `%cpu`, I/O, duration, status, exit, and attempt.
- Nextflow report, timeline, DAG, and `.nextflow.log`.
- AWS Batch `describe-jobs` records for every `native_id`, including resource
  requirements, attempts, status reason, exit code, and log stream.
- CloudWatch task logs for CUDA initialization, OOMs, Spot interruptions, and
  staging failures.

For each process, compare the request with p95 peak RSS and CPU utilization.
The next revision should set memory to approximately p95 RSS plus 30–50%
for normal tasks, retaining a larger floor for GPU graph construction and
retrying only known transient/Spot failures. CPU requests should follow the
measured p95 parallelism rather than the old 6/12-core generic tiers.
