# AWS Batch optimization todo

## Blocked prerequisite

- [x] Use the approved S3 bucket/prefix
      `s3://ginflow-nextflow-work-us-east-1/ginflow/`.
- [x] Verify the launch IAM user has bucket location/list/object
      read-write-delete access with a temporary preflight object.
- [x] Add `s3:GetObjectTagging`, `s3:PutObjectTagging`, and
      `s3:AbortMultipartUpload` to the existing bucket policy so Nextflow can
      publish temporary-tagged work objects with S3 `CopyObject`.
- [ ] Confirm the Batch `ecsInstanceRole` has the same object permissions if
      Fusion task logs report an instance-role staging failure.
- [x] Confirm the Fusion/Wave token in `.env` is available to the launch shell.

## Overnight execution

- [x] Submit the GPU smoke command with `-profile awsbatch,fusion,smoke_test`.
- [x] Confirm `BUILD_RNA_GRAPHS` and `GENERATE_WINDOWS` use the CPU Spot queue.
- [x] Confirm GPU embedding, PQ-CAGRA build, and PQ-CAGRA search use
      `ginflow-gpu-queue-virginia`.
- [x] Confirm the Batch job `resourceRequirements` match the new trace and
      that the g4dn GPU path completes successfully.
- [x] Record Spot interruptions, `maxSpotAttempts`, retries, and wait time:
      no interruption/retry; startup/image provisioning dominated duration.
- [x] Start the 30k `--index cagra --quantize pq --pq_m 16 --pq_nbits 4`
      reproduction with the R4RNA plotting backend.
- [x] Record the first PQ-CAGRA Spot host termination and confirm Batch
      requeued the task on a replacement Spot host.
- [ ] Sync the result to
      `/mnt/ssd_samsung/ginflow-benchmarks/pq-cagra-r4rna-30k`.
- [ ] Compare runtime, p95 RSS, CPU, I/O, queue wait, retries, and output
      sizes against the existing 6k benchmark.

## Resource tuning

- [x] Create homogeneous GPU Spot and On-Demand queues for the 16, 24, 48,
      and 96 GB-class tiers.
- [x] Cap every new GPU On-Demand compute environment at 8 vCPUs; use only
      `.2xlarge` instances there.
- [x] Keep larger `.4xlarge`/`.8xlarge` GPU variants available only in Spot
      environments, subject to the 32-vCPU Spot quota.
- [x] Replace the unsupported AWS Batch `g7` tier with the supported `g6` L4
      24-GB-class tier.

- [x] Apply the local-trace CPU reductions for GPU embedding, reranking, and
      EVD estimation in `conf/awsbatch.config`.
- [x] Replace remote-manifest parsing in dynamic closures with workflow `val`
      counts for node quantization, FAISS, and cuVS; retain a safe 8 GB
      quantizer floor after the observed 4 GB OOM.
- [x] Replace the PQ-CAGRA build/search formulas with allocation-derived host
      and GPU estimators plus a 25% safety margin.
- [ ] Measure the 30k PQ-CAGRA GPU peak against the 11-GiB estimate on the
      selected g4dn instance; keep host RSS and GPU memory separate.
- [ ] If larger Spot host RAM is needed, add `g4dn.8xlarge` to the GPU Spot
      compute environment; note that it still has only one 16-GiB T4.
- [ ] Do not submit the default 150k PQ-CAGRA build to `g4dn.2xlarge`:
      its 16-GiB T4 is below the 49-GiB device-memory request.
- [ ] Check whether GPU embedding can use 2 rather than 4 vCPUs without
      reducing throughput.
- [ ] Check whether PQ-CAGRA build is GPU- or input-staging-bound before
      increasing CPU requests.
- [ ] Tune `FIT_NODE_QUANTIZER` CPU count against BLAS utilization and memory.
- [ ] Tune `RERANK_CANDIDATES`, `ESTIMATE_EVD`, and `ALIGN_CLUSTERS` from p95
      `%cpu`, not from the old generic labels.
- [x] Validate the dynamic build/search requests and GPU-tier routing at 30k
      and 150k with standalone Nextflow tasks.
- [ ] Run smaller controlled comparisons for FAISS Flat, FAISS IVF, cuVS
      CAGRA/IVF, and PQ-HNSWLIB.
- [ ] Reject or resize GPU builds above the largest compatible tier; the
      default 150k build routes to the 96-GB-class queue, whose Blackwell
      hosts still require a compatible PQ-CAGRA build and runtime validation.
- [ ] Revisit CPU Spot versus CPU On-Demand routing after the quota decision.

## Cost and operations

- [ ] Obtain authorized On-Demand pricing or an AWS Cost Explorer export for
      exact dollar comparisons; the current IAM user cannot call
      `pricing:GetProducts`.
- [ ] Add an S3 lifecycle rule for abandoned Nextflow work prefixes.
- [ ] Keep result artifacts and reports, but expire temporary work data after
      the agreed retention period.
- [ ] Record AWS Batch queue wait and instance type for each task so cost per
      process can be estimated from actual runtime and pricing.
- [ ] Switch GPU routing to
      `ginflow-gpu-ondemand-queue-virginia` only after quota approval or when
      Spot interruption risk is unacceptable.
