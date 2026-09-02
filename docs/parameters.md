# Parameters

All pipeline parameters are declared in `nextflow.config` and validated
by `nextflow_schema.json` (`nf-schema@2.7.2`). Pass them on the CLI or
in a Nextflow `-params-file`. Do not put them in a `-c` config file.

String choices are lowercase (`flatip`, not `FlatIP`). Booleans accept
`true`/`false`, `1`/`0`, `yes`/`no`.

Flags that do not apply to the chosen `--index` / `--faiss_index` warn
at launch and are ignored.

## Input / output

| Parameter | Default | Description |
|---|---|---|
| `--input` | — | Structures table used to **build** the database. Omit when searching an existing `--database`. |
| `--query` | — | Structures table to **search**. Omit for a build-only run. |
| `--database` | — | Existing window database directory (`index/`). Requires `--query`. Do not combine with `--input`. |
| `--outdir` | — | **Required.** Output directory. |
| `--save_graphs` | `false` | Publish graph safetensors and metadata under `graphs_shards/`. |
| `--save_embeddings` | `false` | Publish per-shard embeddings and manifests under `embeddings_shards/`. The packed `index/embeddings.npz` is unaffected. |
| `--save_windows` | `false` | Publish raw per-shard sliding-window vectors and manifests under `windows_shards/`. |
| `--save_quantized_windows` | `false` | Publish intermediate quantized window shards under `windows_quantized/`. |

Mode rules: [Run modes](usage.md#run-modes).

## Graph construction

| Parameter | Default | Description |
|---|---|---|
| `--shard_size` | `50` | Max input TSV **rows** per `BUILD_RNA_GRAPHS` job. |
| `--id_column` | `transcript_id` | Molecule identifier column. |
| `--sequence_column` | `sequence` | RNA sequence column (`T` is converted to `U`). |
| `--structure_column` | `secondary_structure` | Dot-bracket column, same length as the sequence. |
| `--start_column` | `start` | Optional 0-based half-open window start. |
| `--end_column` | `end` | Optional window end; required if start is present. |
| `--no_slices` | `false` | Ignore start/end and encode full molecules. |
| `--keep_paired_neighbours` | `true` | Keep crossing-pair partners outside a slice in the graph. |
| `--context_hops` | `4` | Expansion around those partners. Does not change embedding length. |
| `--delimiter` | tab | Field delimiter. |

## Embedding

| Parameter | Default | Description |
|---|---|---|
| `--embed_device` | `cpu` | `cpu` or `gpu`. GPU uses the CUDA-backed  |
| `--max_batch_nodes` | `60000` | Max nodes per embedding microbatch. |
| `--max_batch_edges` | `300000` | Max edges per embedding microbatch. |

## Window search

| Parameter | Default | Description |
|---|---|---|
| `--window_size` | `11` | Residues per window. |
| `--window_stride` | `1` | Step between window starts. |
| `--seed_k` | `50` | Neighbours kept per query window after rerank. Increase with DB size. |
| `--candidate_k` | `200` | ANN pool before rerank. Must be ≥ `--seed_k`. |
| `--seed_min_similarity` | `0.8` | Minimum original-window cosine after rerank. Not applied to PQ/OPQ ADC scores when `--exact_rerank false`. |
| `--search_shard_size` | `--shard_size` | Query records per search task. Ignored when `--input` and `--query` are the same file. |
| `--index` | `faiss` | `faiss`, `cagra`, `ivf`, or `hnswlib`. |
| `--quantize` | `none` | `none`, `sq`, `pq`, `opq`. Node-level, before index windows. |
| `--search_device` | `gpu` | CAGRA search: `gpu` or `cpu`. CPU requires the supported CPU-search representation (`--cagra_to_hnsw` for uncompressed/SQ CAGRA, or the CPU ADC walker for PQ-CAGRA).|
| `--cagra_to_hnsw` | `false` | GPU CAGRA build, persist a CPU-searchable graph. |
| `--exact_rerank` | `true` | Original-window rerank. Skipped for FlatIP/FlatL2. Use `--exact_rerank false` to skip. |
| `--exact_rerank_device` | `cpu` | `cpu` or `gpu`. GPU uses CuPy batches and requires a GPU-capable execution environment. |
| `--exact_rerank_batch_size` | `32` | Query windows per rerank batch. |
| `--exact_rerank_candidate_batch_size` | `2048` | Candidate windows per sub-batch. |
| `--exact_rerank_workers` | `0` | CPU rerank workers; `0` = process CPU count. |
| `--hnswlib_rerank` | `false` | Compatibility alias that forces `--exact_rerank`. |

## Quantization (PQ / OPQ / SQ)

Apply when `--quantize` is not `none`.

| Parameter | Default | Description |
|---|---|---|
| `--pq_m` | `16` | Subquantizers per 128-d node. Must divide 128. |
| `--pq_nbits` | `4` | Bits per subquantizer: `4` or `8`. |
| `--pq_sample_size` | `250000` | Max node vectors sampled to fit the codebook. |
| `--pq_niter` | `12` | k-means iterations per PQ subspace. |
| `--opq_iters` | `10` | OPQ Procrustes iterations. |
| `--pq_seed` | `1` | RNG seed for SQ/PQ/OPQ fitting. |

## FAISS (`--index faiss`)

| Parameter | Default | Description |
|---|---|---|
| `--faiss_index` | `flatip` | `flatip`, `flatl2`, `ivfflat`, `hnsw`. |
| `--faiss_device` | `cpu` | FAISS build and search device: `cpu` or `gpu`. GPU is supported for FlatIP and FlatL2 only and requires a GPU-capable execution environment. |
| `--faiss_nlist` | auto | IVF coarse centroids (`ivfflat`). |
| `--faiss_nprobe` | auto | IVF lists visited (`ivfflat`). |
| `--faiss_hnsw_m` | `32` | FAISS HNSW `M`. CPU only. |
| `--faiss_hnsw_ef_construction` | `200` | FAISS HNSW build exploration. |
| `--faiss_hnsw_ef_search` | `64` | FAISS HNSW search exploration. |

## CAGRA (`--index cagra`)

| Parameter | Default | Description |
|---|---|---|
| `--cagra_intermediate_graph_degree` | `128` | Neighbours before prune. |
| `--cagra_graph_degree` | `64` | Neighbours retained. `32` is practical for compact PQ. |
| `--cagra_build_algo` | `nn_descent` | Stock cuVS builder: `nn_descent`, `ivf_pq`, `iterative_cagra_search`, `ace`. Uncompressed/SQ only. |
| `--cagra_nndescent_iterations` | `6` | NN-Descent iterations for **PQ-CAGRA**. |
| `--cagra_itopk_size` | `64` | Intermediate beam; raised to at least `--candidate_k`. |

## cuVS IVF (`--index ivf`)

| Parameter | Default | Description |
|---|---|---|
| `--cuvs_n_lists` | auto from `n` | Coarse clusters. |
| `--cuvs_n_probes` | auto from `n` | Lists searched per query. |

## hnswlib (`--index hnswlib`)

CPU PQ/OPQ only.

| Parameter | Default | Description |
|---|---|---|
| `--hnswlib_m` | `32` | Graph degree. |
| `--hnswlib_ef_construction` | `200` | Build exploration. |
| `--hnswlib_ef_search` | `200` | Search exploration. |
| `--hnswlib_random_seed` | `1` | Graph RNG. |
| `--hnswlib_num_threads` | `0` | `0` = library default. |

## Seed clustering

| Parameter | Default | Description |
|---|---|---|
| `--cluster_span` | `80` | Max gap from a seed to the current cluster box. |
| `--cluster_min_seeds` | `2` | Drop smaller clusters. |
| `--cluster_diagonal_tolerance` | `12` | Allowed diagonal drift. |
| `--cluster_max_diagonal_span` | `96` | Max diagonal breadth; `0` disables. |
| `--cluster_max_seed_rank` | `10` | Ignore worse seed ranks; `0` keeps all. |

## Alignment and EVD

| Parameter | Default | Description |
|---|---|---|
| `--align_pad` | `32` | Extra nucleotides on each side of the cluster crop. |
| `--align_max_cells` | `16777216` | Max DP cells (crop length product). |
| `--align_cpus` | `8` | Threads for `ALIGN_CLUSTERS` and `ESTIMATE_EVD`. One alignment task loads embeddings once. |
| `--align_mu` | `0.3890` | GINFINITY-SW cosine-to-substitution transform location. |
| `--align_sigma` | `1.0000` | GINFINITY-SW transform scale. |
| `--align_gamma` | `1.5616` | GINFINITY-SW transform exponent. |
| `--align_score_min` | `-2.7734` | Lower bound of the cosine-to-substitution score transform. |
| `--align_score_max` | `6.1810` | Upper bound of the cosine-to-substitution score transform. |
| `--align_gap_open` | `1.6042` | GINFINITY-SW affine gap-open cost. |
| `--align_gap_extend` | `0.1923` | GINFINITY-SW affine gap-extension cost. |
| `--align_score_offset` | `0.0449` | Offset applied by the cosine-to-substitution score transform. |
| `--align_max_alignments` | `16` | Max disjoint local HSPs per crop / EVD sample. |
| `--align_min_score` | `0.0` | Minimum HSP score retained. |
| `--align_min_match_count` | `1` | Minimum matched columns retained. |
| `--evd_samples` | `1000` | Reverse-sequence null alignments. |
| `--evd_max_length` | `400` | Max null-sequence length. |
| `--evd_seed` | `1` | Null-sample RNG. |

## Plots and report

| Parameter | Default | Description |
|---|---|---|
| `--plot_backend` | `none` | `none`, `rnartistcore`, `r4rna`, `both`. |
| `--plot_sw` | `false` | Draw crop cosine + SW score matrices. |
| `--plot_max_pairs` | `25` | Max unique query–target pairs plotted **per query**. |
| `--plot_highlight_colour` | `#24B064` | Aligned span, shared pairs, matching bases, SW traceback. |
| `--report_theme` | `light` | `light` or `dark`. |

## Misc

| Parameter | Default | Description |
|---|---|---|
| `--monochrome_logs` | `false` | Disable ANSI colours in Nextflow logs. |
| `--validate_params` | `true` | Schema validation. |
| `--version` | `false` | Print pipeline version and exit. |

Choosing `--index` / `--quantize`: [Window indexes](indexes.md).
