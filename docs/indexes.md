# Window indexes

GINflow searches sliding windows of GINFINITY embeddings. Two words matter:

| Term | Meaning |
|---|---|
| **Database** | The published directory (`--outdir/faiss/`, later `--database`). It holds the search structure plus `windows.tsv`, `meta.json`, packed residue embeddings, and sequences. |
| **Index** | The nearest-neighbour structure *inside* that database. You choose the **library** (`--index`) and its library-specific structure (`--faiss_index`, `--ngt_index`, or `--cuvs_index`). |

`--database` is a path to an existing database. `--index` is not a path; it selects FAISS, ScaNN, NGT, or cuVS.

If you pass a flag that does not apply to the selected library or FAISS structure, the pipeline **warns at launch** and ignores it.

## Choose a library: `--index`

| `--index` | When to use it | GPU | On disk |
|---|---|---|---|
| `faiss` (default) | Exact search, or a specific FAISS approximation (IVF, HNSW, PQ, LSH, SQ). Best documented path. | Optional (`--faiss_gpu`) for some types | `index.faiss` |
| `scann` | Large window collections where ScaNN’s tree + anisotropic hashing is a better speed/recall trade-off than FAISS IVF/PQ. CPU only. | No | `scann/` |
| `ngt` | NGT proximity graphs, including quantized graph variants. CPU only. | No | `ngt/` |
| `cuvs` | GPU-optimized cuVS CAGRA, IVF-Flat, or IVF-PQ indexes. | Required | `cuvs/` |

Windows are 1408-d (`--window_size 11`) and L2-normalized. All libraries search cosine similarity; FAISS/ScaNN use similarity directly, while NGT and cuVS use their cosine-distance output and convert it back before thresholding.

```bash
# Exact FAISS (default)
nextflow run nicoaira/ginflow -profile docker \
    --index faiss --faiss_index FlatIP \
    --input structures.tsv --outdir results

# Approximate FAISS
nextflow run nicoaira/ginflow -profile docker \
    --index faiss --faiss_index IVFFlat --faiss_nlist 4096 --faiss_nprobe 16 \
    --input structures.tsv --outdir results

# ScaNN
nextflow run nicoaira/ginflow -profile docker \
    --index scann \
    --input structures.tsv --query queries.tsv --outdir results
```

`--faiss_index ScaNN` is invalid. FAISS and ScaNN are different libraries; use `--index scann`.

A later `--query --database results/faiss` run reads `meta.json` and picks the matching library. If you also pass `--index` and it disagrees with the database, the pipeline errors.

## Shared search flags

These always apply. They are not library-specific.

| Flag | Default | Role |
|---|---|---|
| `--window_size` | `11` | Nucleotides per window (concatenated 128-d residues → 1408-d) |
| `--window_stride` | `1` | Step between window starts |
| `--seed_k` | `50` | Neighbours kept per query window before `--seed_min_similarity` |
| `--seed_min_similarity` | `0.8` | Minimum cosine to keep a seed |
| `--search_shard_size` | `--shard_size` | Query records per search task |

`FlatL2`, `LSH`, and NGT’s distance-based modes return different raw distances; the pipeline converts them back to a cosine-like score before the threshold.

## NGT: `--index ngt`

[NGT](https://github.com/NGT-labs/NGT) is a CPU-only graph index. The NGT conda package is provided from the `nicolas.aira` channel because it includes the `qbg` command-line tool needed by QG and QBG.

| `--ngt_index` | Structure | Notes |
|---|---|---|
| `NGT` (default) | Regular NGT graph | Uses normalized cosine distance. |
| `QG` | Quantized graph | Builds a regular graph, then quantizes it with `qbg create-qg` / `qbg build-qg`. |
| `QBG` | Quantized blob graph | Uses NGT’s L2 quantization and probing. This is usually the most memory-efficient option, with an approximate recall trade-off. |

Build an NGT database:

```bash
nextflow run nicoaira/ginflow -profile docker \
    --index ngt --ngt_index NGT \
    --input structures.tsv --outdir results
```

Use QG or QBG by changing only the structure:

```bash
nextflow run nicoaira/ginflow -profile docker \
    --index ngt --ngt_index QBG \
    --input structures.tsv --query queries.tsv --outdir results
```

The database remains under `results/faiss/` for compatibility with the rest of the pipeline. NGT’s native files are in `results/faiss/ngt/`; `meta.json` records `backend: "ngt"` and the selected `index_type`. A later `--query --database results/faiss` run detects NGT automatically. NGT is CPU-only and does not use FAISS or ScaNN parameters.

## cuVS: `--index cuvs`

[cuVS](https://docs.nvidia.com/cuvs/) is NVIDIA’s GPU-optimized vector-search library. It requires `-profile gpu` and a visible NVIDIA GPU for both building and searching. The environment uses RAPIDS `cuvs 24.10.00` with CUDA 11.8 and CuPy 13.0.0; the CUDA 11 runtime is intentional so the image works with NVIDIA 535.x drivers as well as newer drivers.

| `--cuvs_index` | Structure | Use when |
|---|---|---|
| `CAGRA` (default) | GPU graph | You want high-throughput approximate search and the graph fits in GPU memory. |
| `IVF` | IVF-Flat | You want coarse partitioning without product-quantizing stored vectors. |
| `IVF-PQ` | IVF with product quantization | You need a smaller GPU-resident index and can trade some recall for compression. |

Build and search with CAGRA:

```bash
nextflow run nicoaira/ginflow -profile docker,gpu \
    --index cuvs --cuvs_index CAGRA \
    --input structures.tsv --query queries.tsv --outdir results
```

For IVF or IVF-PQ, set the coarse-list and probe parameters as needed:

```bash
nextflow run nicoaira/ginflow -profile docker,gpu \
    --index cuvs --cuvs_index IVF-PQ \
    --cuvs_n_lists 4096 --cuvs_n_probes 32 \
    --cuvs_pq_bits 8 --cuvs_pq_dim 256 \
    --input structures.tsv --outdir results
```

The cuVS database is stored in `results/faiss/cuvs/index.bin` and includes the dataset in the serialized index so a later GPU search can load it. `meta.json` records the cuVS type, metric, list/probe settings, and CAGRA graph settings. Search-only runs detect `backend: "cuvs"` from `meta.json`, but still require `-profile gpu`.

## FAISS: `--index faiss`

Set `--faiss_index` to one of the types below. Only the flags listed for that type are used. Other `--faiss_*` and every `--scann_*` flag is unused (warning if you set it).

### Which FAISS type

| `--faiss_index` | Search | Exact? | GPU | Memory | Use when |
|---|---|---|---|---|---|
| `FlatIP` (default) | Inner product | Yes | Yes | Highest | Default. Small/medium databases, or when you need every seed. |
| `FlatL2` | Squared L2 | Yes | Yes | Highest | Same cost as FlatIP; scores are converted back to cosine. |
| `IVFFlat` | Probe `nprobe` of `nlist` cells | No | Yes | High | First approximation to try. Faster than Flat at large *n*. |
| `HNSW` | Graph walk | No | No | High | High recall without IVF training. |
| `IVFPQ` | IVF + product quantizer | No | Yes | Low | Large databases that must stay in RAM. |
| `IVFSQ` | IVF + scalar quantizer | No | Yes | Medium | Compression with simpler codes than PQ. |
| `PQ` | Product quantizer only | No | No | Low | Compression without IVF. Needs ≥ `2^nbits` windows. |
| `SQ` | Scalar quantizer only | No | No | Medium | Light compression, no IVF. |
| `LSH` | Hamming codes | No | No | Low | Very cheap filter; weakest recall. |
| `IVFPQR` | IVF + PQ + refine | No | No | Medium | Extra PQ codes to rerank. L2-only in FAISS 1.10. |

PQ / IVFPQ / IVFPQR need at least `2^nbits` windows (256 when `--faiss_pq_nbits 8`). `IVFPQR` scores are converted back to a cosine-like similarity. The on-disk `index.faiss` is always a CPU index; GPU is used only inside the task.

`--faiss_gpu` needs `-profile gpu`. Supported types: `FlatIP`, `FlatL2`, `IVFFlat`, `IVFPQ`, `IVFSQ`. `HNSW`, `LSH`, `PQ`, `SQ`, and `IVFPQR` stay on CPU.

```bash
nextflow run nicoaira/ginflow \
    -profile docker,gpu \
    --index faiss \
    --faiss_index IVFFlat \
    --faiss_gpu \
    --input structures.tsv \
    --outdir results
```

### FAISS flags by type

| Flag | Default | Applies to | What it does |
|---|---|---|---|
| `--faiss_index` | `FlatIP` | all FAISS | Structure listed above |
| `--faiss_gpu` | `false` | FlatIP, FlatL2, IVFFlat, IVFPQ, IVFSQ | Train/search on GPU; requires `-profile gpu` |
| `--faiss_nlist` | auto `min(4√n, n/39, n)` | IVFFlat, IVFPQ, IVFSQ, IVFPQR | Coarse centroids (IVF cells) |
| `--faiss_nprobe` | auto `min(8, nlist)` | IVFFlat, IVFPQ, IVFSQ, IVFPQR | Cells visited at search. Higher → more recall, slower |
| `--faiss_pq_m` | `16` | PQ, IVFPQ, IVFPQR | PQ subvectors. Must divide `window_dim` (16 divides 1408) |
| `--faiss_pq_nbits` | `8` | PQ, IVFPQ, IVFPQR | Bits per PQ code (`4`, `8`, `12`, `16`) |
| `--faiss_pq_m_refine` | `4` | IVFPQR | Extra PQ codes used to rerank |
| `--faiss_hnsw_m` | `32` | HNSW | Graph neighbours per node |
| `--faiss_hnsw_ef_construction` | `40` | HNSW | Build-time exploration |
| `--faiss_hnsw_ef_search` | `16` | HNSW | Search-time exploration (also applied when querying) |
| `--faiss_lsh_nbits` | `2 × dim` | LSH | Hamming code length |
| `--faiss_sq_type` | `8bit` | SQ, IVFSQ | Scalar width: `8bit`, `6bit`, `4bit`, `fp16` |

Example of a flag that will warn: `--index faiss --faiss_index FlatIP --faiss_nlist 1024` (`nlist` is IVF-only). `--index faiss --scann_reorder 200` also warns.

## ScaNN: `--index scann`

[ScaNN](https://github.com/google-research/google-research/tree/master/scann) (`scann==1.4.2`) is CPU-only (AVX/FMA on x86). The searcher is serialized under `faiss/scann/`. Every `--faiss_*` flag is unused (warning if you set it).

Plan, from the [ScaNN rules of thumb](https://github.com/google-research/google-research/blob/master/scann/docs/algorithms.md):

| Windows in the database | Plan |
|---|---|
| < 20k | Brute-force (exact) |
| < 100k | Asymmetric hashing + exact reorder |
| ≥ 100k, or `--scann_leaves` set | Spherical tree + AH + reorder |

| Flag | Default | Applies when | What it does |
|---|---|---|---|
| `--scann_leaves` | auto `√n` if n ≥ 100k | optional; setting it **forces** a tree | Number of partitions |
| `--scann_leaves_to_search` | auto `min(8, leaves)` | tree is used | Partitions scored per query. Higher → more recall |
| `--scann_reorder` | `100` | AH is used | Exact-reorder pool after AH. Raised to at least `--seed_k` |
| `--scann_ah_dim` | `2` | AH is used | Dimensions per AH block (official recommendation: 2) |
| `--scann_anisotropic` | `0.2` | AH is used | Anisotropic quantization threshold on unit vectors |
| `--scann_soar` | `false` | tree is used | SOAR spilling (dot-product partitioned indexes) |

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --index scann \
    --scann_leaves 4096 \
    --scann_leaves_to_search 32 \
    --scann_reorder 200 \
    --input structures.tsv \
    --query queries.tsv \
    --outdir results
```

Smoke-test and other small databases stay brute-force unless you set `--scann_leaves`.

## Unused-parameter warnings

At launch the pipeline checks flags you actually set (CLI, or a value different from the default). If a flag does not apply to `--index` / `--faiss_index`, it prints a warning and continues.

```
WARN: Unused index parameters for --index faiss / --faiss_index FlatIP (ignored): --faiss_nlist, --scann_reorder. See docs/indexes.md.
```

Hard errors (not warnings):

- `--index` not `faiss`, `scann`, `ngt`, or `cuvs`
- `--faiss_index ScaNN` (use `--index scann`)
- `--index` disagrees with an existing `--database`
- `--faiss_gpu` without `-profile gpu`
- `--faiss_gpu` with `--index scann`, `ngt`, `cuvs`, or a CPU-only FAISS type
- `--index cuvs` without `-profile gpu`
