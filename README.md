<p align="center">
  <img src="docs/images/ginflow_logo.svg" alt="GINflow" width="150">
</p>

# GINflow

---

## Introduction

**nicoaira/ginflow** is a Nextflow pipeline for BLAST-style search over RNA secondary structures. It uses [GINFINITY](https://github.com/nicoaira/GINFINITY) to build graph shards and node embeddings, then slices those embeddings into windows and searches them with FAISS, [ScaNN](https://github.com/google-research/google-research/tree/master/scann), [NGT](https://github.com/NGT-labs/NGT), [cuVS](https://docs.nvidia.com/cuvs/), or quantized-node [hnswlib](https://github.com/nmslib/hnswlib).

The current steps are:

1. Split a structures TSV into shards
2. Build GINFINITY graph shards (`*.safetensors` + `*.json`)
3. Embed each shard (`*.npz` + manifest)
4. Slice sliding windows (default `w=11`, stride 1)
5. Build a reusable window database (`--index faiss --faiss_index flatip` by default; other FAISS types, `--index scann`, `--index ngt --ngt_index qg|qbg`, GPU-only `--index cuvs --cuvs_index cagra|ivf|ivf-pq`, or `--index hnswlib` with centroid-coded node types) and/or search it for seeds. See [docs/indexes.md](docs/indexes.md). `--faiss_gpu` needs `-profile gpu`; cuVS always needs `-profile gpu`; HNSWLIB uses `-profile conda` or `-profile wave` because it builds the pinned hnswlib 0.8.0 C++ custom-distance driver. The high-recall HNSW profile uses `--node_quantization_k 4096 --hnswlib_candidate_k 5000 --hnswlib_ef_search 5000 --hnswlib_rerank true`.
6. Cluster nearby seeds into HSPs, align each crop with GINFINITY-SW, and rank by database E-value
7. Optionally plot structures (`--plot_backend`: RNArtistCore 2Ds and/or a unified R4RNA alignment arc plot) and/or SW matrices (`--plot_sw`)
8. Write a standalone HTML search report (`report.html`)

## Pipeline overview

<p align="center">
  <img src="docs/images/ginflow_metro.svg" alt="GINflow pipeline metro map" width="100%">
</p>

Three coloured routes match the run modes (`--input`, `--input --query`, `--query --database`). Dashed lines are optional plotters (`--plot_backend`, `--plot_sw`). Open circles are optional steps.

---

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Check your setup with `-profile smoke_test` (10 sequences) or `-profile test` (1200 sequences) before running the workflow on actual data.

Prepare a structures TSV. Three columns are required (`transcript_id`, `sequence`, `secondary_structure`). Extra columns such as `rfam_family` are allowed.

To embed only part of a molecule (a motif, domain, or other interval), add optional `start` and `end` columns. Coordinates are 0-based half-open (`sequence[start:end]`). Empty cells keep the full molecule. Comma-separated lists make several independent subjects from one row:

```tsv
transcript_id	sequence	secondary_structure	rfam_family	start	end
rna-1	ACGUACGUACGUACGU	((....))((....))	RF00005		
rna-2	GGGAAACCCUUUUGGG	......(((....)))	RF00001	9	16
rna-3	…	…	RF00003	34,40	90,95
```

That table yields four subjects: full `rna-1`, `rna-2:9-16`, `rna-3:34-90`, and `rna-3:40-95`. Each is embedded, indexed, and searched on its own. Crossing-pair context is kept by default (`--keep_paired_neighbours`, `--context_hops 4`). Full guide: [Sliced graphs](docs/usage.md#sliced-graphs). A mixed example file is [`tests/data/sliced_structures.tsv`](tests/data/sliced_structures.tsv).

Then:

```bash
# Build a reusable FAISS database
nextflow run nicoaira/ginflow \
   -profile docker \
   --input structures.tsv \
   --outdir <OUTDIR>

# Search an existing database
nextflow run nicoaira/ginflow \
   -profile docker \
   --query queries.tsv \
   --database <OUTDIR>/faiss \
   --outdir <SEARCH_OUTDIR>

# Build and search in one run
nextflow run nicoaira/ginflow \
   -profile docker \
   --input structures.tsv \
   --query queries.tsv \
   --outdir <OUTDIR>
```

### Test database + short queries

From a clone of this repo, build the 1200-sequence test database and search it with four molecules shorter than 200 nt, each from a different Rfam family (`tests/data/example_queries.tsv`: RF00001, RF00003, RF01725, RF01852). Those four records are also in the test table, so the top hit for each query should be itself. To try `start` / `end` windows, use `tests/data/sliced_structures.tsv` as `--input` and/or `--query`.

```bash
# Build the test database
nextflow run . \
    -profile docker,test \
    --outdir testdb

# Search it
nextflow run . \
    -profile docker \
    --query tests/data/example_queries.tsv \
    --database testdb/faiss \
    --outdir test_search
```

Or both steps in one run:

```bash
nextflow run . \
    -profile docker,test \
    --query tests/data/example_queries.tsv \
    --outdir test_search
```

Ranked alignments are in `test_search/alignments.tsv`. Open `test_search/report.html` in a browser for the searchable hit report.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).
