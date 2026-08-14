# GINflow 🧬

---

## Introduction

**nicoaira/ginflow** is a Nextflow pipeline for BLAST-style search over RNA secondary structures. It uses [GINFINITY](https://github.com/nicoaira/GINFINITY) to build graph shards and node embeddings, then slices those embeddings into windows and searches them with FAISS.

The current steps are:

1. Split a structures TSV into shards
2. Build GINFINITY graph shards (`*.safetensors` + `*.json`)
3. Embed each shard (`*.npz` + manifest)
4. Slice sliding windows (default `w=11`, stride 1)
5. Build a reusable FAISS database and/or search it for seeds
6. Cluster nearby seeds into HSPs, align each crop with GINFINITY-SW, and rank by database E-value
7. Optionally plot query/target 2Ds (`--plot_backend rnartistcore`, `r4rna`, or `both`)
8. Write a standalone HTML search report (`report.html`)

---

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Check your setup with `-profile smoke_test` (10 sequences) or `-profile test` (1200 sequences) before running the workflow on actual data.

Prepare a structures TSV:

```tsv
transcript_id	sequence	secondary_structure	rfam_family
rna-1	ACGUACGU	((....))	RF00005
```

Extra columns (e.g `rfam_family` ) are allowed. Then:

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

From a clone of this repo, build the 1200-sequence test database and search it with four molecules shorter than 200 nt, each from a different Rfam family (`tests/data/example_queries.tsv`: RF00001, RF00003, RF01725, RF01852). Those four records are also in the test table, so the top hit for each query should be itself.

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
