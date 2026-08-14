# GINflow 🧬

---

## Introduction

**nicoaira/ginflow** is a Nextflow pipeline for BLAST-style search over RNA secondary structures. It uses [GINFINITY](https://github.com/nicoaira/GINFINITY) to build graph shards and node embeddings from a structures table. FAISS-based seeding and optional banded Smith–Waterman alignments will be added later.

The current steps are:

1. Split the input structures TSV into shards
2. Build GINFINITY graph shards (`*.safetensors` + `*.json`)
3. Embed each shard (`*.npz` + manifest)

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
nextflow run nicoaira/ginflow \
   -profile docker \
   --input structures.tsv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).
