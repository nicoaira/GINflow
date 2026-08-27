# GINflow documentation

**GINflow** ([nicoaira/ginflow](https://github.com/nicoaira/ginflow))
is a Nextflow pipeline for BLAST-style search over RNA secondary
structures. It uses [GINFINITY](https://github.com/nicoaira/GINFINITY)
to turn sequence + dot-bracket into per-nucleotide graph embeddings,
indexes sliding windows of those embeddings, clusters seed hits into
HSPs, aligns each crop with
[GINFINITY-SW](https://github.com/nicoaira/GINFINITY-SW), and ranks
pairs by Karlin–Altschul E-values.

<p align="center">
  <img src="images/ginflow_metro.svg" alt="GINflow pipeline metro map">
</p>

## Start here

| Page | What it covers |
|---|---|
| [Getting started](getting-started.md) | Nextflow install, smoke test, first search, profiles |
| [How it works](how-it-works.md) | Graphs → embeddings → windows → index → seeds → SW → E-values |
| [Usage](usage.md) | Structures table, [sliced graphs](usage.md#sliced-graphs), run modes, common flags |

## Using the pipeline

| Page | What it covers |
|---|---|
| [Window indexes](indexes.md) | FAISS, CAGRA, IVF, hnswlib, `--quantize`, GPU vs CPU |
| [Clustering and alignment](clustering-and-alignment.md) | Seed clustering, GINFINITY-SW crops, pair collapse |
| [E-values](statistics.md) | Reverse-sequence null, λ and K, how to read E and bits |
| [Output](output.md) | Directory layout and how to interpret each file |
| [Plotting and report](plotting.md) | `report.html`, RNArtistCore, R4RNA, SW matrices |
| [Parameters](parameters.md) | Every CLI flag with defaults |
| [Profiles and hardware](profiles.md) | Docker / Singularity / conda, GPU, resources |

## Reference

| Page | What it covers |
|---|---|
| [Development](development.md) | Repo layout, tests, adding a process, metro maps |
| [FAQ](faq.md) | Launch errors, OOM, empty hits, windows vs slices |

## Schematics

- [Pipeline metro map](images/ginflow_metro.svg) — embeddings → one **Build index** station → search → SW → report. Source: [`images/ginflow_metro.mmd`](images/ginflow_metro.mmd). Rebuild with:

  ```bash
  nf-metro render docs/images/ginflow_metro.mmd -o docs/images/ginflow_metro.svg --embed-font
  ```

- [Index schematic](images/ginflow_metro_index.svg) — quantize × index backends and matching search. Hand-edited SVG.

## Quick command

```bash
# Build
nextflow run nicoaira/ginflow \
    -profile docker \
    --input structures.tsv \
    --outdir results

# Search
nextflow run nicoaira/ginflow \
    -profile docker \
    --query queries.tsv \
    --database results/index \
    --outdir search_results
```

Requires Nextflow ≥ 25.10.4. Parameters go on the CLI or in
`-params-file`, not in a `-c` config file.
