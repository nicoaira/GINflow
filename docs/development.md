# Development

GINflow is a Nextflow DSL2 pipeline with local modules and Python (plus
a little C++) under `bin/`. Layout follows nf-core conventions
(module `environment.yml`, Wave containers, `nf-schema`) but the
pipeline is **not** an nf-core pipeline: no institutional configs, no
FASTQ samplesheet leftovers, no nf-core sharing wrappers.

## Layout

```
main.nf                 # thin entry: validate, call GINFLOW, publish
workflows/ginflow.nf    # named analysis workflow
workflows/prepare_windows.nf
modules/<name>/         # one process (sometimes two) per directory
  main.nf
  environment.yml
  environment.gpu.yml   # when the process has a CUDA image
  tests/
bin/                    # tools invoked by processes
assets/alignment.json   # GINFINITY-SW scoring parameters
conf/                   # base + test / smoke_test / pq_cagra resources
docs/                   # this wiki
tests/                  # pipeline-level nf-test + pytest
vendor/                 # hnswlib sources, PQ-CAGRA native code
packaging/conda/        # pq-cagra-adc recipes
scripts/                # container bump, test-table rebuild
```

## Conventions

- **Thin entry workflow.** `main.nf` includes `GINFLOW`, joins
  channels, and declares `output { }` / `publish:`. Do not wrap it in
  another named workflow.
- **Params at the edge.** Named workflows take explicit inputs. Do not
  thread `params.outdir` into `GINFLOW`.
- **Workflow outputs, not `publishDir`.** `outputDir = params.outdir`;
  `workflow.output.mode = 'copy'`.
- **Lowercase `channel`.** `Channel.fromPath` is deprecated.
- **Modules live in `modules/<name>/`**, not `modules/local/`.
- Each process has a `label`, `conda "${moduleDir}/environment.yml"`,
  and the nf-core-style container ternary (Singularity HTTPS SIF vs
  Docker Wave tag).
- Extra TSV columns are allowed and ignored. Optional `start`/`end`
  are not extra: they build sliced graphs.
- Do not invent Nextflow config keys. Dump hashes with
  `nextflow run … -dump-hashes`.

## Tests

Pipeline and module tests use [nf-test](https://code.askimed.com/nf-test/).
Default profile in `nf-test.config` is `smoke_test`. GPU tests
(`*.gpu.nf.test`) are ignored unless you run them explicitly.

```bash
nf-test test tests/default.nf.test \
  tests/search.nf.test \
  modules/build_rna_graphs/tests/main.nf.test \
  modules/embed_rna_graphs/tests/main.nf.test \
  modules/generate_windows/tests/main.nf.test \
  modules/build_faiss_index/tests/main.nf.test \
  modules/search_faiss/tests/main.nf.test \
  modules/cluster_seeds/tests/main.nf.test \
  modules/align_clusters/tests/main.nf.test \
  modules/estimate_evd/tests/main.nf.test \
  modules/draw_rnartistcore/tests/main.nf.test \
  modules/draw_r4rna/tests/main.nf.test \
  modules/draw_sw/tests/main.nf.test \
  modules/split_alignments/tests/main.nf.test \
  modules/merge_alignments/tests/main.nf.test \
  modules/write_report/tests/main.nf.test \
  --profile +docker
```

Python unit tests (index constructors, quantization, rerank, report):

```bash
python3 -m pytest tests/test_faiss_indexes.py \
    tests/test_node_quantization.py \
    tests/test_rerank_candidates.py \
    tests/test_merge_alignments.py \
    tests/test_write_report.py \
    tests/test_parameter_values.py
```

FAISS constructors can also run inside the pinned CPU image:

```bash
docker run --rm -v "$PWD":/work -w /work \
  community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e \
  python3 tests/test_faiss_indexes.py
```

Graph `*.safetensors` / `*.json` checksums are **not** stable across
runs. Snapshot `versions.yml` and assert metadata (record count,
identifiers, node count) instead of hashing those artifacts.

Test inputs:

| File | Role |
|---|---|
| `tests/data/smoke_test_structures.tsv` | 10 full molecules |
| `tests/data/test_structures.tsv` | 1200 full molecules |
| `tests/data/example_queries.tsv` | 4 short queries (self-hits in the test DB) |
| `tests/data/sliced_structures.tsv` | Mixed full / single slice / overlapping slices |
| `tests/data/queries_rouskin_structures.tsv` | 512 Rouskin queries |

## Adding a process

1. `modules/<name>/main.nf` + `environment.yml` + `tests/main.nf.test`.
2. `include` it from `workflows/ginflow.nf` (or `prepare_windows.nf`).
3. Emit a channel; publish from `main.nf` `output { }` if users should
   see the file.
4. If you add a parameter: `nextflow.config` `params { }`,
   `nextflow_schema.json`, and a row in [Parameters](parameters.md).
5. If the flag is index-specific, add it to `warn_unused_index_params`
   in `main.nf`.

## Index backends

Each library has a **build** process and a **search** process, with
its own conda env and container:

| Library | Build | Search |
|---|---|---|
| FAISS | `BUILD_FAISS_INDEX` | `SEARCH_FAISS` |
| cuVS CAGRA / IVF | `BUILD_CUVS_INDEX` | `SEARCH_CUVS` |
| PQ-CAGRA | `BUILD_PQ_CAGRA_INDEX` | `SEARCH_PQ_CAGRA` |
| hnswlib | `BUILD_HNSWLIB_INDEX` | `SEARCH_HNSWLIB` |

`main.nf` normalizes `--index` (`cuvs` → `cagra`, `hnsw` → `hnswlib`),
reads `meta.json` on `--database` runs, and errors if the CLI `--index`
disagrees with the database.

## Metro maps

The overview schematic is generated with
[nf-metro](https://github.com/seqeralabs/nf-metro):

```bash
nf-metro render docs/images/ginflow_metro.mmd \
    -o docs/images/ginflow_metro.svg --embed-font
```

`docs/images/ginflow_metro_index.svg` is hand-edited (quantize ×
backend), not generated from a `.mmd`.

## GitHub wiki

Documentation in `docs/` is the source of truth. To publish the same
pages to the GitHub Wiki tab, create the first page once in the GitHub
UI (the `.wiki.git` repo does not exist until then), then:

```bash
python3 scripts/publish_wiki.py
```

That rewrites in-repo links and image paths for Gollum and force-pushes
`https://github.com/nicoaira/ginflow.wiki.git`.

## Credits

Authors: Nicolas Aira, Uciel Chorostecki.
