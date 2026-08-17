# nicoaira/ginflow: Documentation

The nicoaira/ginflow documentation is split into the following pages:

- [Usage](usage.md)
  - How to run the pipeline, the structures table, [sliced graphs](usage.md#sliced-graphs) (`start` / `end` windows), and every command-line flag.
- [Window indexes](indexes.md)
- FAISS, ScaNN, NGT, and cuVS, how to choose a library and structure, which flags apply to each, GPU requirements, and unused-parameter warnings.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
- [Pipeline schematic](images/ginflow_metro.svg)
  - Metro-map diagram of the three run modes (build, build-and-search, search an existing database) and optional plotters. Source: [`images/ginflow_metro.mmd`](images/ginflow_metro.mmd). Rebuild with:

    ```bash
    nf-metro render docs/images/ginflow_metro.mmd -o docs/images/ginflow_metro.svg --embed-font
    ```
