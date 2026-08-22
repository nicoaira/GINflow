# nicoaira/ginflow: Documentation

The nicoaira/ginflow documentation is split into the following pages:

- [Usage](usage.md)
  - How to run the pipeline, the structures table, [sliced graphs](usage.md#sliced-graphs) (`start` / `end` windows), and every command-line flag.
- [Window indexes](indexes.md)
  - FAISS, CAGRA, IVF, and hnswlib, node-level `--quantize`, how to choose a library, GPU requirements, and unused-parameter warnings.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
- [Pipeline schematic](images/ginflow_metro.svg)
  - Overview metro map (embeddings → one **Build index** station → search → SW → report). Source: [`images/ginflow_metro.mmd`](images/ginflow_metro.mmd). Rebuild with:

    ```bash
    nf-metro render docs/images/ginflow_metro.mmd -o docs/images/ginflow_metro.svg --embed-font
    ```

- [Index schematic](images/ginflow_metro_index.svg)
  - Quantize × index backends and matching search. Hand-edited SVG (not generated from a `.mmd`).
