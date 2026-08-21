# CAGRA + node-level PQ with shared-memory ADC

GPU approximate nearest neighbour index for GINflow sliding windows.

## Layout

```text
original 128-d node embeddings
        |
        +--> product quantizer (M subquantizers, nbits per subcode)
                |
                +--> packed/unpacked node codes
                        |
                        +--> sliding windows of node codes
                        +--> CAGRA graph (SDC / reconstructed inner product)
                        +--> GPU search with shared-memory ADC LUT
```

Search keeps the query in the original node space.  For each query window the
kernel writes a half-precision LUT of shape `(window_size * M, 2**nbits)` into
dynamic shared memory, then scores graph neighbours as a sum of LUT lookups.

Default `window_size=11`, `M=16`, `nbits=4` uses a 5.6 KiB LUT.  `M=16`,
`nbits=8` uses 90 KiB and still fits Ampere's 96 KiB opt-in shared memory.
The itopk / candidate buffer is in global memory so pools of 1,000–5,000
candidates do not compete with the LUT.

The graph is built with reconstructed inner-product SDC (equivalent to the
codebook lookup table).  `n ≤ 16,384` uses brute-force kNN; larger sets use
NN-Descent.

CPU search of a saved graph (no NVIDIA HNSW conversion):

```python
index = pq.load_index("cagra.index")  # host only
labels, dists = pq.search(index, queries, k=1000, device="cpu", num_threads=16)
```

Measured results, recommended serving profiles, and reproduction commands:
[docs/pq-cagra-adc-research.md](../../docs/pq-cagra-adc-research.md).

## Build

```bash
export CUDA_HOME=/usr/local/cuda-12.9
cmake -S vendor/cagra-pq-adc -B vendor/cagra-pq-adc/build \
  -DPython_EXECUTABLE=$PWD/.venv/bin/python
cmake --build vendor/cagra-pq-adc/build -j
```

The extension is copied next to `python/pq_cagra_adc/`.

## Smoke

```bash
.venv/bin/python bin/smoke_pq_cagra_adc.py \
  --structures /path/to/rouskin_sample_6k.tsv \
  --queries /path/to/queries_rouskin_structures.tsv \
  --embeddings /path/to/embeddings.float32.npy \
  --sample-frac 0.20
```
