This directory holds the NVIDIA cuVS CAGRA sources that this index is
adapted from.  The full cuVS tree cannot be compiled in GINflow: it depends
on RAFT, RMM, and a generated JIT dataset-descriptor table, and its Python
API still has no hook for a custom device distance.

`hashmap.hpp` is copied from cuVS (`cpp/src/neighbors/detail/cagra/hashmap.hpp`).
The runtime kernel in `src/cagra_pq_adc.cu` uses the same linear-probing
visited map, CAGRA-style random seed pickup, parent expansion, and
reverse-edge + prune graph construction, with `compute_distance` replaced by
a shared-memory ADC lookup table over node-level PQ codes.
