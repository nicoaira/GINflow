# Vendored hnswlib headers

This directory contains the upstream hnswlib v0.8.0 header-only implementation
plus `compact_hnswlib.cpp`, the small GINflow driver that registers a custom
distance over uint16 centroid-code windows.

Upstream project: <https://github.com/nmslib/hnswlib>

The upstream headers and license are retained here so Nextflow tasks can build
the custom C++ space without a network checkout. `compact_hnswlib.cpp` is
GINflow code and is intentionally separate from the upstream headers. The
third-party license is in `LICENSE`.

To rebuild the standalone driver from the repository root:

```bash
g++ -O3 -std=c++11 -fopenmp \
  -I vendor/hnswlib-0.8.0 \
  vendor/hnswlib-0.8.0/compact_hnswlib.cpp \
  -o compact_hnswlib
```
