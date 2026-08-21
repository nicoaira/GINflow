#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPQ_CAGRA_CPU_ONLY=ON \
  -DPython_EXECUTABLE="$PYTHON" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build build --parallel "${CPU_COUNT:-2}"
cmake --install build
