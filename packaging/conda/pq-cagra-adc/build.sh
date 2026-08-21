#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${CUDA_HOME:-}" ]]; then
  if [[ -x "${PREFIX}/bin/nvcc" ]]; then
    export CUDA_HOME="${PREFIX}"
  elif [[ -x /usr/local/cuda-11.8/bin/nvcc ]]; then
    export CUDA_HOME=/usr/local/cuda-11.8
  elif [[ -x /usr/local/cuda/bin/nvcc ]]; then
    export CUDA_HOME=/usr/local/cuda
  fi
fi
if [[ -n "${CUDA_HOME:-}" ]]; then
  export PATH="${CUDA_HOME}/bin:${PATH}"
fi

cmake -S . -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPQ_CAGRA_CPU_ONLY=OFF \
  -DPython_EXECUTABLE="$PYTHON" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  ${CUDA_HOME:+-DCMAKE_CUDA_COMPILER="${CUDA_HOME}/bin/nvcc"} \
  -DCMAKE_CUDA_ARCHITECTURES="70-real;75-real;80-real;86-real"
cmake --build build --parallel "${CPU_COUNT:-2}"
cmake --install build
