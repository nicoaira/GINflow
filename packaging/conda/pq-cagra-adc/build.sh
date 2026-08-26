#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${BUILD_PREFIX:-}" && -x "${BUILD_PREFIX}/bin/nvcc" ]]; then
  export CUDA_HOME="${BUILD_PREFIX}"
elif [[ -n "${CUDA_HOME:-}" && -x "${CUDA_HOME}/bin/nvcc" ]]; then
  export CUDA_HOME
else
  echo "CUDA 11.8 nvcc was not found in the conda build environment" >&2
  exit 1
fi
export CUDACXX="${CUDA_HOME}/bin/nvcc"
export PATH="${CUDA_HOME}/bin:${PATH}"

cmake --fresh -S . -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPQ_CAGRA_CPU_ONLY=OFF \
  -DPython_EXECUTABLE="$PYTHON" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_CUDA_COMPILER="${CUDACXX}" \
  -DCMAKE_CUDA_HOST_COMPILER="${CXX}" \
  -DCUDAToolkit_ROOT="$BUILD_PREFIX" \
  -DCMAKE_CUDA_ARCHITECTURES="70-real;75-real;80-real;86-real;89-real;89-virtual"
cmake --build build --parallel "${CPU_COUNT:-2}"
cmake --install build
