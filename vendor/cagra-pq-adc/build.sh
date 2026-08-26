#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
if [[ -n "${BUILD_PREFIX:-}" ]]; then
  CUDA_COMPILER="${BUILD_PREFIX}/bin/nvcc"
  if [[ ! -x "${CUDA_COMPILER}" ]]; then
    echo "Conda-build CUDA compiler not found: ${CUDA_COMPILER}" >&2
    exit 1
  fi
else
  CUDA_HOME="${CUDA_HOME:-/usr/local/cuda-12.9}"
  CUDA_COMPILER="${CUDA_HOME}/bin/nvcc"
fi
PYTHON="${PYTHON:-${ROOT}/../../.venv/bin/python}"
if [[ ! -x "${PYTHON}" ]]; then
  PYTHON="$(command -v python3)"
fi
CUDA_HOME="$(cd "$(dirname "${CUDA_COMPILER}")/.." && pwd)"
export PATH="${CUDA_HOME}/bin:${PATH}"
CUDA_TOOLKIT_ROOT="${BUILD_PREFIX:-${PREFIX:-${CUDA_HOME}}}"
cmake --fresh -S "${ROOT}" -B "${ROOT}/build" \
  -DPython_EXECUTABLE="${PYTHON}" \
  -DCMAKE_CUDA_COMPILER="${CUDA_COMPILER}" \
  -DCMAKE_CUDA_HOST_COMPILER="${CXX:-}" \
  -DCUDAToolkit_ROOT="${CUDA_TOOLKIT_ROOT}" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "${ROOT}/build" -j"${CMAKE_BUILD_PARALLEL_LEVEL:-4}"
