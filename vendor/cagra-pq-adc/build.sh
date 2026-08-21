#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
CUDA_HOME="${CUDA_HOME:-/usr/local/cuda-12.9}"
PYTHON="${PYTHON:-${ROOT}/../../.venv/bin/python}"
if [[ ! -x "${PYTHON}" ]]; then
  PYTHON="$(command -v python3)"
fi
export PATH="${CUDA_HOME}/bin:${PATH}"
cmake -S "${ROOT}" -B "${ROOT}/build" \
  -DPython_EXECUTABLE="${PYTHON}" \
  -DCMAKE_CUDA_COMPILER="${CUDA_HOME}/bin/nvcc" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "${ROOT}/build" -j"${CMAKE_BUILD_PARALLEL_LEVEL:-4}"
