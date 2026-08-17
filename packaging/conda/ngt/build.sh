#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DNGT_QBG_DISABLED=OFF \
  -DNGT_AVX_DISABLED=ON \
  -DNGT_MARCH_NATIVE_DISABLED=ON \
  -DNGTQG_NO_ROTATION=ON \
  -DNGTQG_ZERO_GLOBAL=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_BUILD_RPATH="$PREFIX/lib" \
  -DCMAKE_INSTALL_RPATH="$PREFIX/lib"
cmake --build build --target install --parallel "${CPU_COUNT:-2}"

sed -i \
  -e "s#/usr/local/include#${PREFIX}/include#g" \
  -e "s#/usr/local/lib64#${PREFIX}/lib#g" \
  -e "s#/usr/local/lib#${PREFIX}/lib#g" \
  -e "s/else ''/else '-march=x86-64'/g" \
  -e "s/'-march=native'/'-march=x86-64'/g" \
  python/setup.py

cd python
python setup.py install \
  --prefix "$PREFIX" \
  --single-version-externally-managed \
  --record "$PREFIX/conda-meta/ngt-files.txt" \
  --shared-library-without-avx
