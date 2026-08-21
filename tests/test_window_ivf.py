#!/usr/bin/env python3
"""Smoke tests for the standalone custom-distance WindowIVF prototype."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "bin" / "window_ivf.cpp"


@unittest.skipUnless(shutil.which("g++"), "g++ is not installed")
class TestWindowIVF(unittest.TestCase):
    def test_custom_score_and_full_probe(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ginflow-window-ivf-") as tmp:
            root = Path(tmp)
            executable = root / "window_ivf"
            subprocess.run(
                ["g++", "-O2", "-std=c++17", "-fopenmp", str(SOURCE), "-o", str(executable)],
                check=True,
            )
            rng = np.random.default_rng(21)
            count, window_size, pq_m, nbits = 36, 3, 2, 4
            ksub = 1 << nbits
            codes = rng.integers(0, ksub, size=(count, window_size, pq_m), dtype=np.uint8)
            packed_bits = ((codes.reshape(-1, pq_m)[:, :, None] >> np.arange(nbits)) & 1).reshape(
                count * window_size, pq_m * nbits
            )
            packed = np.packbits(packed_bits, axis=1, bitorder="little").reshape(count, -1)
            basis = rng.normal(size=(pq_m, ksub, 5)).astype(np.float32)
            similarity = np.einsum("mkd,mld->mkl", basis, basis, optimize=True).astype(np.float32)
            codes_path = root / "codes.bin"
            similarity_path = root / "similarity.bin"
            index_path = root / "index.bin"
            labels_path = root / "labels.bin"
            distances_path = root / "distances.bin"
            packed.tofile(codes_path)
            similarity.tofile(similarity_path)
            common = [
                "--window-size", str(window_size), "--pq-m", str(pq_m), "--nbits", str(nbits),
                "--similarity", str(similarity_path), "--threads", "1",
            ]
            subprocess.run(
                [str(executable), "build", "--codes", str(codes_path), "--count", str(count),
                 "--index", str(index_path), "--nlist", "6", "--nprobe", "6", *common],
                check=True,
            )
            subprocess.run(
                [str(executable), "search", "--codes", str(codes_path), "--query-count", str(count),
                 "--index", str(index_path), "--k", "5", "--nprobe", "6",
                 "--labels-out", str(labels_path), "--distances-out", str(distances_path), *common],
                check=True,
            )
            labels = np.fromfile(labels_path, dtype="<u8").reshape(count, 5)
            distances = np.fromfile(distances_path, dtype="<f4").reshape(count, 5)
            self.assertTrue(np.all(labels < count))
            for row in range(count):
                expected = np.asarray(
                    [similarity[np.arange(pq_m)[None, :], codes[row], codes[target]].sum()
                     for target in labels[row]],
                    dtype=np.float32,
                )
                np.testing.assert_allclose(-distances[row], expected, rtol=1e-5, atol=1e-5)
            self.assertGreaterEqual(float(np.mean(np.any(labels == np.arange(count)[:, None], axis=1))), 0.95)


if __name__ == "__main__":
    unittest.main()
