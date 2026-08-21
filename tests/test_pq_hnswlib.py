#!/usr/bin/env python3
"""Smoke tests for packed node-PQ windows and ADC HNSW search."""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "vendor" / "hnswlib-0.8.0" / "pq_hnswlib.cpp"
INCLUDE = ROOT / "vendor" / "hnswlib-0.8.0"
sys.path.insert(0, str(ROOT / "bin"))
from node_quantization import (  # noqa: E402
    pack_pq_codes,
    sdc_lookup_table,
    train_node_pq,
    unpack_pq_codes,
)


@unittest.skipUnless(shutil.which("g++"), "g++ is not installed")
class TestPqHnswlib(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp = tempfile.TemporaryDirectory(prefix="ginflow-pq-hnsw-")
        cls.work = Path(cls.tmp.name)
        cls.executable = cls.work / "pq_hnswlib"
        subprocess.run(
            [
                "g++", "-O2", "-std=c++11", "-fopenmp", "-I", str(INCLUDE),
                str(SOURCE), "-o", str(cls.executable),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def test_packed_codes_round_trip(self) -> None:
        rng = np.random.default_rng(7)
        for pq_m, nbits in ((4, 4), (8, 4), (16, 8)):
            codes = rng.integers(0, 1 << nbits, size=(13, pq_m), dtype=np.uint8)
            packed = pack_pq_codes(codes, pq_m, nbits)
            np.testing.assert_array_equal(unpack_pq_codes(packed, pq_m, nbits), codes)

    def test_adc_search_returns_self_hits(self) -> None:
        rng = np.random.default_rng(11)
        window_size = 3
        dim = 16
        pq_m = 4
        nbits = 4
        n_windows = 32
        nodes = rng.standard_normal((n_windows * window_size, dim), dtype=np.float32)
        nodes /= np.linalg.norm(nodes, axis=1, keepdims=True).clip(min=1e-12)
        codebook, node_codes = train_node_pq(nodes, pq_m, nbits, sample_size=len(nodes), niter=6, seed=1)
        windows = node_codes.reshape(n_windows, window_size, pq_m)
        packed = pack_pq_codes(windows.reshape(-1, pq_m), pq_m, nbits).reshape(
            n_windows, window_size * pack_pq_codes(windows[0], pq_m, nbits).shape[1]
        )
        work = self.work / "adc"
        work.mkdir(exist_ok=True)
        codes_path = work / "codes.bin"
        sim_path = work / "sdc.bin"
        index_path = work / "index.bin"
        packed.tofile(codes_path)
        sdc_lookup_table(codebook).tofile(sim_path)
        subprocess.run(
            [
                str(self.executable), "build",
                "--codes", str(codes_path),
                "--similarity", str(sim_path),
                "--index", str(index_path),
                "--count", str(n_windows),
                "--window-size", str(window_size),
                "--pq-m", str(pq_m),
                "--nbits", str(nbits),
                "--m", "8",
                "--ef-construction", "40",
                "--ef-search", "40",
                "--threads", "1",
            ],
            check=True,
        )
        queries = nodes.reshape(n_windows, window_size, dim)[:4]
        query_path = work / "queries.bin"
        codebook_path = work / "codebook.bin"
        labels_path = work / "labels.bin"
        distances_path = work / "distances.bin"
        queries.astype(np.float32).tofile(query_path)
        np.ascontiguousarray(codebook, dtype=np.float32).tofile(codebook_path)
        subprocess.run(
            [
                str(self.executable), "search",
                "--queries", str(query_path),
                "--codebook", str(codebook_path),
                "--index", str(index_path),
                "--query-count", "4",
                "--window-size", str(window_size),
                "--pq-m", str(pq_m),
                "--nbits", str(nbits),
                "--dim", str(dim),
                "--k", "1",
                "--ef-search", "40",
                "--labels-out", str(labels_path),
                "--distances-out", str(distances_path),
            ],
            check=True,
        )
        labels = np.fromfile(labels_path, dtype=np.uint64)
        self.assertEqual(labels.shape[0], 4)
        self.assertTrue(np.all(labels < n_windows))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
