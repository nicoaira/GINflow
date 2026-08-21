#!/usr/bin/env python3
"""Smoke tests for packed node-PQ windows and the custom HNSW distance."""
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
from benchmark_hnswlib import Record  # noqa: E402
from benchmark_pq_hnswlib import (  # noqa: E402
    OriginalNodeWindowStore,
    pack_pq_codes,
    pack_window_codes,
    rerank_candidates_allow_missing,
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

    def test_window_packing_keeps_node_then_subcode_order(self) -> None:
        nodes = np.arange(10, dtype=np.uint8).reshape(5, 2)
        record = Record("record", 5, 0, 0, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "windows.codes.bin"
            pack_window_codes(nodes, [record], 3, 1, 2, 8, path)
            actual = np.fromfile(path, dtype=np.uint8).reshape(3, 6)
        expected = np.asarray([nodes[start : start + 3].reshape(-1) for start in range(3)], dtype=np.uint8)
        np.testing.assert_array_equal(actual, expected)

    def test_rerank_masks_short_cpp_result_rows(self) -> None:
        database = np.eye(4, dtype=np.float32)
        queries = np.eye(2, 4, dtype=np.float32)
        labels = np.asarray([[0, 1, np.iinfo(np.uint64).max], [1, np.iinfo(np.uint64).max, 0]], dtype=np.uint64)
        final_labels, scores, _ = rerank_candidates_allow_missing(database, queries, labels, 2, 2)
        np.testing.assert_array_equal(final_labels, np.asarray([[0, 1], [1, 0]], dtype=np.int64))
        np.testing.assert_allclose(scores, np.asarray([[1.0, 0.0], [1.0, 0.0]], dtype=np.float32))

    def test_original_node_window_store_reconstructs_normalized_windows(self) -> None:
        embeddings = np.asarray(
            [[1, 0], [0, 2], [3, 0], [0, 4], [5, 0]], dtype=np.float16
        )
        records = [Record("a", 3, 0, 0, 2), Record("b", 2, 3, 2, 1)]
        store = OriginalNodeWindowStore(embeddings, records, 2, 1)
        expected = np.asarray(
            [
                [1, 0, 0, 2],
                [0, 2, 3, 0],
                [0, 4, 5, 0],
            ],
            dtype=np.float32,
        )
        expected /= np.linalg.norm(expected, axis=1, keepdims=True)
        np.testing.assert_allclose(store[[0, 1, 2]], expected, rtol=1e-6, atol=1e-6)

    def test_cpp_distance_sums_positional_subquantizer_tables(self) -> None:
        rng = np.random.default_rng(13)
        count, window_size, pq_m, nbits = 48, 3, 4, 4
        ksub = 1 << nbits
        codes = rng.integers(0, ksub, size=(count, window_size, pq_m), dtype=np.uint8)
        packed = pack_pq_codes(codes.reshape(-1, pq_m), pq_m, nbits).reshape(count, -1)
        basis = rng.normal(size=(pq_m, ksub, 6)).astype(np.float32)
        similarity = np.einsum("mkd,mld->mkl", basis, basis, optimize=True).astype(np.float32)
        codes_path = self.work / "codes.bin"
        similarity_path = self.work / "similarity.bin"
        index_path = self.work / "index.bin"
        labels_path = self.work / "labels.bin"
        distances_path = self.work / "distances.bin"
        packed.tofile(codes_path)
        similarity.tofile(similarity_path)

        subprocess.run(
            [
                str(self.executable), "build",
                "--codes", str(codes_path), "--similarity", str(similarity_path),
                "--count", str(count), "--window-size", str(window_size),
                "--pq-m", str(pq_m), "--nbits", str(nbits), "--index", str(index_path),
                "--m", "8", "--ef-construction", "64", "--ef-search", "64", "--threads", "1",
            ],
            check=True,
        )
        subprocess.run(
            [
                str(self.executable), "search",
                "--codes", str(codes_path), "--similarity", str(similarity_path),
                "--query-count", str(count), "--window-size", str(window_size),
                "--pq-m", str(pq_m), "--nbits", str(nbits), "--index", str(index_path),
                "--k", "5", "--ef-search", "64", "--threads", "1",
                "--labels-out", str(labels_path), "--distances-out", str(distances_path),
            ],
            check=True,
        )

        labels = np.fromfile(labels_path, dtype="<u8").reshape(count, 5)
        distances = np.fromfile(distances_path, dtype="<f4").reshape(count, 5)
        self.assertTrue(np.all(labels < count))
        for row in range(count):
            expected = np.asarray(
                [
                    similarity[
                        np.arange(pq_m)[None, :],
                        codes[row],
                        codes[target],
                    ].sum()
                    for target in labels[row]
                ],
                dtype=np.float32,
            )
            np.testing.assert_allclose(-distances[row], expected, rtol=1e-5, atol=1e-5)
        self.assertGreaterEqual(
            float(np.mean(np.any(labels == np.arange(count)[:, None], axis=1))), 0.9
        )


if __name__ == "__main__":
    unittest.main()
