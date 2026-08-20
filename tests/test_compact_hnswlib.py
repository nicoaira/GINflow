#!/usr/bin/env python3
"""Smoke tests for the vendored compact C++ hnswlib space."""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "vendor" / "hnswlib-0.8.0" / "compact_hnswlib.cpp"
INCLUDE = ROOT / "vendor" / "hnswlib-0.8.0"
sys.path.insert(0, str(ROOT / "bin"))
from benchmark_compact_hnswlib import threshold_stats  # noqa: E402
from search_compact_hnswlib import rerank_candidates  # noqa: E402


@unittest.skipUnless(shutil.which("g++"), "g++ is not installed")
class TestCompactHnswlib(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp = tempfile.TemporaryDirectory(prefix="ginflow-compact-hnsw-")
        cls.work = Path(cls.tmp.name)
        cls.executable = cls.work / "compact_hnswlib"
        subprocess.run(
            [
                "g++",
                "-O2",
                "-std=c++11",
                "-fopenmp",
                "-I",
                str(INCLUDE),
                str(SOURCE),
                "-o",
                str(cls.executable),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def test_index_uses_codes_and_returns_registered_scores(self) -> None:
        rng = np.random.default_rng(19)
        n_elements, window_size, n_centroids = 48, 4, 7
        codes = rng.integers(0, n_centroids, size=(n_elements, window_size), dtype=np.uint16)
        basis = rng.normal(size=(n_centroids, 6)).astype(np.float32)
        basis /= np.linalg.norm(basis, axis=1, keepdims=True)
        similarity = np.ascontiguousarray(basis @ basis.T, dtype=np.float32)
        codes_path = self.work / "codes.bin"
        similarity_path = self.work / "similarity.bin"
        index_path = self.work / "index.bin"
        labels_path = self.work / "labels.bin"
        distances_path = self.work / "distances.bin"
        codes.tofile(codes_path)
        similarity.tofile(similarity_path)

        subprocess.run(
            [
                str(self.executable),
                "build",
                "--codes", str(codes_path),
                "--similarity", str(similarity_path),
                "--count", str(n_elements),
                "--window-size", str(window_size),
                "--n-centroids", str(n_centroids),
                "--index", str(index_path),
                "--m", "8",
                "--ef-construction", "64",
                "--ef-search", "64",
                "--threads", "1",
            ],
            check=True,
        )
        subprocess.run(
            [
                str(self.executable),
                "search",
                "--codes", str(codes_path),
                "--similarity", str(similarity_path),
                "--query-count", str(n_elements),
                "--window-size", str(window_size),
                "--n-centroids", str(n_centroids),
                "--index", str(index_path),
                "--k", "5",
                "--ef-search", "64",
                "--threads", "1",
                "--labels-out", str(labels_path),
                "--distances-out", str(distances_path),
            ],
            check=True,
        )

        labels = np.fromfile(labels_path, dtype="<u8").reshape(n_elements, 5)
        distances = np.fromfile(distances_path, dtype="<f4").reshape(n_elements, 5)
        self.assertTrue(np.all(labels < n_elements))
        for row in range(n_elements):
            explicit = np.array(
                [similarity[codes[row], codes[target]].sum() for target in labels[row]],
                dtype=np.float32,
            )
            np.testing.assert_allclose(-distances[row], explicit, rtol=1e-6, atol=1e-6)
        self.assertGreaterEqual(
            float(np.mean(np.any(labels == np.arange(n_elements)[:, None], axis=1))),
            0.9,
        )
        self.assertLess(index_path.stat().st_size, n_elements * window_size * 128 * 4)

    def test_rerank_uses_original_windows_after_code_candidates(self) -> None:
        labels = np.asarray([[1, np.iinfo(np.uint64).max, 0]], dtype=np.uint64)
        query_mapping = [("query", 0, 2)]
        targets = [("target0", 0, 2), ("target1", 0, 2), ("target2", 0, 2)]
        query_embeddings = {"query": np.asarray([[1, 0], [0, 1]], dtype=np.float16)}
        database_embeddings = {
            "target0": np.asarray([[1, 0], [0, 1]], dtype=np.float16),
            "target1": np.asarray([[1, 0], [1, 0]], dtype=np.float16),
            "target2": np.asarray([[0, 1], [0, 1]], dtype=np.float16),
        }

        result_labels, result_scores = rerank_candidates(
            labels,
            query_mapping,
            targets,
            query_embeddings,
            database_embeddings,
            2,
        )

        np.testing.assert_array_equal(result_labels, [[0, 1]])
        np.testing.assert_allclose(result_scores[0, 0], 1.0, rtol=1e-6, atol=1e-6)
        self.assertGreater(float(result_scores[0, 1]), 0.0)

    def test_threshold_stats_compare_seed_sets_in_score_space(self) -> None:
        exact_labels = np.asarray([[0, 1, 2]], dtype=np.int64)
        exact_scores = np.asarray([[0.95, 0.85, 0.75]], dtype=np.float32)
        approximate_labels = np.asarray([[1, 0, 3]], dtype=np.int64)
        approximate_scores = np.asarray([[0.85, 0.95, 0.80]], dtype=np.float32)

        rows = threshold_stats(
            exact_labels,
            exact_scores,
            approximate_labels,
            approximate_scores,
            [0.9, 0.8],
        )

        self.assertEqual(rows[0]["exact_hits"], 1)
        self.assertEqual(rows[0]["approximate_hits"], 1)
        self.assertEqual(rows[0]["overlap_hits"], 1)
        self.assertEqual(rows[1]["exact_hits"], 2)
        self.assertEqual(rows[1]["approximate_hits"], 3)
        self.assertAlmostEqual(float(rows[1]["micro_recall"]), 1.0)
        self.assertAlmostEqual(float(rows[1]["micro_precision"]), 2.0 / 3.0)


if __name__ == "__main__":
    unittest.main()
