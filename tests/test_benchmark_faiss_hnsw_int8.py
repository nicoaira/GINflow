#!/usr/bin/env python3
"""Unit tests for the FAISS HNSW int8 research benchmark helpers."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from benchmark_faiss_hnsw_int8 import (  # noqa: E402
    parse_args,
    quantize,
    recall_at,
    rerank_candidates,
)


class TestFaissHnswInt8Benchmark(unittest.TestCase):
    def test_quantize_uses_rounding_and_float32_input_contract(self) -> None:
        values = np.asarray([[-0.0012, 0.0, 0.0012]], dtype=np.float16)
        encoded = quantize(values, 850.0)
        self.assertEqual(encoded.dtype, np.float32)
        np.testing.assert_array_equal(encoded, [[-1.0, 0.0, 1.0]])

    def test_quantize_rejects_clipping(self) -> None:
        with self.assertRaises(ValueError):
            quantize(np.asarray([[0.2]], dtype=np.float32), 850.0)

    def test_recall_and_original_window_rerank(self) -> None:
        database = np.asarray(
            [[1.0, 0.0], [0.8, 0.6], [0.0, 1.0]], dtype=np.float16
        )
        queries = np.asarray([[0.8, 0.6]], dtype=np.float32)
        candidates = np.asarray([[2, 0, 1]], dtype=np.int64)
        labels, scores, elapsed = rerank_candidates(
            database, queries, candidates, output_k=2, batch_size=1
        )
        np.testing.assert_array_equal(labels, [[1, 0]])
        np.testing.assert_allclose(scores, [[1.0, 0.8]], rtol=2e-4, atol=2e-4)
        self.assertGreaterEqual(elapsed, 0.0)
        self.assertAlmostEqual(recall_at(np.asarray([[1, 2]]), labels, 1), 1.0)

    def test_load_index_is_exposed_for_reuse_sweeps(self) -> None:
        args = parse_args(
            [
                "--database-windows", "database.npy",
                "--queries", "queries.npy",
                "--outdir", "results",
                "--load-index", "index.faiss",
            ]
        )
        self.assertEqual(args.load_index, Path("index.faiss"))


if __name__ == "__main__":
    unittest.main()
