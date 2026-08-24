#!/usr/bin/env python3
"""Batch SW fallback matches serial GINFINITY-SW when the library is present."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

try:
    from ginfinity_sw import ScoringParameters, align_scores
except ImportError:  # pragma: no cover
    ScoringParameters = None  # type: ignore[misc, assignment]


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))


@unittest.skipIf(ScoringParameters is None, "ginfinity_sw is not installed")
class TestSwBatch(unittest.TestCase):
    def test_fallback_matches_serial(self) -> None:
        from sw_batch import _fallback_align_scores_many, align_score_matrices

        params = ScoringParameters(mu=0.2, gamma=5.0, gap_open=6.0, gap_extend=1.0)
        rng = np.random.default_rng(3)
        matrices = [rng.standard_normal((18, 21)) for _ in range(9)]
        serial = [align_scores(matrix, params) for matrix in matrices]
        batched = _fallback_align_scores_many(
            matrices, params, traceback=True, max_cells=1_000_000, workers=4
        )
        self.assertEqual(len(batched), len(serial))
        for left, right in zip(batched, serial):
            self.assertEqual(left.score, right.score)
            self.assertEqual(left.columns, right.columns)
        library = align_score_matrices(matrices, params, workers=4)
        for left, right in zip(library, serial):
            self.assertEqual(left.score, right.score)
            self.assertEqual(left.columns, right.columns)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
