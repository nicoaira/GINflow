#!/usr/bin/env python3
"""Unit tests for selecting whole query-target pairs for plotting."""
from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

# The selection helper has no SW dependency, but plot_sw imports the package at
# module load time. Keep this unit test runnable in the lightweight test env.
if "ginfinity_sw" not in sys.modules:
    fake_sw = types.ModuleType("ginfinity_sw")
    fake_sw.ScoringParameters = object
    fake_sw.align = object()
    fake_sw.similarity_matrix = object()
    fake_sw.transform_scores = object()
    sys.modules["ginfinity_sw"] = fake_sw

from plot_r4rna import selected_pair_rows as r4_selected_pair_rows  # noqa: E402
from plot_rnartistcore import selected_pair_rows as rnartist_selected_pair_rows  # noqa: E402
from plot_sw import selected_pair_rows as sw_selected_pair_rows  # noqa: E402


class TestPlotPairSelection(unittest.TestCase):
    def test_limit_is_unique_pairs_and_retains_every_hsp(self) -> None:
        rows = [
            {"query_id": "q1", "target_id": "t1", "cluster_id": "a"},
            {"query_id": "q1", "target_id": "t1", "cluster_id": "b"},
            {"query_id": "q1", "target_id": "t1", "cluster_id": "c"},
            {"query_id": "q2", "target_id": "t2", "cluster_id": "d"},
            {"query_id": "q3", "target_id": "t3", "cluster_id": "e"},
        ]
        expected = [(0, "a"), (1, "b"), (2, "c"), (3, "d")]

        for selector in (
            rnartist_selected_pair_rows,
            r4_selected_pair_rows,
            sw_selected_pair_rows,
        ):
            selected = selector(rows, 2)
            self.assertEqual(
                [(index, row["cluster_id"]) for index, row in selected],
                expected,
            )


if __name__ == "__main__":
    raise SystemExit(unittest.main())
