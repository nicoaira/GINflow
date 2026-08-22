#!/usr/bin/env python3
"""Unit tests for BLAST-style query-target aggregation."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from merge_alignments import aggregate_pair, aggregate_rows  # noqa: E402


def hsp(cluster_id: str, score: float, query_start: int, target_start: int) -> dict[str, str]:
    return {
        "cluster_id": cluster_id,
        "query_id": "query",
        "target_id": "target",
        "score": str(score),
        "query_start": str(query_start),
        "query_end": str(query_start + 10),
        "target_start": str(target_start),
        "target_end": str(target_start + 10),
        "query_length": "50",
        "target_length": "80",
        "match_count": "8",
        "aligned_columns": "10",
        "ungapped_columns": "10",
        "seed_count": "2",
    }


class TestMergeAlignments(unittest.TestCase):
    def test_pair_scores_and_evalues_use_total_score(self) -> None:
        evd = {"lambda": 2.0, "K": 0.5, "database_residues": 100}
        result = aggregate_pair(
            [
                hsp("a", 10, 2, 3),
                hsp("b", 6, 20, 30),
                hsp("c", 4, 38, 55),
            ],
            evd,
        )

        self.assertEqual(result["score"], "20")
        self.assertEqual(result["total_score"], "20")
        self.assertEqual(result["max_score"], "10")
        self.assertEqual(result["alignment_count"], "3")
        self.assertEqual(result["cluster_ids"], "a,b,c")
        self.assertEqual(result["hsp_scores"], "[10.0,6.0,4.0]")
        expected_database_e = 0.5 * 50 * 100 * math.exp(-2.0 * 20.0)
        expected_pair_e = 0.5 * 50 * 80 * math.exp(-2.0 * 20.0)
        self.assertAlmostEqual(float(result["evalue"]), expected_database_e)
        self.assertAlmostEqual(float(result["evalue_pair"]), expected_pair_e)

    def test_rows_group_by_query_and_target(self) -> None:
        rows = [hsp("b", 6, 20, 30), hsp("a", 10, 2, 3)]
        rows.append({**hsp("c", 4, 1, 1), "target_id": "other"})
        merged = aggregate_rows(rows, {"lambda": 1.0, "K": 1.0, "database_residues": 100})

        self.assertEqual(len(merged), 2)
        query_target = next(row for row in merged if row["target_id"] == "target")
        self.assertEqual(query_target["cluster_ids"], "a,b")
        self.assertEqual(query_target["hsp_scores"], "[10.0,6.0]")


if __name__ == "__main__":
    raise SystemExit(unittest.main())
