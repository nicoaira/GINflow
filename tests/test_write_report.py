#!/usr/bin/env python3
"""Unit tests for pair-level report rendering."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from write_report import build_hits, plot_panel  # noqa: E402


class TestWriteReport(unittest.TestCase):
    def test_renders_every_hsp_plot_for_a_collapsed_pair(self) -> None:
        rows = [{
            "cluster_id": "0",
            "cluster_ids": "0,1,2",
            "query_id": "query",
            "target_id": "target",
            "score": "16",
            "total_score": "16",
            "max_score": "10",
            "alignment_count": "3",
            "evalue": "1e-6",
            "evalue_pair": "1e-7",
            "query_start": "0",
            "query_end": "20",
            "target_start": "0",
            "target_end": "20",
            "query_length": "20",
            "target_length": "20",
            "match_count": "10",
            "aligned_columns": "20",
        }]
        plots = {
            "0_query__target_alignment.svg": '<svg id="same"></svg>',
            "1_query__target_alignment.svg": '<svg id="same"></svg>',
            "2_query__target_alignment.svg": '<svg id="same"></svg>',
        }

        hit = build_hits(rows, {}, {}, plots)[0]
        self.assertEqual(len(hit["plot_r4"]), 3)
        panel = plot_panel(hit)
        self.assertEqual(panel.count('<svg id="same"></svg>'), 3)
        self.assertIn("aligned structures · HSP 1", panel)
        self.assertIn("aligned structures · HSP 2", panel)
        self.assertIn("aligned structures · HSP 3", panel)

if __name__ == "__main__":
    raise SystemExit(unittest.main())
