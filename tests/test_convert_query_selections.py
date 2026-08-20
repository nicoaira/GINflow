#!/usr/bin/env python3
"""Tests for converting window selections into sliced structures inputs."""
from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from convert_query_selections import convert  # noqa: E402


class TestConvertQuerySelections(unittest.TestCase):
    def test_writes_pipeline_slices_and_provenance(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ginflow-query-convert-") as directory:
            work = Path(directory)
            structures = work / "structures.tsv"
            selections = work / "selections.tsv"
            output = work / "queries.tsv"
            structures.write_text(
                "transcript_id\tsequence\tsecondary_structure\n"
                "rna1\tACGTACGTACGT\t............\n"
                "rna2\tAAAAACCCCCCC\t............\n"
            )
            selections.write_text(
                "query_ordinal\ttranscript_id\twindow_offset\tselection_cycle\n"
                "0\trna1\t2\t0\n"
                "1\trna2\t1\t3\n"
            )

            result = convert(structures, selections, output, 4)

            self.assertEqual(result["source_records"], 2)
            self.assertEqual(result["query_rows"], 2)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["sequence"], "ACGUACGUACGU")
            self.assertEqual(rows[0]["start"], "2")
            self.assertEqual(rows[0]["end"], "6")
            self.assertEqual(rows[0]["query_ordinal"], "0")
            self.assertEqual(rows[1]["start"], "1")
            self.assertEqual(rows[1]["end"], "5")

    def test_rejects_out_of_range_window(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ginflow-query-convert-") as directory:
            work = Path(directory)
            structures = work / "structures.tsv"
            selections = work / "selections.tsv"
            output = work / "queries.tsv"
            structures.write_text(
                "transcript_id\tsequence\tsecondary_structure\n"
                "rna1\tAAAA\t....\n"
            )
            selections.write_text("transcript_id\twindow_offset\nrna1\t2\n")

            with self.assertRaisesRegex(ValueError, "outside"):
                convert(structures, selections, output, 4)

    def test_can_write_deduplicated_full_molecules(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ginflow-query-convert-") as directory:
            work = Path(directory)
            structures = work / "structures.tsv"
            selections = work / "selections.tsv"
            output = work / "queries.tsv"
            structures.write_text(
                "transcript_id\tsequence\tsecondary_structure\n"
                "rna1\tACGTACGTACGT\t............\n"
                "rna2\tAAAAACCCCCCC\t............\n"
            )
            selections.write_text(
                "transcript_id\twindow_offset\n"
                "rna1\t2\n"
                "rna1\t5\n"
                "rna2\t1\n"
            )

            result = convert(structures, selections, output, 4, full_molecules=True)

            self.assertEqual(result["query_rows"], 3)
            self.assertEqual(result["output_records"], 2)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["transcript_id"] for row in rows], ["rna1", "rna2"])
            self.assertNotIn("start", rows[0])
            self.assertEqual(rows[0]["sequence"], "ACGUACGUACGU")


if __name__ == "__main__":
    unittest.main()
