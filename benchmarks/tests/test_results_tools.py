#!/usr/bin/env python3
"""Focused unit tests for benchmark validation and static reporting."""
from __future__ import annotations

import copy
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks"))

from plot_results import main as plot_main, render_scope_svg  # noqa: E402
from results_common import LocatedRecord, read_records  # noqa: E402
from validate_results import main as validate_main, validate_records  # noqa: E402


FIXTURE = Path(__file__).resolve().parent / "fixtures" / "results_valid.jsonl"


class TestBenchmarkResultTools(unittest.TestCase):
    def _records(self) -> list[LocatedRecord]:
        records, issues, ignored = read_records([FIXTURE])
        self.assertFalse(issues)
        self.assertFalse(ignored)
        return records

    def test_valid_fixture_aggregates_repeats_and_preserves_skip(self) -> None:
        records = self._records()
        # Exercise the exact schema identifier emitted by benchmark_utils.py.
        for record in records:
            record.payload["schema_version"] = "ginflow-benchmark-v1"
        report = validate_records(records)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["valid_successful_records"], 6)
        self.assertEqual(report["successful_configurations"], 3)
        self.assertEqual(len(report["unavailable"]), 1)
        faiss = next(row for row in report["aggregates"] if row["backend"] == "faiss")
        self.assertEqual(faiss["qps_median"], 105.0)
        self.assertEqual(faiss["qps_min"], 100.0)
        self.assertEqual(faiss["qps_max"], 110.0)

    def test_invalid_recall_and_missing_provenance_are_rejected(self) -> None:
        records = self._records()
        payload = copy.deepcopy(records[0].payload)
        payload["recall_at_100"] = 1.2
        payload["provenance"].pop("query_selection_id")
        report = validate_records([LocatedRecord(payload, "bad.json")])
        self.assertFalse(report["valid"])
        message = "\n".join(issue["message"] for issue in report["errors"])
        self.assertIn("recall_at_100", message)
        self.assertIn("query_selection_id", message)

    def test_incomparable_ground_truth_is_rejected(self) -> None:
        records = self._records()
        mutated = copy.deepcopy(records[2].payload)
        mutated["ground_truth_ids_sha256"] = "c" * 64
        report = validate_records(records[:2] + [LocatedRecord(mutated, "different-ground-truth.json")])
        self.assertFalse(report["valid"])
        self.assertTrue(any(issue["code"] == "incomparable_scope" for issue in report["errors"]))

    def test_cli_generates_svg_and_markdown_without_plotting_threshold_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "report"
            exit_code = plot_main([str(FIXTURE), "--output-dir", str(output)])
            self.assertEqual(exit_code, 0)
            svgs = sorted(output.glob("*.svg"))
            self.assertEqual(len(svgs), 2)
            all_svg = next(path.read_text(encoding="utf-8") for path in svgs if "_all_" in path.name)
            focused_svg = next(path.read_text(encoding="utf-8") for path in svgs if "_gt_" in path.name)
            self.assertIn("All validated points are shown.", all_svg)
            self.assertIn("ivfpq-nprobe-32", all_svg)
            self.assertIn("ah-leaves-512", all_svg)
            self.assertIn("qbg-edge-64", all_svg)
            self.assertIn("ivfpq-nprobe-32", focused_svg)
            self.assertIn("ah-leaves-512", focused_svg)
            self.assertNotIn("qbg-edge-64", focused_svg)
            summary = (output / "benchmark_results.md").read_text(encoding="utf-8")
            self.assertIn("Below or equal to the focused recall threshold", summary)
            self.assertIn("All measured points", summary)
            self.assertIn("GPU backend deliberately unavailable", summary)
            validation = json.loads((output / "validation.json").read_text(encoding="utf-8"))
            self.assertTrue(validation["valid"])

    def test_custom_threshold_controls_svg_scale_and_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "report"
            exit_code = plot_main(
                [str(FIXTURE), "--output-dir", str(output), "--recall-threshold", "0.90"]
            )
            self.assertEqual(exit_code, 0)
            svgs = sorted(output.glob("*.svg"))
            self.assertEqual(len(svgs), 2)
            svg = next(path.read_text(encoding="utf-8") for path in svgs if "_gt_" in path.name)
            self.assertIn("Only recall@100 &gt; 0.900 is shown.", svg)
            self.assertNotIn("Only recall@100 &gt; 0.80 is shown.", svg)
            self.assertNotRegex(svg.lower(), r"\b(?:nan|inf)\b")
            self.assertIn("ivfpq-nprobe-32", svg)
            self.assertNotIn("ah-leaves-512", svg)

    def test_svg_is_deterministic_for_one_eligible_point(self) -> None:
        report = validate_records(self._records())
        scope = [row for row in report["aggregates"] if row["recall_median"] > 0.94]
        self.assertEqual(len(scope), 1)
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "first.svg"
            second = Path(tmp) / "second.svg"
            render_scope_svg(scope, first, "synthetic source", recall_threshold=0.94)
            render_scope_svg(scope, second, "synthetic source", recall_threshold=0.94)
            svg = first.read_text(encoding="utf-8")
            self.assertEqual(svg, second.read_text(encoding="utf-8"))
            self.assertIn("Only recall@100 &gt; 0.94 is shown.", svg)
            self.assertNotIn("generated ", svg)

    def test_no_svg_is_created_when_every_point_is_below_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "report"
            exit_code = plot_main(
                [str(FIXTURE), "--output-dir", str(output), "--recall-threshold", "0.99"]
            )
            self.assertEqual(exit_code, 0)
            svgs = list(output.glob("*.svg"))
            self.assertEqual(len(svgs), 1)
            self.assertIn("_all_", svgs[0].name)
            self.assertIn("All validated points are shown.", svgs[0].read_text(encoding="utf-8"))
            summary = (output / "benchmark_results.md").read_text(encoding="utf-8")
            self.assertIn("Below or equal to the focused recall threshold", summary)

    def test_multiple_slug_colliding_scopes_receive_distinct_svg_names(self) -> None:
        records = self._records()
        duplicated = [copy.deepcopy(record.payload) for record in records]
        for payload in duplicated:
            payload["dataset_id"] = "synthetic-6k"
        with tempfile.TemporaryDirectory() as tmp:
            inputs = Path(tmp) / "results.jsonl"
            inputs.write_text(
                "".join(json.dumps(record.payload) + "\n" for record in records)
                + "".join(json.dumps(payload) + "\n" for payload in duplicated),
                encoding="utf-8",
            )
            output = Path(tmp) / "report"
            self.assertEqual(plot_main([str(inputs), "--output-dir", str(output)]), 0)
            svgs = sorted(output.glob("*.svg"))
            self.assertEqual(len(svgs), 4)
            self.assertEqual(len({path.name for path in svgs}), 4)

    def test_markdown_metadata_is_escaped_without_breaking_table_rows(self) -> None:
        records = self._records()
        for record in records[:2]:
            record.payload["parameter_label"] = "ivf|*unsafe*\nnext"
        records[-1].payload["error"] = "not available | verify\nGPU memory"
        with tempfile.TemporaryDirectory() as tmp:
            inputs = Path(tmp) / "results.jsonl"
            inputs.write_text("".join(json.dumps(record.payload) + "\n" for record in records), encoding="utf-8")
            output = Path(tmp) / "report"
            self.assertEqual(plot_main([str(inputs), "--output-dir", str(output)]), 0)
            summary = (output / "benchmark_results.md").read_text(encoding="utf-8")
            self.assertIn("ivf\\|\\*unsafe\\*<br>next", summary)
            self.assertIn("not available \\| verify<br>GPU memory", summary)

    def test_measurement_contract_and_scope_comparability_are_enforced(self) -> None:
        records = self._records()
        malformed = copy.deepcopy(records[0].payload)
        malformed["measurement"]["timed_queries"] = 128
        malformed["latency_ms"]["p50"] = 12.0
        malformed["latency_ms"]["p95"] = 11.0
        report = validate_records([LocatedRecord(malformed, "bad-measurement.json")])
        self.assertFalse(report["valid"])
        codes = {issue["code"] for issue in report["errors"]}
        self.assertIn("invalid_measurement", codes)
        self.assertIn("invalid_latency_percentiles", codes)

        mismatch = copy.deepcopy(records[0].payload)
        mismatch["measurement"]["query_batch_size"] = 128
        mismatch["measurement"]["timed_batch_count"] = 2
        report = validate_records(records + [LocatedRecord(mismatch, "different-batch.json")], min_repeats=1)
        self.assertFalse(report["valid"])
        self.assertTrue(any(issue["code"] == "incomparable_scope" for issue in report["errors"]))

    def test_recursive_scan_deduplicates_overlapping_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            copied = Path(tmp) / "nested" / "results.jsonl"
            copied.parent.mkdir(parents=True)
            copied.write_text(FIXTURE.read_text(encoding="utf-8"), encoding="utf-8")
            records, issues, ignored = read_records([Path(tmp), copied.resolve()])
            self.assertFalse(issues)
            self.assertFalse(ignored)
            self.assertEqual(len(records), 7)

    def test_generated_validation_report_is_ignored_on_a_repeat_scan(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "report"
            self.assertEqual(plot_main([str(FIXTURE), "--output-dir", str(output)]), 0)
            records, issues, ignored = read_records([output, FIXTURE])
            self.assertFalse(issues)
            self.assertEqual(len(records), 7)
            self.assertTrue(any(path.endswith("validation.json") for path in ignored))

    def test_validator_cli_returns_nonzero_for_bad_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.json"
            payload = copy.deepcopy(self._records()[0].payload)
            payload["qps"] = 0
            bad.write_text(json.dumps(payload), encoding="utf-8")
            output = Path(tmp) / "validation.json"
            self.assertEqual(validate_main([str(bad), "--output", str(output)]), 2)
            report = json.loads(output.read_text(encoding="utf-8"))
            self.assertFalse(report["valid"])


if __name__ == "__main__":
    raise SystemExit(unittest.main())
