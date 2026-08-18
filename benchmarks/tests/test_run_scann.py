"""Contract and pinned-image integration tests for the ScaNN benchmark runner."""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


BENCHMARKS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCHMARKS))

from benchmark_utils import sha256_array  # noqa: E402
from results_common import read_records  # noqa: E402
from run_scann import frontier_specs, parse_args, run  # noqa: E402
from validate_results import validate_records  # noqa: E402


def _unit_rows(n_rows: int, dimension: int, seed: int = 29) -> np.ndarray:
    generator = np.random.default_rng(seed)
    values = generator.standard_normal((n_rows, dimension), dtype=np.float32)
    values /= np.linalg.norm(values, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(values, dtype=np.float32)


def write_prepared_cache(cache_dir: Path, *, dataset: str = "fixture") -> Path:
    """Make a small cache following the backend-neutral benchmark contract."""

    root = cache_dir / dataset
    flat = root / "flat"
    queries_dir = root / "queries"
    truth_dir = root / "ground-truth"
    for directory in (flat, queries_dir, truth_dir):
        directory.mkdir(parents=True)
    vectors = _unit_rows(768, 16)
    query_vector_ids = np.arange(0, 768, 24, dtype=np.int64)[:32]
    queries = np.ascontiguousarray(vectors[query_vector_ids], dtype=np.float32)
    scores = queries @ vectors.T
    identifiers = np.argsort(-scores, axis=1, kind="stable")[:, :100].astype(np.int64)
    embedding_cache_id = "embedding-fixture"
    query_selection_id = "query-fixture"
    ground_truth_cache_id = "truth-fixture"
    np.save(flat / "vectors.npy", vectors, allow_pickle=False)
    np.save(queries_dir / "queries.npy", queries, allow_pickle=False)
    np.savez_compressed(
        truth_dir / "ground-truth.npz",
        ids=identifiers,
        scores=np.take_along_axis(scores, identifiers, axis=1),
    )
    (flat / "flatten-manifest.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset,
                "embedding_cache_id": embedding_cache_id,
                "vectors": {"path": "vectors.npy"},
            }
        ),
        encoding="utf-8",
    )
    (queries_dir / "query-selection.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset,
                "embedding_cache_id": embedding_cache_id,
                "query_selection_id": query_selection_id,
                "query_ids_sha256": sha256_array(query_vector_ids),
                "query_vectors": {"path": "queries.npy"},
            }
        ),
        encoding="utf-8",
    )
    (truth_dir / "ground-truth.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset,
                "embedding_cache_id": embedding_cache_id,
                "query_selection_id": query_selection_id,
                "ground_truth_cache_id": ground_truth_cache_id,
                "metric": "cosine",
                "k": 100,
                "database_window_count": int(vectors.shape[0]),
                "dimension": int(vectors.shape[1]),
                "query_ids_sha256": sha256_array(query_vector_ids),
                "ground_truth": {"path": "ground-truth.npz", "ids_sha256": sha256_array(identifiers)},
            }
        ),
        encoding="utf-8",
    )
    return root


class TestScannRunnerPlanning(unittest.TestCase):
    def test_large_auto_grid_uses_pipeline_tree_and_query_ladders(self) -> None:
        specs = frontier_specs(839_188, plan="frontier")
        self.assertTrue(specs)
        self.assertTrue(all(spec.mode == "tree_ah" for spec in specs))
        self.assertIn(916, {spec.num_leaves for spec in specs})
        self.assertIn(512, {spec.num_leaves for spec in specs})
        self.assertIn(2048, {spec.num_leaves for spec in specs})
        self.assertTrue(any(spec.leaves_to_search == 128 for spec in specs))
        self.assertTrue(any(spec.reorder == 800 for spec in specs))
        self.assertTrue(all(spec.ah_dim == 2 for spec in specs))

    def test_nonproduction_modes_are_refused_above_their_pipeline_size_ranges(self) -> None:
        with self.assertRaisesRegex(ValueError, "not exposed"):
            frontier_specs(839_188, plan="frontier", mode="brute-force")
        with self.assertRaisesRegex(ValueError, "production-representable"):
            frontier_specs(839_188, plan="frontier", mode="ah")

    def test_dry_run_validates_cache_and_forced_tree_without_scann(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary)
            write_prepared_cache(cache_dir)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                code = run(
                    parse_args(
                        [
                            "--cache-dir",
                            str(cache_dir),
                            "--dataset",
                            "fixture",
                            "--plan",
                            "smoke",
                            "--mode",
                            "tree-ah",
                            "--tree-leaves",
                            "16",
                            "--query-batch-size",
                            "16",
                            "--warmup-queries",
                            "16",
                            "--dry-run",
                        ]
                    )
                )
            self.assertEqual(code, 0)
            payload = json.loads(output.getvalue())
            self.assertEqual(payload["schema_version"], "ginflow-scann-benchmark-plan-v1")
            self.assertEqual(payload["dataset_window_count"], 768)
            self.assertTrue(payload["settings"])
            self.assertTrue(all(item["spec"]["mode"] == "tree_ah" for item in payload["settings"]))

    def test_rouskin_runs_require_three_repeats_even_for_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary)
            write_prepared_cache(cache_dir, dataset="rouskin_6k")
            with self.assertRaisesRegex(ValueError, "require --repeats >= 3"):
                run(
                    parse_args(
                        [
                            "--cache-dir",
                            str(cache_dir),
                            "--dataset",
                            "rouskin_6k",
                            "--plan",
                            "smoke",
                            "--warmup-queries",
                            "16",
                            "--query-batch-size",
                            "16",
                            "--repeats",
                            "2",
                        ]
                    )
                )


@unittest.skipUnless(importlib.util.find_spec("scann"), "run in the pinned ScaNN image")
class TestPinnedScannIntegration(unittest.TestCase):
    def test_smoke_tree_and_soar_write_valid_repeated_results(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary) / "cache"
            cache_root = write_prepared_cache(cache_dir)
            output_dir = Path(temporary) / "results"
            code = run(
                parse_args(
                    [
                        "--cache-dir",
                        str(cache_dir),
                        "--dataset",
                        "fixture",
                        "--output-dir",
                        str(output_dir),
                        "--plan",
                        "smoke",
                        "--mode",
                        "tree-ah",
                        "--tree-leaves",
                        "16",
                        "--include-soar",
                        "--repeats",
                        "2",
                        "--warmup-queries",
                        "16",
                        "--query-batch-size",
                        "16",
                        "--threads",
                        "1",
                        "--max-ram-fraction",
                        "1.0",
                    ]
                )
            )
            self.assertEqual(code, 0)
            records, issues, _ = read_records([output_dir])
            report = validate_records(records, read_issues=issues, min_repeats=2)
            self.assertTrue(report["valid"], report)
            self.assertGreaterEqual(report["successful_configurations"], 3, report)
            self.assertTrue(all(record.payload["backend"] == "scann" for record in records))
            self.assertTrue(all(record.payload["peak_vram_bytes"] is None for record in records))
            self.assertTrue(
                all(record.payload["provenance"]["scann_version"] == "1.4.2" for record in records)
            )
            scratch = cache_root / "scratch" / "scann"
            self.assertFalse(any(scratch.glob("scann-*")), list(scratch.glob("*")))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
