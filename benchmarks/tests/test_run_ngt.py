"""Contract tests for the pinned NGT benchmark runner."""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


BENCHMARKS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCHMARKS))

from benchmark_utils import sha256_array  # noqa: E402
from results_common import read_records  # noqa: E402
from run_ngt import NgtSpec, capacity_check, frontier_specs, parse_args, run  # noqa: E402
from validate_results import validate_records  # noqa: E402


def _unit_rows(n_rows: int, dimension: int, seed: int = 37) -> np.ndarray:
    generator = np.random.default_rng(seed)
    values = generator.standard_normal((n_rows, dimension), dtype=np.float32)
    values /= np.linalg.norm(values, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(values, dtype=np.float32)


def write_prepared_cache(cache_dir: Path, *, dataset: str = "fixture") -> Path:
    """Write a small cache following the backend-neutral cache contract."""

    root = cache_dir / dataset
    flat = root / "flat"
    queries_dir = root / "queries"
    truth_dir = root / "ground-truth"
    for directory in (flat, queries_dir, truth_dir):
        directory.mkdir(parents=True)
    vectors = _unit_rows(1_024, 64)
    query_vector_ids = np.arange(0, 1_024, 32, dtype=np.int64)[:32]
    queries = np.ascontiguousarray(vectors[query_vector_ids], dtype=np.float32)
    scores = queries @ vectors.T
    identifiers = np.argsort(-scores, axis=1, kind="stable")[:, :100].astype(np.int64)
    embedding_cache_id = "embedding-ngt-fixture"
    query_selection_id = "query-ngt-fixture"
    ground_truth_cache_id = "truth-ngt-fixture"
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
        )
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
        )
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
        )
    )
    return root


class TestNgtRunnerPlanning(unittest.TestCase):
    def test_pipeline_specs_are_exactly_the_three_exposed_structures(self) -> None:
        specs = frontier_specs(839_188, 1_408, plan="frontier")
        self.assertEqual([(spec.mode, spec.kind) for spec in specs], [("pipeline", "NGT"), ("pipeline", "QG"), ("pipeline", "QBG")])
        self.assertTrue(all(spec.parameters(query_batch_size=32)["pipeline_faithful"] for spec in specs))

    def test_native_frontier_is_opt_in_and_labelled(self) -> None:
        native = frontier_specs(839_188, 1_408, plan="smoke", include_native_frontier=True)[3:]
        self.assertTrue(native)
        self.assertTrue(all(spec.mode == "native" for spec in native))
        self.assertTrue(any(spec.kind == "NGT" and spec.epsilon is not None for spec in native))
        self.assertTrue(any(spec.kind == "QG" and spec.result_expansion is not None for spec in native))
        self.assertTrue(any(spec.kind == "QBG" and spec.exploration_size is not None for spec in native))

    def test_pipeline_spec_rejects_hidden_parameter_override(self) -> None:
        with self.assertRaises(ValueError):
            NgtSpec("NGT", epsilon=0.1)
        with self.assertRaises(ValueError):
            NgtSpec("QG", "native", qg_max_edges=15)

    def test_capacity_gate_refuses_qbg_when_workspace_is_not_safe(self) -> None:
        spec = NgtSpec("QBG")
        result = capacity_check(
            spec,
            n_vectors=839_188,
            dimension=1_408,
            available_ram_bytes=32 * 1024**3,
            free_disk_bytes=1_000 * 1024**3,
            max_ram_fraction=0.70,
            max_disk_fraction=0.80,
        )
        self.assertFalse(result.feasible)
        self.assertIn("MemAvailable", str(result.reason))

    def test_dry_run_validates_cache_without_ngt(self) -> None:
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
                            "--warmup-queries",
                            "16",
                            "--query-batch-size",
                            "16",
                            "--dry-run",
                        ]
                    )
                )
            self.assertEqual(code, 0)
            payload = json.loads(output.getvalue())
            self.assertEqual(payload["schema_version"], "ginflow-ngt-benchmark-plan-v1")
            self.assertEqual(payload["dataset_window_count"], 1_024)
            self.assertEqual(len(payload["settings"]), 3)


NGT_READY = shutil.which("qbg") is not None and importlib.util.find_spec("ngtpy") is not None


@unittest.skipUnless(NGT_READY, "run in the pinned NGT 2.3.12 image")
class TestPinnedNgtIntegration(unittest.TestCase):
    def test_smoke_pipeline_frontier_writes_valid_repeats(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary) / "cache"
            write_prepared_cache(cache_dir)
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
                        "--repeats",
                        "2",
                        "--warmup-queries",
                        "16",
                        "--query-batch-size",
                        "16",
                        "--max-ram-fraction",
                        "1.0",
                        "--max-disk-fraction",
                        "1.0",
                    ]
                )
            )
            self.assertEqual(code, 0)
            records, issues, _ = read_records([output_dir])
            report = validate_records(records, read_issues=issues, min_repeats=2)
            self.assertTrue(report["valid"], report)
            # NGT and QG should be runnable on the tiny fixture.  QBG may be
            # deliberately skipped on a constrained host due its 24 GiB gate.
            self.assertGreaterEqual(report["successful_configurations"], 2, report)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
