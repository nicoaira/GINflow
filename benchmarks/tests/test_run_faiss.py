"""Contract and pinned-image integration tests for the FAISS benchmark runner."""
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
from run_faiss import frontier_specs, parse_args, run  # noqa: E402
from validate_results import validate_records  # noqa: E402


def _unit_rows(n_rows: int, dimension: int, seed: int = 19) -> np.ndarray:
    generator = np.random.default_rng(seed)
    values = generator.standard_normal((n_rows, dimension), dtype=np.float32)
    values /= np.linalg.norm(values, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(values, dtype=np.float32)


def write_prepared_cache(cache_dir: Path, *, dataset: str = "fixture") -> Path:
    """Make a small cache that follows the shared on-disk benchmark contract."""

    root = cache_dir / dataset
    flat = root / "flat"
    queries_dir = root / "queries"
    truth_dir = root / "ground-truth"
    for directory in (flat, queries_dir, truth_dir):
        directory.mkdir(parents=True)
    vectors = _unit_rows(10_240, 16)
    queries = np.ascontiguousarray(vectors[np.arange(0, 10_240, 320)[:32]], dtype=np.float32)
    scores = queries @ vectors.T
    identifiers = np.argsort(-scores, axis=1, kind="stable")[:, :100].astype(np.int64)
    query_vector_ids = np.arange(0, 10_240, 320, dtype=np.int64)[:32]
    embedding_cache_id = "embedding-fixture"
    query_selection_id = "query-fixture"
    ground_truth_cache_id = "truth-fixture"
    np.save(flat / "vectors.npy", vectors, allow_pickle=False)
    np.save(queries_dir / "queries.npy", queries, allow_pickle=False)
    np.savez_compressed(truth_dir / "ground-truth.npz", ids=identifiers, scores=np.take_along_axis(scores, identifiers, axis=1))
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


class TestFaissRunnerPlanning(unittest.TestCase):
    def test_frontier_uses_current_pipeline_families_and_divisible_pq(self) -> None:
        specs = frontier_specs(839_188, 1_408, plan="frontier", device="cpu")
        self.assertEqual(specs[0].index_type, "FlatIP")
        self.assertTrue(any(spec.index_type == "IVFFlat" for spec in specs))
        self.assertTrue(any(spec.index_type == "IVFPQ" for spec in specs))
        self.assertTrue(any(spec.index_type == "HNSW" for spec in specs))
        for spec in specs:
            if spec.index_type == "IVFPQ":
                self.assertEqual(1_408 % int(spec.pq_m), 0)
                self.assertIn(spec.pq_nbits, {4, 8})
            if spec.nlist is not None:
                self.assertLessEqual(spec.nlist, 839_188)

    def test_dry_run_validates_prepared_cache_without_faiss(self) -> None:
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
            self.assertEqual(payload["schema_version"], "ginflow-faiss-benchmark-plan-v1")
            self.assertEqual(payload["dataset_window_count"], 10_240)
            self.assertTrue(payload["settings"])

    def test_rouskin_runs_require_three_repeats(self) -> None:
        with self.assertRaisesRegex(ValueError, "require --repeats >= 3"):
            run(parse_args(["--dataset", "rouskin_6k", "--repeats", "2"]))


@unittest.skipUnless(importlib.util.find_spec("faiss"), "run in the pinned FAISS image")
class TestPinnedFaissIntegration(unittest.TestCase):
    def test_smoke_frontier_writes_two_valid_repeats_per_setting(self) -> None:
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
                        "--threads",
                        "1",
                        "--add-batch-size",
                        "64",
                        "--train-sample-size",
                        "9984",
                        "--max-ram-fraction",
                        "1.0",
                    ]
                )
            )
            self.assertEqual(code, 0)
            records, issues, _ = read_records([output_dir])
            report = validate_records(records, read_issues=issues, min_repeats=2)
            self.assertTrue(report["valid"], report)
            self.assertGreaterEqual(report["successful_configurations"], 8)


def gpu_faiss_ready() -> bool:
    if not importlib.util.find_spec("faiss"):
        return False
    import faiss

    try:
        return bool(hasattr(faiss, "StandardGpuResources") and faiss.get_num_gpus() > 0)
    except RuntimeError:
        return False


@unittest.skipUnless(gpu_faiss_ready(), "run in the pinned FAISS GPU image with a visible CUDA device")
class TestPinnedGpuFaissIntegration(unittest.TestCase):
    def test_smoke_gpu_frontier_writes_valid_repeats(self) -> None:
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
                        "--device",
                        "gpu",
                        "--plan",
                        "smoke",
                        "--repeats",
                        "2",
                        "--warmup-queries",
                        "16",
                        "--query-batch-size",
                        "16",
                        "--threads",
                        "1",
                        "--add-batch-size",
                        "64",
                        "--train-sample-size",
                        "9984",
                        "--resource-temp-mib",
                        "32",
                        "--max-vram-fraction",
                        "0.9",
                    ]
                )
            )
            self.assertEqual(code, 0)
            records, issues, _ = read_records([output_dir])
            report = validate_records(records, read_issues=issues, min_repeats=2)
            self.assertTrue(report["valid"], report)
            self.assertGreaterEqual(report["successful_configurations"], 5, report)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
