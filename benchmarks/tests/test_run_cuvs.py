"""Planning and pinned-image integration checks for the cuVS benchmark runner."""
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
from run_cuvs import (  # noqa: E402
    GIB,
    capacity_check,
    frontier_specs,
    group_build_specs,
    parse_args,
    run,
)
from validate_results import validate_records  # noqa: E402


def _unit_rows(n_rows: int, dimension: int, seed: int = 37) -> np.ndarray:
    generator = np.random.default_rng(seed)
    values = generator.standard_normal((n_rows, dimension), dtype=np.float32)
    values /= np.linalg.norm(values, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(values, dtype=np.float32)


def write_prepared_cache(cache_dir: Path, *, dataset: str = "fixture", rows: int = 10_240) -> Path:
    """Make the shared cache contract without requiring a GINFINITY run."""

    root = cache_dir / dataset
    flat = root / "flat"
    queries_dir = root / "queries"
    truth_dir = root / "ground-truth"
    for directory in (flat, queries_dir, truth_dir):
        directory.mkdir(parents=True)
    vectors = _unit_rows(rows, 32)
    query_vector_ids = np.arange(0, rows, max(1, rows // 32), dtype=np.int64)[:32]
    queries = np.ascontiguousarray(vectors[query_vector_ids], dtype=np.float32)
    scores = queries @ vectors.T
    identifiers = np.argsort(-scores, axis=1, kind="stable")[:, :100].astype(np.int64)
    embedding_cache_id = "embedding-fixture"
    query_selection_id = "query-fixture"
    ground_truth_cache_id = "truth-fixture"
    np.save(flat / "vectors.npy", vectors, allow_pickle=False)
    np.save(queries_dir / "queries.npy", queries, allow_pickle=False)
    np.savez_compressed(truth_dir / "ground-truth.npz", ids=identifiers)
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


class TestCuvsRunnerPlanning(unittest.TestCase):
    def test_frontier_covers_exposed_build_and_query_controls(self) -> None:
        specs = frontier_specs(839_188, 1_408, plan="frontier")
        cagra = [spec for spec in specs if spec.index_type == "CAGRA"]
        ivf = [spec for spec in specs if spec.index_type == "IVF"]
        pq = [spec for spec in specs if spec.index_type == "IVF-PQ"]
        self.assertEqual({spec.itopk_size for spec in cagra}, {100, 128, 256})
        self.assertTrue({(spec.intermediate_graph_degree, spec.graph_degree) for spec in cagra}.issuperset({(128, 64), (256, 64), (256, 128)}))
        self.assertEqual({spec.n_lists for spec in ivf}, {2_048, 4_096})
        self.assertTrue({spec.n_probes for spec in ivf}.issuperset({8, 16, 32, 64, 128}))
        self.assertTrue({spec.pq_bits for spec in pq}.issuperset({4, 6, 8}))
        self.assertIn(256, {spec.pq_dim for spec in pq})
        for build, searches in group_build_specs(ivf):
            self.assertTrue(all(spec.build_key() == build.build_key() for spec in searches))
            self.assertEqual({spec.build_n_probes for spec in searches}, {8})
            self.assertGreaterEqual(len(searches), 5)
        for spec in cagra:
            parameters = spec.parameters(query_batch_size=32)
            self.assertFalse(parameters["search_override_controls"])
            self.assertIn("cuvs_itopk_size", parameters["build_persisted_controls"])

    def test_30k_raw_matrix_gate_is_explicit_on_an_8_gib_gpu(self) -> None:
        spec = frontier_specs(4_115_576, 1_408, plan="frontier")[0]
        capacity = capacity_check(
            spec,
            n_vectors=4_115_576,
            dimension=1_408,
            available_vram_bytes=8 * GIB,
            available_ram_bytes=64 * GIB,
            max_vram_fraction=0.9,
            max_ram_fraction=0.8,
        )
        self.assertFalse(capacity.feasible)
        self.assertGreater(capacity.raw_vector_bytes, 8 * GIB)
        self.assertIn("complete raw float32 matrix", str(capacity.reason))

    def test_dry_run_validates_cache_without_cuvs(self) -> None:
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
            self.assertEqual(payload["schema_version"], "ginflow-cuvs-benchmark-plan-v1")
            self.assertEqual(payload["dataset_window_count"], 10_240)
            self.assertTrue(payload["settings"])

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


def cuvs_gpu_ready() -> bool:
    if importlib.util.find_spec("cuvs") is None or importlib.util.find_spec("cupy") is None:
        return False
    try:
        import cupy as cp

        return int(cp.cuda.runtime.getDeviceCount()) > 0
    except RuntimeError:
        return False


@unittest.skipUnless(cuvs_gpu_ready(), "run in the pinned cuVS image with a visible CUDA device")
class TestPinnedCuvsIntegration(unittest.TestCase):
    def test_smoke_ivf_writes_valid_fixed_batch_repeats(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary) / "cache"
            write_prepared_cache(cache_dir, rows=2_048)
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
                        "--only",
                        "ivf",
                        "--repeats",
                        "2",
                        "--warmup-queries",
                        "16",
                        "--query-batch-size",
                        "16",
                        "--max-ram-fraction",
                        "1.0",
                        "--max-vram-fraction",
                        "1.0",
                    ]
                )
            )
            self.assertEqual(code, 0)
            records, issues, _ = read_records([output_dir])
            report = validate_records(records, read_issues=issues, min_repeats=2)
            self.assertTrue(report["valid"], report)
            self.assertGreaterEqual(report["successful_configurations"], 2, report)

    def test_smoke_cagra_and_ivf_pq_use_production_roundtrips(self) -> None:
        """Exercise the two other production adapter families through the runner."""

        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary) / "cache"
            write_prepared_cache(cache_dir, rows=2_048)
            output_dir = Path(temporary) / "results"
            for family in ("cagra", "ivf-pq"):
                code = run(
                    parse_args(
                        [
                            "--cache-dir",
                            str(cache_dir),
                            "--dataset",
                            "fixture",
                            "--output-dir",
                            str(output_dir / family),
                            "--plan",
                            "smoke",
                            "--only",
                            family,
                            "--repeats",
                            "1",
                            "--warmup-queries",
                            "16",
                            "--query-batch-size",
                            "16",
                            "--max-ram-fraction",
                            "1.0",
                            "--max-vram-fraction",
                            "1.0",
                        ]
                    )
                )
                self.assertEqual(code, 0, family)
                records, issues, _ = read_records([output_dir / family])
                report = validate_records(records, read_issues=issues, min_repeats=1)
                self.assertTrue(report["valid"], (family, report))
                self.assertGreaterEqual(report["successful_configurations"], 1, (family, report))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
