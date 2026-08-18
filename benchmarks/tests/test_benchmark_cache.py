"""Small-fixture tests for the backend-neutral benchmark cache."""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

BENCHMARKS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCHMARKS))

from benchmark_utils import (  # noqa: E402
    LATENCY_UNIT,
    MEASUREMENT_PROTOCOL,
    QPS_SCOPE,
    RUNNER_VERSION,
    make_result_record,
    measure_search_repeats,
    recall_at_k,
    sha256_file,
    validate_result_record,
)
from prepare_cache import (  # noqa: E402
    PREPARATION_IMPLEMENTATION_PATHS,
    cache_windows_request,
    compute_exact_topk,
    flatten_window_shards,
    make_ground_truth_cache,
    make_query_cache,
    preparation_implementation_fingerprint,
    run_window_cache,
)


def unit_rows(rows: list[list[float]]) -> np.ndarray:
    values = np.asarray(rows, dtype=np.float32)
    values /= np.linalg.norm(values, axis=1, keepdims=True)
    return values


def write_fixture_shards(directory: Path) -> None:
    first = unit_rows([[1, 0, 0, 0], [0.9, 0.1, 0, 0], [0, 1, 0, 0]])
    second = unit_rows([[0, 0, 1, 0], [0, 0, 0.8, 0.2], [0, 0, 0, 1]])
    for name, arrays, records in (
        (
            "0000",
            {"short": first[:1], "long": first[1:]},
            [
                {"identifier": "short", "n_windows": 1},
                {"identifier": "long", "n_windows": 2},
            ],
        ),
        (
            "0001",
            {"third": second},
            [{"identifier": "third", "n_windows": 3}],
        ),
    ):
        np.savez_compressed(directory / f"{name}.windows.npz", **arrays)
        manifest = {
            "status": "complete",
            "window_size": 1,
            "stride": 1,
            "embedding_dim": 4,
            "window_dim": 4,
            "l2_normalized": True,
            "ginfinity_version": "fixture",
            "model_version": "fixture-model",
            "checkpoint_sha256": "a" * 64,
            "records": [
                {**record, "shape": [record["n_windows"], 4]}
                for record in records
            ],
            "skipped_short": [],
        }
        (directory / f"{name}.windows.manifest.json").write_text(json.dumps(manifest))


class TestCachePreparation(unittest.TestCase):
    def test_root_window_entrypoint_is_the_canonical_executable(self) -> None:
        repo_root = Path(__file__).resolve().parents[2]
        root_entrypoint = repo_root / "bin" / "slice_windows.py"
        self.assertTrue(root_entrypoint.is_file())
        self.assertFalse(root_entrypoint.is_symlink())
        self.assertTrue(os.access(root_entrypoint, os.X_OK))

    def test_flatten_then_record_balanced_selection_then_ground_truth(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            windows = root / "windows"
            windows.mkdir()
            write_fixture_shards(windows)
            flat = root / "flat"
            flatten = flatten_window_shards(windows, flat, dataset_id="fixture", row_chunk=1)
            self.assertEqual(flatten["vectors"]["shape"], [6, 4])
            vectors = np.load(flat / "vectors.npy", mmap_mode="r")
            np.testing.assert_allclose(np.linalg.norm(vectors, axis=1), np.ones(6), atol=1e-6)

            selection = make_query_cache(
                flat,
                root / "queries",
                dataset_id="fixture",
                query_count=3,
                seed="fixed",
            )
            selected_rows = (root / "queries" / "queries.tsv").read_text().splitlines()[1:]
            selected_records = {row.split("\t")[3] for row in selected_rows}
            self.assertEqual(selected_records, {"short", "long", "third"})
            self.assertEqual(selection["query_count"], 3)

            truth = make_ground_truth_cache(
                flat,
                root / "queries",
                root / "ground-truth",
                dataset_id="fixture",
                k=3,
                database_chunk=2,
                query_chunk=1,
                engine="cpu",
            )
            self.assertEqual(truth["engine"]["used"], "cpu")
            payload = np.load(root / "ground-truth" / "ground-truth.npz")
            self.assertEqual(payload["ids"].shape, (3, 3))
            self.assertTrue(np.all(payload["ids"][:, 0] >= 0))

    def test_chunked_topk_matches_direct_cosine(self) -> None:
        vectors = unit_rows(
            [
                [1, 0, 0, 0],
                [0.8, 0.2, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [0.1, 0.2, 0.3, 0.4],
                [0.2, 0.6, 0.1, 0.2],
            ]
        )
        queries = vectors[[0, 3, 6]]
        scores, identifiers, engine = compute_exact_topk(
            vectors,
            queries,
            k=4,
            database_chunk=3,
            query_chunk=2,
            engine="cpu",
        )
        expected_scores = queries @ vectors.T
        expected_ids = np.argsort(-expected_scores, axis=1, kind="stable")[:, :4]
        self.assertEqual(engine, "cpu")
        np.testing.assert_array_equal(identifiers, expected_ids)
        np.testing.assert_allclose(scores, np.take_along_axis(expected_scores, expected_ids, axis=1), atol=1e-6)

    def test_chunked_topk_keeps_lowest_ids_for_ties_across_blocks(self) -> None:
        vectors = unit_rows(
            [
                [1, 0, 0, 0],
                [1, 0, 0, 0],
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
            ]
        )
        scores, identifiers, engine = compute_exact_topk(
            vectors,
            vectors[:1],
            k=3,
            database_chunk=1,
            query_chunk=1,
            engine="cpu",
        )
        self.assertEqual(engine, "cpu")
        np.testing.assert_array_equal(identifiers, np.array([[0, 1, 2]], dtype=np.int64))
        np.testing.assert_allclose(scores, np.ones((1, 3), dtype=np.float32), atol=1e-6)

    def test_cache_ids_are_repeatable_and_relocatable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifests: list[tuple[dict[str, object], dict[str, object], dict[str, object]]] = []
            for cache_name in ("first", "relocated"):
                windows = root / cache_name / "windows"
                windows.mkdir(parents=True)
                write_fixture_shards(windows)
                flat = root / cache_name / "flat"
                flatten = flatten_window_shards(windows, flat, dataset_id="fixture", row_chunk=1)
                queries = make_query_cache(flat, root / cache_name / "queries", dataset_id="fixture", query_count=3, seed="fixed")
                truth = make_ground_truth_cache(
                    flat,
                    root / cache_name / "queries",
                    root / cache_name / "ground-truth",
                    dataset_id="fixture",
                    k=3,
                    database_chunk=2,
                    query_chunk=1,
                    engine="cpu",
                )
                manifests.append((flatten, queries, truth))

            first, relocated = manifests
            self.assertEqual(first[0]["embedding_cache_id"], relocated[0]["embedding_cache_id"])
            self.assertEqual(first[1]["query_selection_id"], relocated[1]["query_selection_id"])
            self.assertEqual(first[2]["ground_truth_cache_id"], relocated[2]["ground_truth_cache_id"])

            # Regenerating the same content must not let timestamps or elapsed
            # time alter any compatibility key.
            flat = root / "first" / "flat"
            queries_dir = root / "first" / "queries"
            truth_dir = root / "first" / "ground-truth"
            repeated_flatten = flatten_window_shards(root / "first" / "windows", flat, dataset_id="fixture", row_chunk=1, overwrite=True)
            repeated_queries = make_query_cache(flat, queries_dir, dataset_id="fixture", query_count=3, seed="fixed", overwrite=True)
            repeated_truth = make_ground_truth_cache(
                flat,
                queries_dir,
                truth_dir,
                dataset_id="fixture",
                k=3,
                database_chunk=2,
                query_chunk=1,
                engine="cpu",
                overwrite=True,
            )
            self.assertEqual(first[0]["embedding_cache_id"], repeated_flatten["embedding_cache_id"])
            self.assertEqual(first[1]["query_selection_id"], repeated_queries["query_selection_id"])
            self.assertEqual(first[2]["ground_truth_cache_id"], repeated_truth["ground_truth_cache_id"])

    def test_preparation_fingerprint_tracks_live_runtime_files_not_checkout_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repo = root / "repo"
            for relative in PREPARATION_IMPLEMENTATION_PATHS:
                path = repo / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"fixture for {relative}\n", encoding="utf-8")
            runtime_script = repo / "bin" / "slice_windows.py"
            runtime_script.chmod(0o755)
            source_a = root / "source-a.tsv"
            source_b = root / "elsewhere" / "source-b.tsv"
            source_b.parent.mkdir()
            source_a.write_text("transcript_id\tsequence\tsecondary_structure\na\tA\t.\n", encoding="utf-8")
            source_b.write_bytes(source_a.read_bytes())

            initial = preparation_implementation_fingerprint(repo, "docker,gpu")
            request_a = cache_windows_request(
                repo, source_a, "fixture", profile="docker,gpu", shard_size=50, window_size=11, window_stride=1
            )
            request_b = cache_windows_request(
                repo, source_b, "fixture", profile="docker,gpu", shard_size=50, window_size=11, window_stride=1
            )
            self.assertEqual(request_a["identity"], request_b["identity"])
            self.assertNotEqual(request_a["input_path"], request_b["input_path"])

            runtime_script.write_text("changed implementation\n", encoding="utf-8")
            changed = preparation_implementation_fingerprint(repo, "docker,gpu")
            self.assertNotEqual(initial["id"], changed["id"])
            changed_request = cache_windows_request(
                repo, source_a, "fixture", profile="docker,gpu", shard_size=50, window_size=11, window_stride=1
            )
            self.assertNotEqual(request_a["identity"], changed_request["identity"])

    def test_cache_windows_dry_run_is_side_effect_free(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "fixture.tsv"
            source.write_text("transcript_id\tsequence\tsecondary_structure\na\tA\t.\n")
            repo_root = Path(__file__).resolve().parents[2]
            response = run_window_cache(
                repo_root=repo_root,
                cache_dir=root / "cache",
                dataset_id="fixture",
                input_path=source,
                profile="docker,gpu",
                shard_size=50,
                window_size=11,
                window_stride=1,
                nextflow="nextflow",
                dry_run=True,
            )
            self.assertEqual(response["status"], "planned")
            command = response["execution"]["command"]
            self.assertEqual(command[:4], ["nextflow", "run", str(repo_root / "benchmarks" / "cache_windows.nf"), "-c"])
            self.assertNotIn("-main-script", command)
            self.assertFalse((root / "cache" / "fixture" / "window-cache-request.json").exists())

CACHE_WINDOW_DOCKER_INTEGRATION_READY = (
    os.environ.get("GINFLOW_RUN_NEXTFLOW_CACHE_INTEGRATION") == "1"
    and shutil.which("nextflow") is not None
    and shutil.which("docker") is not None
)


@unittest.skipUnless(CACHE_WINDOW_DOCKER_INTEGRATION_READY, "set GINFLOW_RUN_NEXTFLOW_CACHE_INTEGRATION=1 with Docker")
class TestPinnedCacheWindowIntegration(unittest.TestCase):
    def test_root_window_script_runs_and_completed_cache_is_reused(self) -> None:
        repo_root = Path(__file__).resolve().parents[2]
        source = repo_root / "benchmarks" / "tests" / "fixtures" / "cache_windows_tiny.tsv"
        with tempfile.TemporaryDirectory() as temporary:
            cache_dir = Path(temporary) / "cache"
            common = {
                "repo_root": repo_root,
                "cache_dir": cache_dir,
                "dataset_id": "tiny-window-script",
                "input_path": source,
                "profile": "docker,gpu",
                "shard_size": 50,
                "window_size": 11,
                "window_stride": 1,
                "nextflow": "nextflow",
            }
            created = run_window_cache(**common)
            self.assertEqual(created["status"], "created")
            windows = cache_dir / "tiny-window-script" / "artifacts" / "windows"
            self.assertEqual(len(list(windows.glob("*.windows.npz"))), 1)
            reused = run_window_cache(**common)
            self.assertEqual(reused["status"], "reused")


class TestResultSchema(unittest.TestCase):
    def test_fixed_batch_measurement_reports_real_latency_distribution(self) -> None:
        queries = np.arange(16, dtype=np.float32).reshape(8, 2)
        truth = np.tile(np.array([[0, 1]], dtype=np.int64), (8, 1))
        calls: list[int] = []

        def search(batch: np.ndarray, k: int) -> np.ndarray:
            calls.append(batch.shape[0])
            time.sleep(0.001 + float(batch[0, 0]) * 0.0003)
            return np.tile(np.arange(k, dtype=np.int64), (batch.shape[0], 1))

        measurement = next(
            measure_search_repeats(
                search,
                queries,
                truth,
                k=2,
                warmup_queries=2,
                repeats=1,
                query_batch_size=2,
            )
        )
        self.assertEqual(calls, [2, 2, 2, 2, 2])  # one warm-up plus four timed batches
        self.assertEqual(measurement.measurement["protocol"], MEASUREMENT_PROTOCOL)
        self.assertEqual(measurement.measurement["query_batch_size"], 2)
        self.assertEqual(measurement.measurement["timed_batch_count"], 4)
        self.assertEqual(measurement.measurement["latency_unit"], LATENCY_UNIT)
        self.assertEqual(measurement.measurement["qps_scope"], QPS_SCOPE)
        self.assertGreater(measurement.latency_ms["p95"], measurement.latency_ms["p50"])
        self.assertAlmostEqual(
            measurement.qps,
            measurement.measurement["timed_queries"] / measurement.measurement["timed_seconds"],
            places=6,
        )

    def test_recall_and_successful_result_contract(self) -> None:
        labels = np.array([[0, 1, 2], [5, 4, -1]], dtype=np.int64)
        truth = np.array([[2, 1, 0], [4, 5, 3]], dtype=np.int64)
        self.assertAlmostEqual(recall_at_k(labels, truth, k=3), (1.0 + 2 / 3) / 2)
        record = make_result_record(
            backend="fixture",
            dataset_id="fixture",
            dataset_window_count=6,
            dimension=4,
            parameter_label="fixture",
            parameters={"example": True},
            run_id="fixture-run",
            repeat=0,
            warmup_queries=1,
            timed_queries=3,
            query_ids_sha256="a" * 64,
            ground_truth_ids_sha256="b" * 64,
            provenance={
                "git_commit": "deadbeef",
                "runner_version": RUNNER_VERSION,
                "hardware_id": "fixture-machine",
                "embedding_cache_id": "embedding-fixture",
                "query_selection_id": "queries-fixture",
                "ground_truth_cache_id": "ground-truth-fixture",
            },
            qps=100.0,
            latency_ms={"mean": 10.0, "p50": 10.0, "p95": 10.0},
            recall=0.9,
            index_bytes=123,
            build_seconds=0.2,
            peak_rss_bytes=456,
            peak_vram_bytes=None,
            measurement={
                "protocol": MEASUREMENT_PROTOCOL,
                "query_batch_size": 1,
                "timed_batch_count": 3,
                "timed_queries": 3,
                "latency_unit": LATENCY_UNIT,
                "qps_scope": QPS_SCOPE,
            },
        )
        self.assertEqual(validate_result_record(record), [])


if __name__ == "__main__":
    raise SystemExit(unittest.main())
