#!/usr/bin/env python3
"""Tests for exact original-window reranking."""
from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from rerank_candidates import main  # noqa: E402
from rerank_core import rerank_label_matrix  # noqa: E402


class TestRerankCandidates(unittest.TestCase):
    def test_cpu_workers_and_candidate_blocks_match_exact_order(self) -> None:
        database = np.asarray(
            [[1.0, 0.0], [0.8, 0.6], [0.0, 1.0], [-1.0, 0.0]], dtype=np.float32
        )
        queries = np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        labels = np.asarray([[3, 1, 0], [0, 1, 2]], dtype=np.int64)
        result_labels, result_scores, _ = rerank_label_matrix(
            database,
            queries,
            labels,
            2,
            batch_size=1,
            candidate_batch_size=1,
            workers=2,
        )
        np.testing.assert_array_equal(result_labels, [[0, 1], [2, 1]])
        np.testing.assert_allclose(result_scores, [[1.0, 0.8], [1.0, 0.6]])

    def test_cli_uses_original_residue_embeddings_and_writes_thresholded_seeds(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ginflow-rerank-test-") as tmp:
            root = Path(tmp)
            database = root / "database"
            database.mkdir()
            np.savez(
                database / "embeddings.npz",
                target=np.asarray([[1, 0], [0, 1], [1, 1]], dtype=np.float16),
            )
            (database / "windows.tsv").write_text(
                "faiss_id\ttranscript_id\tstart\tend\n"
                "0\ttarget\t0\t2\n"
                "1\ttarget\t1\t3\n"
            )
            query_windows = root / "query.windows.npz"
            np.savez(
                query_windows,
                query=np.asarray([[1, 0, 0, 1]], dtype=np.float32),
            )
            query_manifest = root / "query.windows.manifest.json"
            query_manifest.write_text(
                '{"stride": 1, "window_size": 2, "records": '
                '[{"identifier":"query", "n_windows":1}]}'
            )
            candidates = root / "candidates.tsv"
            with candidates.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["query_id", "query_start", "query_end", "target_id", "target_start", "target_end", "score", "rank"]
                )
                writer.writerow(["query", 0, 2, "target", 1, 3, 0.0, 1])
                writer.writerow(["query", 0, 2, "target", 0, 2, 0.0, 2])
            output = root / "seeds.tsv"
            metrics = root / "metrics.json"
            code = main(
                [
                    "--candidates", str(candidates),
                    "--database", str(database),
                    "--query-windows", str(query_windows),
                    "--query-manifests", str(query_manifest),
                    "--output", str(output),
                    "--metrics", str(metrics),
                    "--output-k", "2",
                    "--min-similarity", "0.9",
                    "--batch-size", "1",
                    "--candidate-batch-size", "1",
                    "--workers", "2",
                    "--device", "cpu",
                ]
            )
            self.assertEqual(code, 0)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["target_start"] for row in rows], ["0"])
            self.assertEqual(rows[0]["rank"], "1")
            self.assertEqual(json_load(metrics)["backend"], "exact_original_window_rerank")


def json_load(path: Path) -> dict:
    import json

    return json.loads(path.read_text())


if __name__ == "__main__":
    unittest.main()
