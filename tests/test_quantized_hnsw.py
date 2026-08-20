#!/usr/bin/env python3
"""Unit tests for quantized candidate windows and the HNSW score encoding."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from generate_quantized_windows import generate_windows, slice_codes  # noqa: E402
from build_hnswlib import build_database  # noqa: E402
from hnswlib_index import encoded_score_from_codes  # noqa: E402
from quantize_node_embeddings import (  # noqa: E402
    assign_nodes,
    centroid_similarity,
    fit_centroids,
    sample_training_nodes,
    write_quantization,
)
from search_hnswlib import search_windows  # noqa: E402


class TestQuantization(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(7)
        self.embeddings = rng.standard_normal((24, 8), dtype=np.float32)
        self.embeddings /= np.linalg.norm(self.embeddings, axis=1, keepdims=True)

    def test_centroids_and_similarity_are_stable(self) -> None:
        first, effective_k = fit_centroids(self.embeddings, requested_k=6, niter=4, seed=19)
        second, second_k = fit_centroids(self.embeddings, requested_k=6, niter=4, seed=19)
        self.assertEqual(effective_k, 6)
        self.assertEqual(second_k, effective_k)
        self.assertEqual(first.dtype, np.float16)
        self.assertEqual(first.shape, (6, 8))
        np.testing.assert_array_equal(first, second)
        similarity = centroid_similarity(first)
        self.assertEqual(similarity.dtype, np.float32)
        self.assertEqual(similarity.shape, (6, 6))
        np.testing.assert_allclose(similarity, similarity.T, atol=1e-6)
        np.testing.assert_allclose(np.diag(similarity), 1.0, atol=2e-3)

    def test_effective_k_and_codes(self) -> None:
        centroids, effective_k = fit_centroids(self.embeddings[:3], requested_k=8, niter=2, seed=2)
        self.assertEqual(effective_k, 3)
        codes = assign_nodes(self.embeddings, centroids)
        self.assertEqual(codes.shape, (24,))
        self.assertGreaterEqual(int(codes.min()), 0)
        self.assertLess(int(codes.max()), effective_k)

    def test_written_quantization_and_windows_keep_original_dtype(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            embedding_path = root / "shard.npz"
            manifest_path = root / "shard.manifest.json"
            source = self.embeddings.astype(np.float16)
            np.savez_compressed(embedding_path, record=source)
            manifest_path.write_text(
                json.dumps(
                    {
                        "records": [{"identifier": "record", "core_length": source.shape[0]}],
                        "checkpoint_sha256": "test",
                    }
                )
            )
            shards = [(embedding_path, manifest_path)]
            training, seen = sample_training_nodes(shards, sample_size=100, seed=3)
            self.assertEqual(seen, source.shape[0])
            centroids, _ = fit_centroids(training, requested_k=4, niter=3, seed=3)
            quantization_dir = root / "quantization"
            write_quantization(shards, quantization_dir, centroids, {"n_centroids": 4})
            with np.load(embedding_path) as arrays:
                self.assertEqual(arrays["record"].dtype, np.float16)
            windows_dir = root / "windows"
            result = generate_windows(quantization_dir, windows_dir, window_size=4, stride=2)
            self.assertEqual(result["n_windows"], 11)
            with np.load(windows_dir / "shard.windows.npz") as arrays:
                self.assertEqual(arrays["record"].shape, (11, 4))
                self.assertEqual(arrays["record"].dtype, np.uint16)

    def test_window_slicing(self) -> None:
        codes = np.arange(7, dtype=np.uint16)
        np.testing.assert_array_equal(
            slice_codes(codes, window_size=3, stride=2),
            np.array([[0, 1, 2], [2, 3, 4], [4, 5, 6]], dtype=np.uint16),
        )

    def test_encoded_score_is_the_registered_sum(self) -> None:
        rng = np.random.default_rng(11)
        centroids = rng.standard_normal((5, 8), dtype=np.float32)
        centroids /= np.linalg.norm(centroids, axis=1, keepdims=True)
        centroids = centroids.astype(np.float16)
        similarity = centroid_similarity(centroids)
        query = np.array([[0, 1, 2], [3, 4, 0]], dtype=np.uint16)
        target = np.array([[0, 2, 1], [4, 4, 0]], dtype=np.uint16)
        encoded, explicit = encoded_score_from_codes(query, target, centroids, similarity).T
        np.testing.assert_allclose(encoded, explicit, atol=2e-4)


@unittest.skipUnless(__import__("importlib.util").util.find_spec("hnswlib"), "hnswlib is not installed")
class TestHnswRoundTrip(unittest.TestCase):
    def test_database_build_and_search_preserve_original_embeddings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            embedding_path = root / "shard.npz"
            manifest_path = root / "shard.manifest.json"
            sidecar_path = root / "shard.json"
            source = self._source_embeddings()
            np.savez_compressed(embedding_path, record=source)
            manifest_path.write_text(
                json.dumps(
                    {
                        "records": [{"identifier": "record", "core_length": source.shape[0]}],
                        "checkpoint_sha256": "test",
                    }
                )
            )
            sidecar_path.write_text(
                json.dumps(
                    {
                        "identifiers": ["record"],
                        "sequences": ["A" * source.shape[0]],
                        "structures": ["." * source.shape[0]],
                    }
                )
            )
            centroids, _ = fit_centroids(source.astype(np.float32), requested_k=4, niter=3, seed=3)
            quantization_dir = root / "quantization"
            write_quantization(
                [(embedding_path, manifest_path)],
                quantization_dir,
                centroids,
                {"n_centroids": 4, "checkpoint_sha256": "test"},
            )
            windows_dir = root / "windows"
            generate_windows(quantization_dir, windows_dir, window_size=4, stride=1)
            database = root / "faiss"
            meta = build_database(
                windows_dir,
                quantization_dir,
                database,
                m=8,
                ef_construction=40,
                ef_search=20,
                random_seed=1,
                num_threads=1,
                embeddings=[embedding_path],
                graph_metadata=[sidecar_path],
            )
            self.assertEqual(meta["backend"], "hnswlib")
            self.assertTrue(meta["has_residue_embeddings"])
            with np.load(database / "embeddings.npz") as arrays:
                self.assertEqual(arrays["record"].dtype, np.float16)
                self.assertEqual(arrays["record"].shape, source.shape)
            hits = search_windows(windows_dir, database, 2, 0.0, 20, 1)
            self.assertGreater(len(hits), 0)
            self.assertEqual(hits[0][0], "record")
            self.assertEqual(hits[0][3], "record")
            self.assertGreater(float(hits[0][6]), 0.0)

    def _source_embeddings(self) -> np.ndarray:
        rng = np.random.default_rng(23)
        values = rng.standard_normal((24, 8), dtype=np.float32)
        values /= np.linalg.norm(values, axis=1, keepdims=True)
        return values.astype(np.float16)

    def test_index_round_trip(self) -> None:
        import hnswlib

        rng = np.random.default_rng(21)
        vectors = rng.standard_normal((12, 16), dtype=np.float32)
        index = hnswlib.Index(space="ip", dim=16)
        index.init_index(max_elements=len(vectors), M=8, ef_construction=40, random_seed=1)
        index.add_items(vectors, np.arange(len(vectors)))
        index.set_ef(20)
        labels, distances = index.knn_query(vectors[:2], k=2)
        self.assertEqual(labels.shape, (2, 2))
        np.testing.assert_allclose(1.0 - distances[:, 0], np.sum(vectors[:2] * vectors[labels[:, 0]], axis=1), atol=1e-5)


if __name__ == "__main__":
    unittest.main()
