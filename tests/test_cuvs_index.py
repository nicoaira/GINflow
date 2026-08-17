#!/usr/bin/env python3
"""GPU tests for the cuVS CAGRA, IVF-Flat, and IVF-PQ indexes."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import cuvs_index  # noqa: E402
from faiss_index import distances_to_similarity  # noqa: E402
from cuvs_index import (  # noqa: E402
    CUVS_INDEX_TYPES,
    build_populated_index,
    is_cuvs_type,
    load_index,
    normalize_index_type,
    serialize_index,
)


def unit_vectors(n: int, dim: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vectors = rng.standard_normal((n, dim), dtype=np.float32)
    vectors /= np.linalg.norm(vectors, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(vectors)


class TestCuvsHelpers(unittest.TestCase):
    def test_aliases(self) -> None:
        self.assertEqual(normalize_index_type("cagra"), "CAGRA")
        self.assertEqual(normalize_index_type("ivf-flat"), "IVF")
        self.assertEqual(normalize_index_type("ivfpq"), "IVF-PQ")
        self.assertTrue(is_cuvs_type("ivf"))
        self.assertFalse(is_cuvs_type("flatip"))
        self.assertEqual(set(CUVS_INDEX_TYPES), {"CAGRA", "IVF", "IVF-PQ"})

    def test_bounds(self) -> None:
        self.assertEqual(cuvs_index.resolve_n_lists(100, 4), 4)
        self.assertEqual(cuvs_index.resolve_n_probes(4, 20), 4)
        with self.assertRaises(ValueError):
            cuvs_index.resolve_n_lists(4, 5)
        with self.assertRaises(ValueError):
            cuvs_index.resolve_degree(1, 4, "cuvs_graph_degree")

    def test_cosine_distance_conversion(self) -> None:
        np.testing.assert_allclose(
            distances_to_similarity(np.array([[0.0, 0.25]], dtype=np.float32), "cuvs_cosine"),
            [[1.0, 0.75]],
        )


def gpu_ready() -> bool:
    return cuvs_index.gpu_available()


@unittest.skipUnless(gpu_ready(), "cuVS and an NVIDIA GPU are required")
class TestCuvsIndexes(unittest.TestCase):
    def setUp(self) -> None:
        self.vectors = unit_vectors(256, 32, seed=4)

    def _roundtrip(self, kind: str, **options: int | str) -> None:
        index, details = build_populated_index(self.vectors, kind, **options)
        self.assertEqual(details["backend"], "cuvs")
        self.assertEqual(details["index_type"], normalize_index_type(kind))
        distances, labels = index.search(self.vectors[:3], 5)
        self.assertEqual(labels.shape, (3, 5))
        self.assertEqual(labels[:, 0].tolist(), [0, 1, 2])
        self.assertTrue(np.all(distances[:, 0] < 1e-5))
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "cuvs"
            serialize_index(index, destination)
            loaded = load_index(
                destination,
                {
                    "index_type": kind,
                    "n_windows": self.vectors.shape[0],
                    "n_lists": details["n_lists"],
                    "n_probes": details["n_probes"],
                    "itopk_size": details["itopk_size"],
                },
            )
            loaded_distances, loaded_labels = loaded.search(self.vectors[:1], 5)
            self.assertEqual(int(loaded_labels[0, 0]), 0)
            self.assertTrue(np.isfinite(loaded_distances[0, 0]))
            loaded.close()

    def test_cagra(self) -> None:
        self._roundtrip(
            "cagra",
            intermediate_graph_degree=64,
            graph_degree=32,
            build_algo="nn_descent",
            itopk_size=4,
        )

    def test_ivf_flat(self) -> None:
        self._roundtrip("ivf", n_lists=4, n_probes=4)

    def test_ivf_pq(self) -> None:
        self._roundtrip("ivf-pq", n_lists=4, n_probes=4, pq_bits=8, pq_dim=8)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
