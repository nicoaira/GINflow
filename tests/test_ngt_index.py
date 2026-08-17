#!/usr/bin/env python3
"""CPU tests for the NGT index variants."""
from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from ngt_index import (  # noqa: E402
    NGT_INDEX_TYPES,
    build_populated_index,
    is_ngt_type,
    load_index,
    normalize_index_type,
    serialize_index,
)


def unit_vectors(n: int, dim: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vectors = rng.standard_normal((n, dim), dtype=np.float32)
    vectors /= np.linalg.norm(vectors, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(vectors)


class TestNgtHelpers(unittest.TestCase):
    def test_aliases(self) -> None:
        self.assertEqual(normalize_index_type("ngt"), "NGT")
        self.assertEqual(normalize_index_type("qg"), "QG")
        self.assertEqual(normalize_index_type("qbg"), "QBG")
        self.assertEqual(normalize_index_type("anng"), "NGT")
        self.assertEqual(normalize_index_type("quantized-blob-graph"), "QBG")
        self.assertTrue(is_ngt_type("qg"))
        self.assertFalse(is_ngt_type("FlatIP"))
        self.assertEqual(set(NGT_INDEX_TYPES), {"NGT", "QG", "QBG"})


NGT_READY = shutil.which("qbg") is not None
try:
    import ngtpy  # noqa: F401
except ImportError:
    NGT_READY = False


@unittest.skipUnless(NGT_READY, "ngtpy and qbg are required")
class TestNgtIndexes(unittest.TestCase):
    def setUp(self) -> None:
        self.vectors = unit_vectors(100, 4, seed=4)

    def test_variants_search_and_roundtrip(self) -> None:
        for kind in NGT_INDEX_TYPES:
            with self.subTest(kind=kind):
                index, details = build_populated_index(self.vectors, kind.lower())
                self.assertEqual(details["backend"], "ngt")
                self.assertEqual(details["index_type"], kind)
                distances, labels = index.search(self.vectors[:3], 5)
                self.assertEqual(labels.shape, (3, 5))
                self.assertTrue(np.all(labels[:, 0] >= 0))
                self.assertTrue(np.all(np.isfinite(distances)))
                with tempfile.TemporaryDirectory() as tmp:
                    destination = Path(tmp) / "ngt"
                    serialize_index(index, destination)
                    loaded = load_index(destination, {"index_type": kind, "n_windows": 100})
                    loaded_distances, loaded_labels = loaded.search(self.vectors[:1], 5)
                    self.assertEqual(loaded_labels.shape, (1, 5))
                    self.assertEqual(int(loaded_labels[0, 0]), int(labels[0, 0]))
                    self.assertTrue(np.isfinite(loaded_distances[0, 0]))
                    loaded.close()


if __name__ == "__main__":
    raise SystemExit(unittest.main())
