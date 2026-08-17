#!/usr/bin/env python3
"""Unit tests for the ScaNN window index (skipped when scann is not installed)."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from faiss_index import IndexOptions  # noqa: E402
from scann_index import (  # noqa: E402
    apply_search_params,
    build_populated_searcher,
    is_scann_database,
    is_scann_type,
    load_index,
    resolve_scann_plan,
    serialize_index,
)

try:
    import scann  # noqa: F401

    SCANN_AVAILABLE = True
except ImportError:
    SCANN_AVAILABLE = False


def unit_vectors(n: int, dim: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    xb = rng.standard_normal((n, dim), dtype=np.float32)
    xb /= np.linalg.norm(xb, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(xb)


class TestScannHelpers(unittest.TestCase):
    def test_type_aliases(self) -> None:
        self.assertTrue(is_scann_type("ScaNN"))
        self.assertTrue(is_scann_type("scann"))
        self.assertTrue(is_scann_type("scan"))
        self.assertFalse(is_scann_type("FlatIP"))
        self.assertFalse(is_scann_type(None))

    def test_database_detection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.assertFalse(is_scann_database(root, {"index_type": "FlatIP"}))
            self.assertTrue(is_scann_database(root, {"index_type": "ScaNN"}))
            self.assertTrue(is_scann_database(root, {"backend": "scann"}))
            (root / "scann").mkdir()
            self.assertTrue(is_scann_database(root, {"index_type": "FlatIP"}))

    def test_auto_plan(self) -> None:
        small = resolve_scann_plan(500, IndexOptions())
        self.assertEqual(small["scoring"], "brute_force")
        self.assertFalse(small["partitioned"])
        forced = resolve_scann_plan(500, IndexOptions(scann_leaves=4, scann_leaves_to_search=2))
        self.assertEqual(forced["scoring"], "ah")
        self.assertTrue(forced["partitioned"])
        self.assertEqual(forced["num_leaves"], 4)
        self.assertEqual(forced["leaves_to_search"], 2)
        mid = resolve_scann_plan(50_000, IndexOptions())
        self.assertEqual(mid["scoring"], "ah")
        self.assertFalse(mid["partitioned"])
        large = resolve_scann_plan(200_000, IndexOptions())
        self.assertEqual(large["scoring"], "ah")
        self.assertTrue(large["partitioned"])
        self.assertGreaterEqual(large["num_leaves"], 1)

    def test_gpu_rejected(self) -> None:
        xb = unit_vectors(8, 16)
        with self.assertRaises(ValueError) as ctx:
            build_populated_searcher(xb, IndexOptions(index_type="scann", gpu=True))
        self.assertIn("not supported", str(ctx.exception))


@unittest.skipUnless(SCANN_AVAILABLE, "scann is not installed")
class TestScannBuild(unittest.TestCase):
    def setUp(self) -> None:
        self.xb = unit_vectors(320, 32)

    def test_brute_force_self_hit(self) -> None:
        index, details = build_populated_searcher(
            self.xb, IndexOptions(index_type="scann", scann_num_neighbors=10)
        )
        self.assertEqual(details["index_type"], "ScaNN")
        self.assertEqual(details["scann_scoring"], "brute_force")
        self.assertEqual(details["metric"], "inner_product")
        self.assertEqual(index.ntotal, self.xb.shape[0])
        scores, labels = index.search(self.xb[:5], 1)
        np.testing.assert_array_equal(labels[:, 0], np.arange(5))
        self.assertTrue(np.all(scores[:, 0] > 0.99))
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "build" / "scann"
            serialize_index(index, src)
            dest = Path(tmp) / "search" / "scann"
            dest.parent.mkdir(parents=True)
            src.rename(dest)
            loaded = load_index(dest, ntotal=index.ntotal)
            scores2, labels2 = loaded.search(self.xb[:5], 1)
            np.testing.assert_array_equal(labels2[:, 0], np.arange(5))
            self.assertTrue(np.all(scores2[:, 0] > 0.99))

    def test_ah_tree_roundtrip(self) -> None:
        options = IndexOptions(
            index_type="scann",
            scann_leaves=8,
            scann_leaves_to_search=4,
            scann_reorder=20,
            scann_num_neighbors=10,
            scann_ah_dim=2,
        )
        index, details = build_populated_searcher(self.xb, options)
        self.assertEqual(details["scann_scoring"], "ah")
        self.assertTrue(details["scann_partitioned"])
        self.assertEqual(details["nlist"], 8)
        self.assertEqual(details["nprobe"], 4)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "build" / "scann"
            serialize_index(index, path)
            assets = (path / "scann_assets.pbtxt").read_text()
            self.assertIn("dataset.npy", assets)
            self.assertNotIn(str(path), assets)
            moved = Path(tmp) / "search" / "scann"
            moved.parent.mkdir(parents=True)
            path.rename(moved)
            loaded = load_index(
                moved,
                ntotal=index.ntotal,
                leaves_to_search=details["nprobe"],
                reorder=details["scann_reorder"],
            )
            apply_search_params(loaded, nprobe=4, reorder=20)
            scores, labels = loaded.search(self.xb[:3], 5)
            self.assertEqual(labels.shape, (3, 5))
            self.assertTrue(np.all(labels[:, 0] >= 0))
            self.assertEqual(scores.shape, (3, 5))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
