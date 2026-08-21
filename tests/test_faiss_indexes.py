#!/usr/bin/env python3
"""CPU unit tests for FAISS index construction and score conversion."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import faiss  # noqa: E402
from faiss_index import (  # noqa: E402
    GPU_INDEX_TYPES,
    INDEX_TYPES,
    IndexOptions,
    build_populated_index,
    cuda_device_visible,
    distances_to_similarity,
    gpu_device_count,
    gpu_runtime_available,
    normalize_index_type,
    prepare_search_index,
    require_gpu,
    resolve_nlist,
    resolve_pq_m,
    supports_gpu,
)


def unit_vectors(n: int, dim: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    xb = rng.standard_normal((n, dim), dtype=np.float32)
    xb /= np.linalg.norm(xb, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(xb)


class TestHelpers(unittest.TestCase):
    def test_aliases(self) -> None:
        expected = {
            "flatip": "FlatIP",
            "flatl2": "FlatL2",
            "hnsw": "HNSW",
            "ivfflat": "IVFFlat",
        }
        for public_value, canonical in expected.items():
            with self.subTest(public_value=public_value):
                self.assertEqual(normalize_index_type(public_value), canonical)
        self.assertEqual(normalize_index_type("IndexIVFFlat"), "IVFFlat")
        self.assertEqual(normalize_index_type("hnsw_flat"), "HNSW")
        with self.assertRaises(ValueError):
            normalize_index_type("CAGRA")

    def test_gpu_matrix(self) -> None:
        for kind in ("flatip", "flatl2"):
            self.assertTrue(supports_gpu(kind))
        for kind in ("hnsw", "ivfflat"):
            self.assertFalse(supports_gpu(kind))
        self.assertTrue(GPU_INDEX_TYPES.issubset(set(INDEX_TYPES)))

    def test_require_gpu_rejects_hnsw(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            require_gpu("hnsw")
        self.assertIn("not supported", str(ctx.exception))

    def test_nlist_auto_and_bounds(self) -> None:
        self.assertEqual(resolve_nlist(10, 4), 4)
        self.assertGreaterEqual(resolve_nlist(1000, None), 1)
        self.assertLessEqual(resolve_nlist(1000, None), 1000)
        with self.assertRaises(ValueError):
            resolve_nlist(10, 11)

    def test_pq_m_must_divide_dim(self) -> None:
        self.assertEqual(resolve_pq_m(1408, 16), 16)
        with self.assertRaises(ValueError) as ctx:
            resolve_pq_m(1408, 10)
        self.assertIn("1408", str(ctx.exception))

    def test_l2_and_hamming_conversion(self) -> None:
        # Unit vectors: ||x-x||^2 = 0 -> cosine 1; ||x-y||^2 = 2 -> cosine 0
        l2 = np.array([[0.0, 2.0]], dtype=np.float32)
        sim = distances_to_similarity(l2, "l2")
        np.testing.assert_allclose(sim, [[1.0, 0.0]], atol=1e-6)
        _ = distances_to_similarity


class TestBuildIndexes(unittest.TestCase):
    def setUp(self) -> None:
        self.xb = unit_vectors(320, 32)

    def _roundtrip(self, options: IndexOptions) -> None:
        index, details = build_populated_index(self.xb, options)
        self.assertEqual(index.ntotal, self.xb.shape[0])
        self.assertEqual(details["index_type"], normalize_index_type(options.index_type))
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / "index.faiss")
            faiss.write_index(index, path)
            loaded = faiss.read_index(path)
            self.assertEqual(loaded.ntotal, index.ntotal)
            scores, labels = loaded.search(self.xb[:3], 5)
            self.assertEqual(labels.shape, (3, 5))
            self.assertTrue(np.all(labels[:, 0] >= 0))

    def test_flatip_self_hit(self) -> None:
        index, details = build_populated_index(self.xb, IndexOptions(index_type="flatip"))
        self.assertEqual(details["metric"], "inner_product")
        scores, labels = index.search(self.xb[:5], 1)
        np.testing.assert_array_equal(labels[:, 0], np.arange(5))
        self.assertTrue(np.all(scores[:, 0] > 0.99))

    def test_flatl2(self) -> None:
        self._roundtrip(IndexOptions(index_type="flatl2"))

    def test_hnsw(self) -> None:
        self._roundtrip(IndexOptions(index_type="hnsw", hnsw_m=8, hnsw_ef_search=16))

    def test_ivfflat(self) -> None:
        self._roundtrip(IndexOptions(index_type="ivfflat", nlist=4, nprobe=4))

    def test_unknown_pq_type_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            normalize_index_type("ivfpq")

    def test_gpu_without_runtime(self) -> None:
        if gpu_runtime_available() and gpu_device_count() > 0:
            self.skipTest("GPU FAISS is available in this environment")
        with self.assertRaises(ValueError) as ctx:
            build_populated_index(self.xb, IndexOptions(index_type="flatip", gpu=True))
        message = str(ctx.exception)
        self.assertTrue("faiss-cpu" in message or "CUDA" in message or "gpu" in message.lower())


def gpu_ready() -> bool:
    return gpu_runtime_available() and cuda_device_visible() and gpu_device_count() > 0


@unittest.skipUnless(gpu_ready(), "FAISS GPU runtime and a CUDA device are required")
class TestGpuIndexes(unittest.TestCase):
    def setUp(self) -> None:
        self.xb = unit_vectors(320, 32)

    def _roundtrip_gpu(self, options: IndexOptions) -> None:
        index, details = build_populated_index(self.xb, options)
        self.assertTrue(details["gpu"])
        self.assertEqual(index.ntotal, self.xb.shape[0])
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / "index.faiss")
            faiss.write_index(index, path)
            loaded = faiss.read_index(path)
            search_opts = IndexOptions(gpu=True, nprobe=options.nprobe, hnsw_ef_search=options.hnsw_ef_search)
            loaded, metric, lsh_nbits = prepare_search_index(
                loaded, {"index_type": details["index_type"], "metric": details["metric"]}, search_opts
            )
            scores, labels = loaded.search(self.xb[:3], 5)
            self.assertEqual(labels.shape, (3, 5))
            self.assertTrue(np.all(labels[:, 0] >= 0))
            _ = metric, lsh_nbits

    def test_gpu_flatip(self) -> None:
        self._roundtrip_gpu(IndexOptions(index_type="flatip", gpu=True))

    def test_gpu_flatl2(self) -> None:
        self._roundtrip_gpu(IndexOptions(index_type="flatl2", gpu=True))

    def test_gpu_ivfflat_is_rejected(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            build_populated_index(
                self.xb, IndexOptions(index_type="ivfflat", nlist=4, nprobe=4, gpu=True)
            )
        self.assertIn("--index ivf", str(ctx.exception))

    def test_gpu_rejected_for_cpu_only_types(self) -> None:
        for kind in sorted(set(INDEX_TYPES) - set(GPU_INDEX_TYPES)):
            with self.subTest(kind=kind):
                with self.assertRaises(ValueError) as ctx:
                    build_populated_index(self.xb, IndexOptions(index_type=kind, gpu=True, nlist=4, pq_m=8))
                self.assertIn("not supported", str(ctx.exception))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
