#!/usr/bin/env python3
"""Unit tests for node-level PQ CAGRA with shared-memory ADC."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor" / "cagra-pq-adc" / "python"))

try:
    import pq_cagra_adc as pq
except ImportError as exc:  # pragma: no cover - extension not built
    pq = None
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def require_ext() -> None:
    if pq is None:
        raise unittest.SkipTest(f"pq_cagra_adc extension is unavailable: {_IMPORT_ERROR}")


class TestPqCagraAdc(unittest.TestCase):
    def test_window_slicing_is_node_major(self) -> None:
        require_ext()
        nodes = np.arange(20, dtype=np.uint8).reshape(10, 2)
        windows = pq.slice_code_windows(nodes, [5, 5], window_size=3, stride=1)
        self.assertEqual(windows.shape, (6, 3, 2))
        np.testing.assert_array_equal(windows[0], nodes[0:3])
        np.testing.assert_array_equal(windows[3], nodes[5:8])

    def test_opq_rotation_is_orthogonal_and_reduces_error(self) -> None:
        require_ext()
        rng = np.random.default_rng(5)
        mix = rng.normal(size=(32, 32)).astype(np.float32)
        nodes = pq.normalize_rows(rng.normal(size=(400, 32)).astype(np.float32) @ mix)
        codebook_pq, codes_pq = pq.train_node_pq(nodes, pq_m=8, nbits=4, sample_size=400, niter=6, seed=5)
        rotation, codebook_opq, codes_opq = pq.train_node_opq(
            nodes, pq_m=8, nbits=4, sample_size=400, niter=6, opq_iters=5, seed=5
        )
        np.testing.assert_allclose(rotation @ rotation.T, np.eye(32), atol=1e-5)
        self.assertAlmostEqual(float(np.linalg.det(rotation)), 1.0, places=4)
        mse_pq = pq.reconstruction_mse(nodes, codes_pq, codebook_pq)
        mse_opq = pq.reconstruction_mse(pq.rotate_nodes(nodes, rotation), codes_opq, codebook_opq)
        self.assertLess(mse_opq, mse_pq)
        rotated_windows = pq.rotate_windows(nodes[:6].reshape(2, 3, 32), rotation)
        self.assertEqual(rotated_windows.shape, (2, 3, 32))

    def test_adc_matches_reconstructed_inner_product(self) -> None:
        require_ext()
        rng = np.random.default_rng(3)
        pq_m, ksub, dsub, window_size = 4, 16, 8, 3
        codebook = rng.normal(size=(pq_m, ksub, dsub)).astype(np.float32)
        queries = pq.normalize_rows(rng.normal(size=(5, window_size * pq_m * dsub))).reshape(5, window_size, pq_m * dsub)
        codes = rng.integers(0, ksub, size=(7, window_size, pq_m), dtype=np.uint8)
        reconstructed = codebook[np.arange(pq_m)[None, None, :], codes, :]
        reconstructed = reconstructed.reshape(7, window_size, pq_m * dsub)
        expected = -np.einsum("qwd,nwd->qn", queries, reconstructed)
        actual = pq.adc_distances(queries, codebook, codes)
        np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-5)

    def test_gpu_search_matches_exact_adc_on_tiny_complete_graph(self) -> None:
        require_ext()
        rng = np.random.default_rng(11)
        n, window_size, pq_m, nbits, dsub = 48, 3, 4, 4, 8
        nodes = rng.normal(size=(n + 20, pq_m * dsub)).astype(np.float32)
        codebook, node_codes = pq.train_node_pq(nodes, pq_m, nbits, sample_size=n + 20, niter=8, seed=11)
        lengths = [n + 20]
        codes = pq.slice_code_windows(node_codes, lengths, window_size, 1)[:n]
        queries = pq.slice_float_windows(nodes, lengths, window_size, 1)[:8]
        index = pq.build_index(codes, codebook, graph_degree=n - 1, intermediate_graph_degree=n - 1)
        labels, distances = pq.search(index, queries, k=5, itopk_size=16, search_width=1, min_iterations=2, max_iterations=8)
        exact = pq.adc_distances(queries, codebook, codes)
        exact_labels = np.argsort(exact, axis=1)[:, :5]
        np.testing.assert_array_equal(labels, exact_labels)
        np.testing.assert_allclose(distances, np.take_along_axis(exact, labels, axis=1), rtol=1e-3, atol=1e-3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "index.bin"
            pq.save_index(index, path)
            loaded = pq.load_index(path)
        labels2, _distances2 = pq.search(loaded, queries, k=5, itopk_size=16, max_iterations=8)
        np.testing.assert_array_equal(labels2, labels)

    def test_cpu_search_matches_exact_adc_on_tiny_complete_graph(self) -> None:
        require_ext()
        rng = np.random.default_rng(11)
        n, window_size, pq_m, nbits, dsub = 48, 3, 4, 4, 8
        nodes = rng.normal(size=(n + 20, pq_m * dsub)).astype(np.float32)
        codebook, node_codes = pq.train_node_pq(nodes, pq_m, nbits, sample_size=n + 20, niter=8, seed=11)
        lengths = [n + 20]
        codes = pq.slice_code_windows(node_codes, lengths, window_size, 1)[:n]
        queries = pq.slice_float_windows(nodes, lengths, window_size, 1)[:8]
        index = pq.build_index(codes, codebook, graph_degree=n - 1, intermediate_graph_degree=n - 1)
        labels, distances = pq.search(
            index, queries, k=5, itopk_size=16, search_width=1, min_iterations=2, max_iterations=8, device="cpu"
        )
        exact = pq.adc_distances(queries, codebook, codes)
        exact_labels = np.argsort(exact, axis=1)[:, :5]
        np.testing.assert_array_equal(labels, exact_labels)
        np.testing.assert_allclose(distances, np.take_along_axis(exact, labels, axis=1), rtol=1e-5, atol=1e-5)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "index.bin"
            pq.save_index(index, path)
            loaded = pq.load_index(path)
        labels2, _ = pq.search(
            loaded, queries, k=5, itopk_size=16, max_iterations=8, device="cpu", num_threads=2
        )
        np.testing.assert_array_equal(labels2, labels)


if __name__ == "__main__":
    unittest.main()
