#!/usr/bin/env python3
"""Unit tests for node-level SQ / PQ / OPQ and SDC artifacts."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from node_quantization import (  # noqa: E402
    encode_pq,
    pack_pq_codes,
    reconstruct_pq,
    rotate_nodes,
    sdc_lookup_table,
    train_node_opq,
    train_node_pq,
    train_node_sq,
    unpack_pq_codes,
)


def l2_nodes(n: int, dim: int = 128, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    xb = rng.standard_normal((n, dim), dtype=np.float32)
    xb /= np.linalg.norm(xb, axis=1, keepdims=True).clip(min=1e-12)
    return np.ascontiguousarray(xb)


class TestNodeQuantization(unittest.TestCase):
    def test_pack_roundtrip(self) -> None:
        rng = np.random.default_rng(3)
        for pq_m, nbits in ((8, 4), (16, 4), (16, 8)):
            codes = rng.integers(0, 1 << nbits, size=(20, pq_m), dtype=np.uint8)
            packed = pack_pq_codes(codes, pq_m, nbits)
            np.testing.assert_array_equal(unpack_pq_codes(packed, pq_m, nbits), codes)

    def test_sdc_table_is_symmetric_inner_product(self) -> None:
        codebook = np.zeros((2, 4, 3), dtype=np.float32)
        codebook[0, 0] = [1, 0, 0]
        codebook[0, 1] = [0, 1, 0]
        table = sdc_lookup_table(codebook)
        self.assertEqual(table.shape, (2, 4, 4))
        self.assertAlmostEqual(float(table[0, 0, 0]), 1.0)
        self.assertAlmostEqual(float(table[0, 0, 1]), 0.0)
        np.testing.assert_allclose(table[0], table[0].T, atol=1e-6)

    def test_pq_encode_reconstructs(self) -> None:
        nodes = l2_nodes(64, 16, seed=4)
        codebook, codes = train_node_pq(nodes, pq_m=4, nbits=4, sample_size=64, niter=5, seed=1)
        rebuilt = reconstruct_pq(codes, codebook)
        self.assertEqual(codes.shape, (64, 4))
        self.assertEqual(rebuilt.shape, nodes.shape)
        self.assertLess(float(np.mean((nodes - rebuilt) ** 2)), 0.5)

    def test_opq_rotation_is_orthogonal(self) -> None:
        nodes = l2_nodes(80, 16, seed=5)
        rotation, codebook, codes = train_node_opq(
            nodes, pq_m=4, nbits=4, sample_size=80, niter=4, opq_iters=3, seed=2
        )
        gram = rotation @ rotation.T
        np.testing.assert_allclose(gram, np.eye(16), atol=1e-5)
        rotated = rotate_nodes(nodes, rotation)
        np.testing.assert_array_equal(encode_pq(rotated, codebook), codes)

    def test_sq_roundtrip_range(self) -> None:
        nodes = l2_nodes(40, 8, seed=6)
        scale, zero, codes = train_node_sq(nodes)
        self.assertEqual(codes.dtype, np.uint8)
        self.assertEqual(codes.shape, nodes.shape)
        self.assertTrue(np.all(scale > 0))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
