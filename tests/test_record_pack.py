#!/usr/bin/env python3
"""Packed residue-embedding store: concatenated memmap plus legacy NPZ."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from record_pack import (  # noqa: E402
    load_embedding_files,
    load_residue_embeddings,
    write_residue_embeddings,
)


class TestRecordPack(unittest.TestCase):
    def test_roundtrip_and_selective_load(self) -> None:
        embeddings = {
            "a": np.arange(8, dtype=np.float16).reshape(2, 4),
            "b": np.arange(12, dtype=np.float16).reshape(3, 4) + 10,
            "c": np.arange(4, dtype=np.float16).reshape(1, 4) + 20,
        }
        with tempfile.TemporaryDirectory() as raw:
            outdir = Path(raw)
            write_residue_embeddings(outdir, embeddings)
            self.assertTrue((outdir / "embeddings.vectors.npy").is_file())
            self.assertTrue((outdir / "embeddings.npz").is_file())
            loaded = load_residue_embeddings(outdir)
            self.assertEqual(set(loaded), set(embeddings))
            np.testing.assert_array_equal(loaded["b"], embeddings["b"])
            subset = load_residue_embeddings(outdir, ids={"c", "a"})
            self.assertEqual(set(subset), {"a", "c"})
            np.testing.assert_array_equal(subset["a"], embeddings["a"])

    def test_legacy_npz_selective_and_sharded(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            first = directory / "one.npz"
            second = directory / "two.npz"
            np.savez(first, q1=np.ones((2, 3), dtype=np.float32), q2=np.zeros((2, 3), dtype=np.float32))
            np.savez(second, q3=np.full((2, 3), 2.0, dtype=np.float32))
            merged = load_embedding_files([first, second], ids={"q1", "q3"})
            self.assertEqual(set(merged), {"q1", "q3"})
            np.testing.assert_array_equal(merged["q3"], np.full((2, 3), 2.0, dtype=np.float32))
            with self.assertRaises(ValueError):
                load_embedding_files([first], ids={"q1", "missing"})


if __name__ == "__main__":
    raise SystemExit(unittest.main())
