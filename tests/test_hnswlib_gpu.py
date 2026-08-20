#!/usr/bin/env python3
"""Unit tests for the optional GPU HNSWLIB companion representation."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from hnswlib_gpu import (  # noqa: E402
    compatible_tuple,
    quantize_windows,
    pair_window_shards,
    read_window_rows,
    resolve_scale,
)


class TestHnswlibGpuRepresentation(unittest.TestCase):
    def test_auto_scale_uses_int8_range_without_clipping(self) -> None:
        scale = resolve_scale(0.14, None)
        self.assertAlmostEqual(scale, 850.0)
        values = np.asarray([[-0.14, 0.0, 0.14]], dtype=np.float32)
        encoded = quantize_windows(values, scale)
        self.assertEqual(encoded.dtype, np.int8)
        self.assertLessEqual(int(np.abs(encoded).max()), 127)

    def test_explicit_scale_rejects_clipping(self) -> None:
        with self.assertRaises(ValueError):
            resolve_scale(0.5, 300.0)

    def test_original_window_rows_keep_manifest_order_and_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            window_path = root / "shard.windows.npz"
            manifest_path = root / "shard.windows.manifest.json"
            values = np.arange(12, dtype=np.float32).reshape(3, 4)
            np.savez(window_path, record=values)
            manifest = {
                "window_size": 2,
                "stride": 1,
                "window_dim": 4,
                "checkpoint_sha256": "abc",
                "records": [{"identifier": "record", "n_windows": 3}],
            }
            manifest_path.write_text(json.dumps(manifest))
            loaded, mapping = read_window_rows(window_path, manifest)
            np.testing.assert_array_equal(loaded, values)
            self.assertEqual(mapping, [("record", 0, 2), ("record", 1, 3), ("record", 2, 4)])
            self.assertEqual(compatible_tuple(manifest), (2, 1, 4, "abc"))

    def test_window_manifests_can_be_staged_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            windows_dir = root / "windows"
            manifests_dir = root / "manifests"
            windows_dir.mkdir()
            manifests_dir.mkdir()
            np.savez(windows_dir / "shard.windows.npz", record=np.zeros((1, 4), dtype=np.float32))
            (manifests_dir / "shard.windows.manifest.json").write_text("{}")
            self.assertEqual(
                pair_window_shards(windows_dir, manifests_dir),
                [(windows_dir / "shard.windows.npz", manifests_dir / "shard.windows.manifest.json")],
            )


if __name__ == "__main__":
    unittest.main()
