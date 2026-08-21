#!/usr/bin/env python3
"""Build sliding windows from node-level SQ / PQ / OPQ codes."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from node_quantization import (
    decode_sq,
    load_json,
    n_windows,
    node_code_bytes,
    pack_pq_codes,
    normalize_rows,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-size", type=int, required=True)
    parser.add_argument("--stride", type=int, default=1)
    return parser.parse_args()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def pair_code_shards(input_dir: Path) -> list[tuple[Path, Path]]:
    nodes_dir = input_dir / "nodes"
    code_paths = sorted(nodes_dir.glob("*.quantized.npz"))
    manifest_paths = sorted(nodes_dir.glob("*.quantized.manifest.json"))
    code_map = {path.name[: -len(".quantized.npz")]: path for path in code_paths}
    manifest_map = {path.name[: -len(".quantized.manifest.json")]: path for path in manifest_paths}
    missing_codes = sorted(set(manifest_map) - set(code_map))
    missing_manifests = sorted(set(code_map) - set(manifest_map))
    if missing_codes or missing_manifests:
        raise ValueError(
            "quantized node / manifest shard names do not match: "
            f"only-codes={missing_manifests} only-manifests={missing_codes}"
        )
    if not code_map:
        raise ValueError(f"no quantized node shards found in {nodes_dir}")
    return [(code_map[key], manifest_map[key]) for key in sorted(code_map)]


def slice_pq_windows(codes: np.ndarray, window_size: int, stride: int, pq_m: int, nbits: int) -> np.ndarray:
    values = np.asarray(codes)
    if values.ndim != 2:
        raise ValueError(f"PQ node codes must be 2-D, got {values.shape}")
    count = n_windows(int(values.shape[0]), window_size, stride)
    if count == 0:
        width = window_size * node_code_bytes(pq_m, nbits)
        return np.empty((0, width), dtype=np.uint8)
    view = np.lib.stride_tricks.sliding_window_view(values, window_size, axis=0)[::stride][:count]
    stacked = np.ascontiguousarray(np.transpose(view, (0, 2, 1)))
    packed = pack_pq_codes(stacked.reshape(-1, pq_m), pq_m, nbits)
    return packed.reshape(count, window_size * packed.shape[1])


def slice_sq_windows(codes: np.ndarray, scale: np.ndarray, zero: np.ndarray, window_size: int, stride: int) -> np.ndarray:
    reconstructed = normalize_rows(decode_sq(codes, scale, zero))
    count = n_windows(int(reconstructed.shape[0]), window_size, stride)
    if count == 0:
        return np.empty((0, window_size * reconstructed.shape[1]), dtype=np.float32)
    view = np.lib.stride_tricks.sliding_window_view(reconstructed, window_size, axis=0)[::stride][:count]
    stacked = np.ascontiguousarray(np.transpose(view, (0, 2, 1)))
    return stacked.reshape(count, window_size * reconstructed.shape[1]).astype(np.float32, copy=False)


def generate_windows(input_dir: Path, output_dir: Path, window_size: int, stride: int) -> dict[str, Any]:
    if window_size < 1 or stride < 1:
        raise ValueError("window-size and stride must be >= 1")
    quantizer = load_json(input_dir / "quantizer.json")
    mode = str(quantizer["mode"])
    shards = pair_code_shards(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scale = zero = None
    if mode == "sq":
        scale = np.load(input_dir / "sq_scale.npy")
        zero = np.load(input_dir / "sq_zero.npy")
    pq_m = int(quantizer.get("pq_m", 0) or 0)
    nbits = int(quantizer.get("pq_nbits", 0) or 0)
    total_windows = 0
    records_written = 0
    for code_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        windows: dict[str, np.ndarray] = {}
        records: list[dict[str, Any]] = []
        skipped_short: list[dict[str, Any]] = []
        with np.load(code_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                length = int(record["length"])
                codes = np.asarray(arrays[identifier])
                if mode == "sq":
                    sliced = slice_sq_windows(codes, scale, zero, window_size, stride)
                else:
                    sliced = slice_pq_windows(codes, window_size, stride, pq_m, nbits)
                if sliced.shape[0] == 0:
                    skipped_short.append({"identifier": identifier, "length": length})
                    continue
                windows[identifier] = sliced
                records.append(
                    {
                        "identifier": identifier,
                        "length": length,
                        "n_windows": int(sliced.shape[0]),
                    }
                )
                total_windows += int(sliced.shape[0])
                records_written += 1
        prefix = code_path.name[: -len(".quantized.npz")]
        np.savez_compressed(output_dir / f"{prefix}.windows.npz", **windows)
        write_json(
            output_dir / f"{prefix}.windows.manifest.json",
            {
                "quantized_windows": True,
                "mode": mode,
                "window_size": window_size,
                "stride": stride,
                "window_dim": int(next(iter(windows.values())).shape[1]) if windows else 0,
                "pq_m": pq_m or None,
                "pq_nbits": nbits or None,
                "records": records,
                "skipped_short": skipped_short,
            },
        )
    summary = {
        "mode": mode,
        "window_size": window_size,
        "stride": stride,
        "records": records_written,
        "n_windows": total_windows,
    }
    write_json(output_dir / "windows.json", summary)
    return summary


def main() -> int:
    args = parse_args()
    summary = generate_windows(args.input_dir, args.output_dir, args.window_size, args.stride)
    print(
        f"wrote quantized windows: mode={summary['mode']} "
        f"records={summary['records']} n_windows={summary['n_windows']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
