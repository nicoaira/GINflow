#!/usr/bin/env python3
"""Generate sliding windows of centroid codes for the HNSWLIB candidate index."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def n_windows(length: int, window_size: int, stride: int) -> int:
    if length < window_size:
        return 0
    return (length - window_size) // stride + 1


def slice_codes(codes: np.ndarray, window_size: int, stride: int) -> np.ndarray:
    values = np.asarray(codes)
    if values.ndim != 1:
        raise ValueError(f"quantized node codes must be 1-D, got {values.shape}")
    count = n_windows(int(values.shape[0]), window_size, stride)
    if count == 0:
        return np.empty((0, window_size), dtype=values.dtype)
    windows = np.lib.stride_tricks.sliding_window_view(values, window_size)[::stride]
    return np.ascontiguousarray(windows[:count])


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


def generate_windows(input_dir: Path, output_dir: Path, window_size: int, stride: int) -> dict:
    if window_size < 1:
        raise ValueError("--window-size must be >= 1")
    if stride < 1:
        raise ValueError("--stride must be >= 1")
    shards = pair_code_shards(input_dir)
    quantization = load_json(input_dir / "quantization.json")
    output_dir.mkdir(parents=True, exist_ok=True)
    records_written = 0
    total_windows = 0
    for code_path, code_manifest_path in shards:
        code_manifest = load_json(code_manifest_path)
        windows: dict[str, np.ndarray] = {}
        records = []
        skipped_short = []
        with np.load(code_path) as arrays:
            for record in code_manifest.get("records", []):
                identifier = str(record["identifier"])
                if identifier not in arrays.files:
                    raise KeyError(f"{identifier} is in {code_manifest_path} but missing from {code_path}")
                codes = np.asarray(arrays[identifier])
                sliced = slice_codes(codes, window_size, stride)
                length = int(record.get("length", codes.shape[0]))
                if sliced.shape[0] == 0:
                    skipped_short.append({"identifier": identifier, "length": length})
                    continue
                windows[identifier] = sliced
                records.append(
                    {
                        "identifier": identifier,
                        "length": length,
                        "n_windows": int(sliced.shape[0]),
                        "shape": [int(sliced.shape[0]), int(sliced.shape[1])],
                    }
                )
                total_windows += int(sliced.shape[0])
                records_written += 1

        prefix = code_path.name[: -len(".quantized.npz")]
        output_npz = output_dir / f"{prefix}.windows.npz"
        output_manifest = output_dir / f"{prefix}.windows.manifest.json"
        np.savez_compressed(output_npz, **windows)
        output_manifest.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "quantized_nodes": True,
                    "quantized_windows": True,
                    "metric": quantization.get("metric", "cosine"),
                    "normalization": quantization.get("normalization", "per_node_l2"),
                    "n_centroids": int(quantization["n_centroids"]),
                    "embedding_dim": int(quantization["embedding_dim"]),
                    "centroid_dtype": quantization.get("centroid_dtype", "float16"),
                    "similarity_file": quantization.get("similarity_file", "similarity.npy"),
                    "window_size": int(window_size),
                    "stride": int(stride),
                    "code_dim": int(window_size),
                    "window_dim": int(window_size) * int(quantization["embedding_dim"]),
                    "l2_normalized": False,
                    "input": code_path.name,
                    "input_manifest": code_manifest_path.name,
                    "ginfinity_version": code_manifest.get("ginfinity_version"),
                    "model_version": code_manifest.get("model_version"),
                    "checkpoint_sha256": code_manifest.get("checkpoint_sha256"),
                    "records": records,
                    "skipped_short": skipped_short,
                },
                indent=2,
            )
            + "\n"
        )

    result = {
        "status": "complete",
        "n_shards": len(shards),
        "n_records": records_written,
        "n_windows": total_windows,
        "window_size": window_size,
        "stride": stride,
        "n_centroids": int(quantization["n_centroids"]),
        "embedding_dim": int(quantization["embedding_dim"]),
    }
    (output_dir / "windows.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = generate_windows(args.input_dir, args.output_dir, args.window_size, args.stride)
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"outdir": str(args.output_dir), **result}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
