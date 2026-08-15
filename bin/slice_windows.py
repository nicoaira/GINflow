#!/usr/bin/env python3
"""Slice per-nucleotide embeddings into concatenated, L2-normalized windows."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def n_windows(length: int, window_size: int, stride: int) -> int:
    if length < window_size:
        return 0
    return (length - window_size) // stride + 1


def slice_record(emb: np.ndarray, window_size: int, stride: int) -> np.ndarray:
    """Return (n_windows, window_size * dim) unit-normalized float32 windows."""
    if emb.ndim != 2:
        raise ValueError(f"expected a 2-D embedding array, got shape {emb.shape}")
    length, dim = emb.shape
    count = n_windows(length, window_size, stride)
    out_dim = window_size * dim
    if count == 0:
        return np.empty((0, out_dim), dtype=np.float32)
    views = np.lib.stride_tricks.sliding_window_view(emb, window_size, axis=0)[::stride]
    flat = np.ascontiguousarray(views.transpose(0, 2, 1).reshape(count, out_dim), dtype=np.float32)
    norms = np.linalg.norm(flat, axis=1, keepdims=True)
    norms = np.maximum(norms, np.float32(1e-12))
    return flat / norms


def load_embed_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text())
    if not isinstance(manifest, dict):
        raise ValueError(f"{path} is not a JSON object")
    return manifest


def slice_npz(npz_path: Path, manifest: dict, window_size: int, stride: int) -> tuple[dict[str, np.ndarray], dict]:
    arrays = np.load(npz_path)
    records = []
    skipped_short = []
    windows: dict[str, np.ndarray] = {}
    embedding_dim = None

    for record in manifest.get("records", []):
        identifier = record["identifier"]
        if identifier not in arrays.files:
            raise KeyError(f"{identifier} is in the embed manifest but missing from {npz_path}")
        emb = np.asarray(arrays[identifier])
        if embedding_dim is None:
            embedding_dim = int(emb.shape[1])
        elif int(emb.shape[1]) != embedding_dim:
            raise ValueError(
                f"{identifier} has embedding dim {emb.shape[1]}, expected {embedding_dim}"
            )
        sliced = slice_record(emb, window_size, stride)
        length = int(record.get("core_length", emb.shape[0]))
        if sliced.shape[0] == 0:
            skipped_short.append({"identifier": identifier, "length": length})
            continue
        windows[identifier] = sliced
        records.append({
            "identifier": identifier,
            "length": length,
            "n_windows": int(sliced.shape[0]),
            "shape": [int(sliced.shape[0]), int(sliced.shape[1])],
        })

    if embedding_dim is None:
        raise ValueError(f"{npz_path} produced no embedding records")

    window_manifest = {
        "status": "complete",
        "window_size": window_size,
        "stride": stride,
        "embedding_dim": embedding_dim,
        "window_dim": window_size * embedding_dim,
        "l2_normalized": True,
        "ginfinity_version": manifest.get("ginfinity_version"),
        "model_version": manifest.get("model_version"),
        "checkpoint_sha256": manifest.get("checkpoint_sha256"),
        "input": str(npz_path),
        "input_manifest": str(manifest.get("output", "")),
        "records": records,
        "skipped_short": skipped_short,
    }
    return windows, window_manifest


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Embedding NPZ from ginfinity embed-graphs")
    parser.add_argument("--manifest", type=Path, required=True, help="Embedding manifest JSON")
    parser.add_argument("--output", type=Path, required=True, help="Output windows NPZ")
    parser.add_argument("--windows-manifest", type=Path, required=True, help="Output windows manifest JSON")
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.window_size < 1:
        print("error: --window-size must be >= 1", file=sys.stderr)
        return 2
    if args.stride < 1:
        print("error: --stride must be >= 1", file=sys.stderr)
        return 2

    manifest = load_embed_manifest(args.manifest)
    windows, window_manifest = slice_npz(args.input, manifest, args.window_size, args.stride)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **windows)
    args.windows_manifest.write_text(json.dumps(window_manifest, indent=2) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "windows_manifest": str(args.windows_manifest),
        "records": len(window_manifest["records"]),
        "skipped_short": len(window_manifest["skipped_short"]),
        "n_windows": int(sum(item["n_windows"] for item in window_manifest["records"])),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
