#!/usr/bin/env python3
"""I/O helpers for the compact C++ hnswlib backend."""
from __future__ import annotations

import csv
import shutil
from pathlib import Path

import numpy as np

from build_hnswlib import load_json, pair_window_shards, read_code_rows


def pack_code_windows(
    windows_dir: Path,
    codes_path: Path,
) -> tuple[int, int, list[tuple[str, int, int]], list[dict]]:
    """Pack quantized window shards into the raw uint16 C++ input format."""
    shards = pair_window_shards(windows_dir)
    manifests = [load_json(manifest_path) for _window_path, manifest_path in shards]
    reference = manifests[0]
    if not reference.get("quantized_windows"):
        raise ValueError(f"{windows_dir} does not contain quantized window manifests")
    expected = (
        reference.get("window_size"),
        reference.get("stride"),
        reference.get("window_dim"),
        reference.get("n_centroids"),
        reference.get("embedding_dim"),
    )
    for manifest in manifests[1:]:
        actual = (
            manifest.get("window_size"),
            manifest.get("stride"),
            manifest.get("window_dim"),
            manifest.get("n_centroids"),
            manifest.get("embedding_dim"),
        )
        if actual != expected:
            raise ValueError(f"incompatible quantized window manifest: {actual} vs {expected}")

    all_codes: list[np.ndarray] = []
    mapping: list[tuple[str, int, int]] = []
    window_size = int(reference["window_size"])
    n_centroids = int(reference["n_centroids"])
    for window_path, manifest_path in shards:
        codes, local_mapping = read_code_rows(window_path, load_json(manifest_path))
        values = np.asarray(codes)
        if values.dtype.kind not in "ui" or values.ndim != 2 or values.shape[1] != window_size:
            raise ValueError(f"{window_path} does not contain an integer code matrix")
        if n_centroids > np.iinfo(np.uint16).max + 1:
            raise ValueError("compact hnswlib currently supports at most 65536 centroids")
        if values.size and (int(values.min()) < 0 or int(values.max()) >= n_centroids):
            raise ValueError(f"codes in {window_path} are outside [0, {n_centroids})")
        all_codes.append(np.ascontiguousarray(values, dtype=np.uint16))
        mapping.extend(local_mapping)

    if not all_codes:
        raise ValueError(f"{windows_dir} contains no quantized windows")
    packed = np.ascontiguousarray(np.concatenate(all_codes, axis=0), dtype=np.uint16)
    codes_path.parent.mkdir(parents=True, exist_ok=True)
    packed.tofile(codes_path)
    return int(packed.shape[0]), window_size, mapping, manifests


def write_mapping(path: Path, mapping: list[tuple[str, int, int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        writer.writerows((row_id, identifier, start, end) for row_id, (identifier, start, end) in enumerate(mapping))


def write_similarity_binary(quantization_dir: Path, path: Path) -> tuple[np.ndarray, Path]:
    source = quantization_dir / "similarity.npy"
    values = np.asarray(np.load(source))
    if values.ndim != 2 or values.shape[0] != values.shape[1] or values.dtype != np.float32:
        raise ValueError(f"{source} must be a square float32 matrix, got {values.shape} {values.dtype}")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.ascontiguousarray(values, dtype=np.float32).tofile(path)
    return values, source


def copy_quantization(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("centroids.npy", "similarity.npy", "quantization.json"):
        source_file = source / name
        if not source_file.exists():
            raise FileNotFoundError(f"missing quantization artifact {source_file}")
        shutil.copy2(source_file, destination / name)


def read_compact_results(
    labels_path: Path,
    distances_path: Path,
    n_queries: int,
    requested_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.fromfile(labels_path, dtype="<u8")
    distances = np.fromfile(distances_path, dtype="<f4")
    expected = n_queries * requested_k
    if labels.size != expected or distances.size != expected:
        raise ValueError(
            f"compact search returned {labels.size} labels/{distances.size} distances; "
            f"expected {expected} of each"
        )
    return labels.reshape(n_queries, requested_k), distances.reshape(n_queries, requested_k)
