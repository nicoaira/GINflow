#!/usr/bin/env python3
"""Shared quantized-window encoding and HNSWLIB helpers."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def quantization_dir(path: Path) -> Path:
    if (path / "centroids.npy").exists():
        return path
    nested = path / "quantization"
    if (nested / "centroids.npy").exists():
        return nested
    raise FileNotFoundError(f"could not find quantization/centroids.npy below {path}")


def load_quantization(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    directory = quantization_dir(path)
    centroids_path = directory / "centroids.npy"
    similarity_path = directory / "similarity.npy"
    metadata_path = directory / "quantization.json"
    centroids = np.load(centroids_path)
    similarity = np.load(similarity_path)
    metadata = load_json(metadata_path) if metadata_path.exists() else {}
    if centroids.dtype != np.float16:
        raise ValueError(f"{centroids_path} must contain float16 centroids, got {centroids.dtype}")
    if centroids.ndim != 2 or centroids.shape[0] == 0 or centroids.shape[1] == 0:
        raise ValueError(f"invalid centroid shape {centroids.shape}")
    expected = (int(centroids.shape[0]), int(centroids.shape[0]))
    if similarity.shape != expected:
        raise ValueError(f"{similarity_path} shape {similarity.shape} does not match {expected}")
    if similarity.dtype != np.float32:
        raise ValueError(f"{similarity_path} must contain float32 similarities, got {similarity.dtype}")
    if int(metadata.get("n_centroids", centroids.shape[0])) != centroids.shape[0]:
        raise ValueError("quantization metadata n_centroids does not match centroids.npy")
    if int(metadata.get("embedding_dim", centroids.shape[1])) != centroids.shape[1]:
        raise ValueError("quantization metadata embedding_dim does not match centroids.npy")
    return (
        np.ascontiguousarray(centroids),
        np.ascontiguousarray(similarity),
        metadata,
    )


def encode_code_windows(codes: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    """Encode code windows so hnswlib IP equals the raw position-wise score."""
    values = np.asarray(codes)
    if values.ndim != 2:
        raise ValueError(f"code windows must be 2-D, got {values.shape}")
    if values.shape[0] == 0:
        return np.empty((0, int(values.shape[1]) * int(centroids.shape[1])), dtype=np.float32)
    codebook = np.ascontiguousarray(centroids, dtype=np.float32)
    if values.min() < 0 or values.max() >= codebook.shape[0]:
        raise ValueError(
            f"window code range [{int(values.min())}, {int(values.max())}] is outside "
            f"the centroid range [0, {codebook.shape[0]})"
        )
    encoded = codebook[np.asarray(values, dtype=np.int64)]
    flat = np.ascontiguousarray(encoded.reshape(values.shape[0], -1), dtype=np.float32)
    return flat


def raw_pairwise_score(
    query_codes: np.ndarray, target_codes: np.ndarray, similarity: np.ndarray
) -> np.ndarray:
    """Return the explicit sum of registered centroid similarities per window."""
    query = np.asarray(query_codes, dtype=np.int64)
    target = np.asarray(target_codes, dtype=np.int64)
    if query.ndim != 2 or target.ndim != 2 or query.shape != target.shape:
        raise ValueError(f"query/target code windows must have the same 2-D shape, got {query.shape} and {target.shape}")
    return np.sum(similarity[query, target], axis=1)


def similarity_from_distance(distances: np.ndarray) -> np.ndarray:
    """Convert hnswlib's IP distance to the raw registered score."""
    return np.asarray(1.0 - np.asarray(distances, dtype=np.float32), dtype=np.float32)


def require_hnswlib() -> Any:
    try:
        import hnswlib
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ValueError(
            "hnswlib is not installed. Use a Conda/Wave environment containing hnswlib==0.8.0."
        ) from exc
    return hnswlib


def create_index(
    dimension: int,
    max_elements: int,
    m: int,
    ef_construction: int,
    ef_search: int,
    random_seed: int,
    num_threads: int,
) -> Any:
    if dimension < 1 or max_elements < 1:
        raise ValueError("HNSWLIB dimension and max_elements must be >= 1")
    if m < 2:
        raise ValueError("--hnswlib-m must be >= 2")
    if ef_construction < 1 or ef_search < 1:
        raise ValueError("--hnswlib-ef-construction and --hnswlib-ef-search must be >= 1")
    hnswlib = require_hnswlib()
    index = hnswlib.Index(space="ip", dim=int(dimension))
    index.init_index(
        max_elements=int(max_elements),
        M=int(m),
        ef_construction=int(ef_construction),
        random_seed=int(random_seed),
    )
    if num_threads > 0:
        index.set_num_threads(int(num_threads))
    index.set_ef(int(ef_search))
    return index


def load_index(path: Path, dimension: int, ef_search: int, num_threads: int) -> Any:
    hnswlib = require_hnswlib()
    index = hnswlib.Index(space="ip", dim=int(dimension))
    index.load_index(str(path))
    if num_threads > 0:
        index.set_num_threads(int(num_threads))
    index.set_ef(int(ef_search))
    return index


def encoded_score_from_codes(
    query_codes: np.ndarray,
    target_codes: np.ndarray,
    centroids: np.ndarray,
    similarity: np.ndarray,
) -> np.ndarray:
    """Check the encoded inner product against the explicit custom score."""
    query_vectors = encode_code_windows(query_codes, centroids)
    target_vectors = encode_code_windows(target_codes, centroids)
    encoded = np.sum(query_vectors * target_vectors, axis=1)
    raw = raw_pairwise_score(query_codes, target_codes, similarity)
    return np.stack((encoded, raw), axis=1)
