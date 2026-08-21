#!/usr/bin/env python3
"""Node-level SQ / PQ / OPQ for GINflow windows.

Compression is applied to 128-d residue embeddings. Sliding windows are then
formed from the compressed nodes. Product quantization of concatenated 1408-d
windows is not supported.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterator

import numpy as np


NODE_DIM = 128


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def shard_prefix(path: Path, suffix: str) -> str:
    if not path.name.endswith(suffix):
        raise ValueError(f"{path} does not end with {suffix}")
    return path.name[: -len(suffix)]


def pair_shards(embeddings: list[Path], manifests: list[Path]) -> list[tuple[Path, Path]]:
    embedding_map = {shard_prefix(path, ".npz"): path for path in embeddings}
    manifest_map = {shard_prefix(path, ".manifest.json"): path for path in manifests}
    missing_embeddings = sorted(set(manifest_map) - set(embedding_map))
    missing_manifests = sorted(set(embedding_map) - set(manifest_map))
    if missing_embeddings or missing_manifests:
        raise ValueError(
            "embedding / manifest shard names do not match: "
            f"only-embeddings={missing_manifests} only-manifests={missing_embeddings}"
        )
    return [(embedding_map[key], manifest_map[key]) for key in sorted(embedding_map)]


def normalize_rows(values: np.ndarray) -> np.ndarray:
    rows = np.ascontiguousarray(values, dtype=np.float32)
    if rows.ndim != 2:
        raise ValueError(f"expected a 2-D matrix, got {rows.shape}")
    norms = np.linalg.norm(rows, axis=1, keepdims=True)
    return rows / np.maximum(norms, np.float32(1e-12))


def n_windows(length: int, window_size: int, stride: int) -> int:
    if length < window_size:
        return 0
    return (length - window_size) // stride + 1


def node_code_bytes(pq_m: int, nbits: int) -> int:
    return (pq_m * nbits + 7) // 8


def pack_pq_codes(codes: np.ndarray, pq_m: int, nbits: int) -> np.ndarray:
    """Pack one-byte subcodes into little-endian bytes (FAISS-compatible)."""
    values = np.asarray(codes, dtype=np.uint8)
    if values.ndim != 2 or values.shape[1] != pq_m:
        raise ValueError(f"PQ codes must have shape (n, {pq_m}), got {values.shape}")
    if values.size and int(values.max()) >= (1 << nbits):
        raise ValueError("PQ code exceeds the requested subquantizer range")
    bits = ((values[:, :, None] >> np.arange(nbits, dtype=np.uint8)) & 1).reshape(
        values.shape[0], pq_m * nbits
    )
    return np.packbits(bits, axis=1, bitorder="little")


def unpack_pq_codes(packed: np.ndarray, pq_m: int, nbits: int) -> np.ndarray:
    values = np.asarray(packed, dtype=np.uint8)
    if values.ndim != 2:
        raise ValueError(f"packed PQ codes must be 2-D, got {values.shape}")
    bits = np.unpackbits(values, axis=1, bitorder="little")
    needed = pq_m * nbits
    if bits.shape[1] < needed:
        raise ValueError("packed PQ code buffer is too short")
    bit_rows = bits[:, :needed].reshape(values.shape[0], pq_m, nbits)
    weights = (np.uint16(1) << np.arange(nbits, dtype=np.uint16))[None, None, :]
    return np.asarray((bit_rows * weights).sum(axis=2), dtype=np.uint8)


def sdc_lookup_table(codebook: np.ndarray) -> np.ndarray:
    """Symmetric inner-product table, shape (M, ksub, ksub)."""
    book = np.ascontiguousarray(codebook, dtype=np.float32)
    if book.ndim != 3:
        raise ValueError(f"codebook must be (M, ksub, dsub), got {book.shape}")
    return np.einsum("m i d, m j d -> m i j", book, book, optimize=True).astype(np.float32, copy=False)


def _kmeans(vectors: np.ndarray, k: int, niter: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    count = int(vectors.shape[0])
    k = min(k, count)
    centroids = np.array(vectors[rng.choice(count, size=k, replace=False)], dtype=np.float32, copy=True)
    for _ in range(niter):
        dots = vectors @ centroids.T
        sq_v = np.einsum("nd,nd->n", vectors, vectors)
        sq_c = np.einsum("kd,kd->k", centroids, centroids)
        labels = np.argmin(sq_v[:, None] + sq_c[None, :] - 2.0 * dots, axis=1)
        sums = np.zeros_like(centroids)
        np.add.at(sums, labels, vectors)
        counts = np.bincount(labels, minlength=k)
        nonempty = counts > 0
        centroids[nonempty] = sums[nonempty] / counts[nonempty, None]
        if np.any(~nonempty):
            centroids[~nonempty] = vectors[rng.choice(count, size=int((~nonempty).sum()), replace=False)]
    return centroids


def _validate_pq(dimension: int, pq_m: int, nbits: int) -> tuple[int, int]:
    if pq_m < 1 or dimension % pq_m != 0:
        raise ValueError(f"embedding dimension {dimension} must be divisible by pq_m={pq_m}")
    if nbits < 1 or nbits > 8:
        raise ValueError("nbits must be in [1, 8]")
    return dimension // pq_m, 1 << nbits


def _sample_nodes(values: np.ndarray, sample_size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sample_count = min(int(sample_size), int(values.shape[0]))
    if sample_count == values.shape[0]:
        return values
    return np.ascontiguousarray(values[rng.choice(values.shape[0], size=sample_count, replace=False)])


def _fit_pq_codebook(sample: np.ndarray, pq_m: int, ksub: int, dsub: int, niter: int, seed: int) -> np.ndarray:
    codebook = np.empty((pq_m, ksub, dsub), dtype=np.float32)
    for sub in range(pq_m):
        subvectors = np.ascontiguousarray(sample[:, sub * dsub : (sub + 1) * dsub], dtype=np.float32)
        codebook[sub] = _kmeans(subvectors, ksub, niter, seed + sub)
    return codebook


def encode_pq(values: np.ndarray, codebook: np.ndarray) -> np.ndarray:
    book = np.ascontiguousarray(codebook, dtype=np.float32)
    rows = np.ascontiguousarray(values, dtype=np.float32)
    if rows.ndim != 2 or book.ndim != 3 or rows.shape[1] != book.shape[0] * book.shape[2]:
        raise ValueError(f"values shape {rows.shape} does not match codebook {book.shape}")
    pq_m, ksub, dsub = book.shape
    codes = np.empty((rows.shape[0], pq_m), dtype=np.uint8)
    for sub in range(pq_m):
        centroids = book[sub]
        all_sub = rows[:, sub * dsub : (sub + 1) * dsub]
        sq_c = np.einsum("kd,kd->k", centroids, centroids)
        assigned = np.empty(all_sub.shape[0], dtype=np.int64)
        for start in range(0, all_sub.shape[0], 65536):
            stop = min(start + 65536, all_sub.shape[0])
            chunk = all_sub[start:stop]
            dots = chunk @ centroids.T
            sq_v = np.einsum("nd,nd->n", chunk, chunk)
            assigned[start:stop] = np.argmin(sq_v[:, None] + sq_c[None, :] - 2.0 * dots, axis=1)
        codes[:, sub] = assigned.astype(np.uint8)
    return codes


def reconstruct_pq(codes: np.ndarray, codebook: np.ndarray) -> np.ndarray:
    book = np.ascontiguousarray(codebook, dtype=np.float32)
    packed = np.ascontiguousarray(codes)
    if packed.ndim != 2 or book.ndim != 3 or packed.shape[1] != book.shape[0]:
        raise ValueError(f"codes shape {packed.shape} does not match codebook {book.shape}")
    pq_m, _ksub, dsub = book.shape
    reconstructed = np.empty((packed.shape[0], pq_m * dsub), dtype=np.float32)
    for sub in range(pq_m):
        reconstructed[:, sub * dsub : (sub + 1) * dsub] = book[sub, packed[:, sub]]
    return reconstructed


def rotate_nodes(nodes: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    values = np.ascontiguousarray(nodes, dtype=np.float32)
    matrix = np.ascontiguousarray(rotation, dtype=np.float32)
    if values.ndim != 2 or matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("rotation must be a square matrix applied to 2-D nodes")
    if values.shape[1] != matrix.shape[0]:
        raise ValueError(f"node dim {values.shape[1]} does not match rotation {matrix.shape}")
    return np.ascontiguousarray(values @ matrix, dtype=np.float32)


def rotate_windows(windows: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    matrix = np.ascontiguousarray(rotation, dtype=np.float32)
    values = np.ascontiguousarray(windows, dtype=np.float32)
    dim = int(matrix.shape[0])
    if values.ndim == 2:
        if values.shape[1] % dim != 0:
            raise ValueError(f"flat window width {values.shape[1]} is not a multiple of {dim}")
        rotated = rotate_nodes(values.reshape(-1, dim), matrix)
        return rotated.reshape(values.shape)
    if values.ndim == 3:
        if values.shape[2] != dim:
            raise ValueError(f"window dim {values.shape[2]} does not match rotation {matrix.shape}")
        rotated = rotate_nodes(values.reshape(-1, dim), matrix)
        return rotated.reshape(values.shape)
    raise ValueError("windows must be 2-D or 3-D")


def _orthogonal_procrustes(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    gram = source.T @ target
    u_mat, _singular, vt_mat = np.linalg.svd(gram, full_matrices=True)
    rotation = u_mat @ vt_mat
    if np.linalg.det(rotation) < 0:
        u_mat[:, -1] *= -1
        rotation = u_mat @ vt_mat
    return np.ascontiguousarray(rotation, dtype=np.float32)


def train_node_pq(
    nodes: np.ndarray,
    pq_m: int,
    nbits: int,
    sample_size: int = 250_000,
    niter: int = 12,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    values = normalize_rows(nodes)
    dsub, ksub = _validate_pq(int(values.shape[1]), pq_m, nbits)
    sample = _sample_nodes(values, sample_size, seed)
    codebook = _fit_pq_codebook(sample, pq_m, ksub, dsub, niter, seed)
    return codebook, encode_pq(values, codebook)


def train_node_opq(
    nodes: np.ndarray,
    pq_m: int,
    nbits: int,
    sample_size: int = 250_000,
    niter: int = 12,
    opq_iters: int = 10,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = normalize_rows(nodes)
    dimension = int(values.shape[1])
    dsub, ksub = _validate_pq(dimension, pq_m, nbits)
    sample = _sample_nodes(values, sample_size, seed)
    rotation = np.eye(dimension, dtype=np.float32)
    codebook = np.empty((pq_m, ksub, dsub), dtype=np.float32)
    inner_iters = max(2, min(niter, 8))
    for step in range(opq_iters):
        rotated_sample = rotate_nodes(sample, rotation)
        codebook = _fit_pq_codebook(rotated_sample, pq_m, ksub, dsub, inner_iters, seed + 17 * (step + 1))
        reconstructed = reconstruct_pq(encode_pq(rotated_sample, codebook), codebook)
        rotation = _orthogonal_procrustes(sample, reconstructed)
    rotated_sample = rotate_nodes(sample, rotation)
    codebook = _fit_pq_codebook(rotated_sample, pq_m, ksub, dsub, niter, seed + 1009)
    codes = encode_pq(rotate_nodes(values, rotation), codebook)
    return rotation, codebook, codes


def train_node_sq(nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-dimension affine int8. Returns (scale, zero, uint8 codes)."""
    values = normalize_rows(nodes)
    vmin = values.min(axis=0)
    vmax = values.max(axis=0)
    scale = np.maximum(vmax - vmin, np.float32(1e-12)) / np.float32(255.0)
    zero = vmin
    codes = np.clip(np.rint((values - zero) / scale), 0, 255).astype(np.uint8)
    return scale.astype(np.float32), zero.astype(np.float32), codes


def decode_sq(codes: np.ndarray, scale: np.ndarray, zero: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(codes.astype(np.float32) * scale + zero, dtype=np.float32)


def iter_nodes(shards: list[tuple[Path, Path]]) -> Iterator[tuple[Path, dict, str, np.ndarray]]:
    for embedding_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        with np.load(embedding_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                if identifier not in arrays.files:
                    raise KeyError(f"{identifier} is in {manifest_path} but missing from {embedding_path}")
                yield embedding_path, manifest, identifier, normalize_rows(np.asarray(arrays[identifier]))


def collect_nodes(shards: list[tuple[Path, Path]]) -> tuple[np.ndarray, list[dict[str, Any]]]:
    blocks: list[np.ndarray] = []
    records: list[dict[str, Any]] = []
    for embedding_path, manifest, identifier, values in iter_nodes(shards):
        blocks.append(values)
        records.append(
            {
                "identifier": identifier,
                "length": int(values.shape[0]),
                "shard": embedding_path.name,
                "manifest": Path(manifest.get("source", "")).name if isinstance(manifest, dict) else "",
            }
        )
    if not blocks:
        raise ValueError("no node embeddings found")
    return np.concatenate(blocks, axis=0), records


def slice_code_windows(node_codes: np.ndarray, lengths: list[int], window_size: int, stride: int = 1) -> np.ndarray:
    if node_codes.ndim != 2:
        raise ValueError(f"node codes must be 2-D (n_nodes, pq_m), got {node_codes.shape}")
    if window_size < 1 or stride < 1:
        raise ValueError("window_size and stride must be positive")
    pq_m = int(node_codes.shape[1])
    windows: list[np.ndarray] = []
    offset = 0
    for length in lengths:
        codes = node_codes[offset : offset + length]
        offset += length
        if length < window_size:
            continue
        view = np.lib.stride_tricks.sliding_window_view(codes, window_size, axis=0)[::stride]
        windows.append(np.ascontiguousarray(np.transpose(view, (0, 2, 1))))
    if offset != node_codes.shape[0]:
        raise ValueError(f"record lengths sum to {offset} but node_codes has {node_codes.shape[0]} rows")
    if not windows:
        return np.empty((0, window_size, pq_m), dtype=node_codes.dtype)
    return np.ascontiguousarray(np.concatenate(windows, axis=0))


def slice_float_windows(nodes: np.ndarray, lengths: list[int], window_size: int, stride: int = 1) -> np.ndarray:
    values = normalize_rows(nodes)
    windows: list[np.ndarray] = []
    offset = 0
    for length in lengths:
        rows = values[offset : offset + length]
        offset += length
        if length < window_size:
            continue
        view = np.lib.stride_tricks.sliding_window_view(rows, window_size, axis=0)[::stride]
        windows.append(np.ascontiguousarray(np.transpose(view, (0, 2, 1))))
    if offset != values.shape[0]:
        raise ValueError(f"record lengths sum to {offset} but nodes has {values.shape[0]} rows")
    if not windows:
        return np.empty((0, window_size, values.shape[1]), dtype=np.float32)
    flat = np.concatenate(windows, axis=0).astype(np.float32, copy=False)
    return flat.reshape(flat.shape[0], flat.shape[1] * flat.shape[2])


def pack_window_codes(
    node_codes: np.ndarray,
    lengths: list[int],
    window_size: int,
    stride: int,
    pq_m: int,
    nbits: int,
) -> np.ndarray:
    windows = slice_code_windows(node_codes, lengths, window_size, stride)
    if windows.size == 0:
        width = window_size * node_code_bytes(pq_m, nbits)
        return np.empty((0, width), dtype=np.uint8)
    packed_nodes = pack_pq_codes(windows.reshape(-1, pq_m), pq_m, nbits)
    return packed_nodes.reshape(windows.shape[0], window_size * packed_nodes.shape[1])
