"""Node-level PQ + CAGRA search with shared-memory ADC."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from . import pq_cagra_adc_ext as _ext

Index = _ext.Index


def normalize_rows(values: np.ndarray) -> np.ndarray:
    rows = np.ascontiguousarray(values, dtype=np.float32)
    if rows.ndim != 2:
        raise ValueError(f"expected a 2-D matrix, got {rows.shape}")
    norms = np.linalg.norm(rows, axis=1, keepdims=True)
    return rows / np.maximum(norms, np.float32(1e-12))


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
    """Assign L2-nearest subcodes. `values` must already be in codebook space."""
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


def reconstruction_mse(values: np.ndarray, codes: np.ndarray, codebook: np.ndarray) -> float:
    residual = np.ascontiguousarray(values, dtype=np.float32) - reconstruct_pq(codes, codebook)
    return float(np.mean(residual * residual))


def train_node_pq(
    nodes: np.ndarray,
    pq_m: int,
    nbits: int,
    sample_size: int = 50_000,
    niter: int = 20,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a product quantizer on L2-normalized node embeddings.

    Returns `(codebook, codes)` where codebook has shape `(pq_m, 2**nbits, dsub)`
    and codes has shape `(n_nodes, pq_m)` as uint8 subcodes.
    """
    values = normalize_rows(nodes)
    dsub, ksub = _validate_pq(int(values.shape[1]), pq_m, nbits)
    sample = _sample_nodes(values, sample_size, seed)
    codebook = _fit_pq_codebook(sample, pq_m, ksub, dsub, niter, seed)
    return codebook, encode_pq(values, codebook)


def _orthogonal_procrustes(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Return orthogonal R minimizing ||source @ R - target||_F."""
    gram = source.T @ target
    u_mat, _singular, vt_mat = np.linalg.svd(gram, full_matrices=True)
    rotation = u_mat @ vt_mat
    if np.linalg.det(rotation) < 0:
        u_mat[:, -1] *= -1
        rotation = u_mat @ vt_mat
    return np.ascontiguousarray(rotation, dtype=np.float32)


def rotate_nodes(nodes: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    values = np.ascontiguousarray(nodes, dtype=np.float32)
    matrix = np.ascontiguousarray(rotation, dtype=np.float32)
    if values.ndim != 2 or matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("rotation must be a square matrix applied to 2-D nodes")
    if values.shape[1] != matrix.shape[0]:
        raise ValueError(f"node dim {values.shape[1]} does not match rotation {matrix.shape}")
    return np.ascontiguousarray(values @ matrix, dtype=np.float32)


def rotate_windows(windows: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    """Apply a per-node rotation to windows of shape (n, W, D) or (n, W*D)."""
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


def train_node_opq(
    nodes: np.ndarray,
    pq_m: int,
    nbits: int,
    sample_size: int = 50_000,
    niter: int = 12,
    opq_iters: int = 10,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit Optimized PQ (Ge et al.): orthogonal rotation then subspace PQ.

    Returns `(rotation, codebook, codes)` with `rotation` shaped `(D, D)`.
    Database codes are of rotated nodes `x @ R`. Queries must be rotated
    with `rotate_windows` before ADC search. Exact rerank still uses the
    original unrotated windows.
    """
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


def slice_code_windows(node_codes: np.ndarray, lengths: list[int], window_size: int, stride: int = 1) -> np.ndarray:
    """Pack per-record node codes into windows of shape `(n_windows, window_size, pq_m)`."""
    if node_codes.ndim != 2:
        raise ValueError(f"node codes must be 2-D (n_nodes, pq_m), got {node_codes.shape}")
    if window_size < 1 or stride < 1:
        raise ValueError("window_size and stride must be positive")
    pq_m = int(node_codes.shape[1])
    windows: list[np.ndarray] = []
    offset = 0
    for length in lengths:
        if length < 0:
            raise ValueError("record lengths must be non-negative")
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
    return np.concatenate(windows, axis=0).astype(np.float32, copy=False)


def adc_lut(query_windows: np.ndarray, codebook: np.ndarray) -> np.ndarray:
    """Return ADC lookup tables shaped `(nq, window_size * pq_m, ksub)` as float32."""
    queries = np.ascontiguousarray(query_windows, dtype=np.float32)
    book = np.ascontiguousarray(codebook, dtype=np.float32)
    if queries.ndim != 3:
        raise ValueError("query windows must have shape (nq, window_size, dim)")
    if book.ndim != 3:
        raise ValueError("codebook must have shape (pq_m, ksub, dsub)")
    nq, window_size, dim = queries.shape
    pq_m, ksub, dsub = book.shape
    if dim != pq_m * dsub:
        raise ValueError(f"query dim {dim} does not match codebook {pq_m}*{dsub}")
    reshaped = queries.reshape(nq, window_size, pq_m, dsub)
    # lut[q, p, m, c] = -dot(q[p,m], codebook[m,c])
    lut = -np.einsum("q p m d, m c d -> q p m c", reshaped, book, optimize=True)
    return np.ascontiguousarray(lut.reshape(nq, window_size * pq_m, ksub), dtype=np.float32)


def adc_distances(query_windows: np.ndarray, codebook: np.ndarray, codes: np.ndarray) -> np.ndarray:
    """Exact ADC distances, smaller is closer. Shape `(nq, n)`."""
    lut = adc_lut(query_windows, codebook)
    packed = np.ascontiguousarray(codes)
    if packed.ndim == 3:
        packed = packed.reshape(packed.shape[0], packed.shape[1] * packed.shape[2])
    if packed.ndim != 2 or packed.shape[1] != lut.shape[1]:
        raise ValueError(f"codes shape {packed.shape} is incompatible with LUT {lut.shape}")
    positions = np.arange(packed.shape[1])
    distances = np.empty((lut.shape[0], packed.shape[0]), dtype=np.float32)
    for query_id in range(lut.shape[0]):
        distances[query_id] = lut[query_id, positions, packed].sum(axis=1)
    return distances


def build_index(
    codes: np.ndarray,
    codebook: np.ndarray,
    *,
    graph_degree: int = 32,
    intermediate_graph_degree: int = 64,
    nndescent_iterations: int = 0,
) -> Any:
    packed = np.ascontiguousarray(codes, dtype=np.uint8)
    book = np.ascontiguousarray(codebook, dtype=np.float32)
    if packed.ndim == 3:
        window_size, pq_m = int(packed.shape[1]), int(packed.shape[2])
    elif packed.ndim == 2:
        if book.ndim != 3:
            raise ValueError("codebook must be 3-D when codes are 2-D")
        pq_m = int(book.shape[0])
        if packed.shape[1] % pq_m != 0:
            raise ValueError("2-D code width must be divisible by pq_m")
        window_size = int(packed.shape[1] // pq_m)
    else:
        raise ValueError("codes must be 2-D or 3-D")
    index = Index()
    index.build(
        packed,
        book,
        window_size=window_size,
        pq_m=pq_m,
        nbits=int(np.log2(book.shape[1])),
        dsub=int(book.shape[2]),
        graph_degree=graph_degree,
        intermediate_graph_degree=intermediate_graph_degree,
        nndescent_iterations=nndescent_iterations,
    )
    return index


def search(
    index: Any,
    query_windows: np.ndarray,
    k: int = 10,
    *,
    device: str = "cuda",
    rotation: np.ndarray | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, np.ndarray]:
    queries = np.ascontiguousarray(query_windows, dtype=np.float32)
    if queries.ndim == 2:
        dim = int(index.dim)
        window_size = int(index.window_size)
        if queries.shape[1] != window_size * dim:
            raise ValueError(
                f"flat queries have width {queries.shape[1]}, expected {window_size * dim}"
            )
        queries = queries.reshape(queries.shape[0], window_size, dim)
    if rotation is not None:
        queries = rotate_windows(queries, rotation)
    kind = str(device).lower()
    if kind in {"cpu", "host"}:
        return index.search_cpu(queries, k=k, **kwargs)
    if kind in {"cuda", "gpu"}:
        return index.search(queries, k=k, **kwargs)
    raise ValueError("device must be cuda or cpu")


def save_index(index: Any, path: str | Path) -> None:
    index.save(str(path))


def load_index(path: str | Path) -> Any:
    return Index.load(str(path))
