#!/usr/bin/env python3
"""Shared batched exact-reranking kernels for CPU and CuPy backends."""
from __future__ import annotations

import os
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np


INVALID_LABEL = np.int64(-1)


def _top_k(
    labels: np.ndarray,
    scores: np.ndarray,
    output_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return stable descending top-k rows, padding unavailable labels with -1."""
    if labels.ndim != 2 or scores.shape != labels.shape:
        raise ValueError("labels and scores must have the same 2-D shape")
    if output_k < 1 or output_k > labels.shape[1]:
        raise ValueError("output_k must be within the candidate block width")
    width = labels.shape[1]
    order = np.argpartition(-scores, min(output_k - 1, width - 1), axis=1)[:, :output_k]
    partial_scores = np.take_along_axis(scores, order, axis=1)
    stable_order = np.argsort(-partial_scores, axis=1, kind="stable")
    order = np.take_along_axis(order, stable_order, axis=1)
    return (
        np.take_along_axis(labels, order, axis=1).astype(np.int64, copy=False),
        np.take_along_axis(scores, order, axis=1).astype(np.float32, copy=False),
    )


def _merge_top_k(
    current_labels: np.ndarray,
    current_scores: np.ndarray,
    block_labels: np.ndarray,
    block_scores: np.ndarray,
    output_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.concatenate((current_labels, block_labels), axis=1)
    scores = np.concatenate((current_scores, block_scores), axis=1)
    return _top_k(labels, scores, output_k)


def _score_cpu_batch(
    database_windows: Any,
    query_windows: np.ndarray,
    labels: np.ndarray,
    output_k: int,
    candidate_batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Score one query batch while bounding candidate-vector materialization."""
    if labels.ndim != 2 or labels.shape[0] != query_windows.shape[0]:
        raise ValueError("candidate labels and query windows have incompatible shapes")
    if output_k < 1 or output_k > labels.shape[1]:
        raise ValueError("output_k must be within the candidate-label width")
    if candidate_batch_size < 1:
        raise ValueError("candidate_batch_size must be positive")

    batch_count, candidate_count = labels.shape
    final_labels = np.full((batch_count, output_k), INVALID_LABEL, dtype=np.int64)
    final_scores = np.full((batch_count, output_k), -np.inf, dtype=np.float32)
    database_size = int(database_windows.shape[0])
    for start in range(0, candidate_count, candidate_batch_size):
        stop = min(start + candidate_batch_size, candidate_count)
        block_labels = np.asarray(labels[:, start:stop], dtype=np.int64)
        valid = (block_labels >= 0) & (block_labels < database_size)
        safe_labels = np.where(valid, block_labels, 0)
        candidates = np.asarray(database_windows[safe_labels], dtype=np.float32)
        scores = np.einsum(
            "bd,bkd->bk",
            np.asarray(query_windows, dtype=np.float32),
            candidates,
            optimize=True,
        ).astype(np.float32, copy=False)
        scores[~valid] = -np.inf
        block_width = min(output_k, stop - start)
        block_top_labels, block_top_scores = _top_k(block_labels, scores, block_width)
        if block_width < output_k:
            padded_labels = np.full((batch_count, output_k), INVALID_LABEL, dtype=np.int64)
            padded_scores = np.full((batch_count, output_k), -np.inf, dtype=np.float32)
            padded_labels[:, :block_width] = block_top_labels
            padded_scores[:, :block_width] = block_top_scores
            block_top_labels, block_top_scores = padded_labels, padded_scores
        final_labels, final_scores = _merge_top_k(
            final_labels, final_scores, block_top_labels, block_top_scores, output_k
        )
    return final_labels, final_scores


def _score_gpu_batch(
    database_windows: Any,
    query_windows: np.ndarray,
    labels: np.ndarray,
    output_k: int,
    candidate_batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Score one query batch with CuPy, importing the GPU dependency lazily."""
    try:
        import cupy as cp
    except ImportError as exc:  # pragma: no cover - GPU environment specific
        raise RuntimeError("CUDA exact reranking requires CuPy") from exc
    if int(cp.cuda.runtime.getDeviceCount()) < 1:  # pragma: no cover - GPU environment specific
        raise RuntimeError("CUDA exact reranking requires a visible NVIDIA GPU")

    batch_count, candidate_count = labels.shape
    final_labels = np.full((batch_count, output_k), INVALID_LABEL, dtype=np.int64)
    final_scores = np.full((batch_count, output_k), -np.inf, dtype=np.float32)
    database_size = int(database_windows.shape[0])
    query_device = cp.asarray(np.asarray(query_windows, dtype=np.float32))
    for start in range(0, candidate_count, candidate_batch_size):
        stop = min(start + candidate_batch_size, candidate_count)
        block_labels = np.asarray(labels[:, start:stop], dtype=np.int64)
        valid = (block_labels >= 0) & (block_labels < database_size)
        safe_labels = np.where(valid, block_labels, 0)
        candidates = cp.asarray(np.asarray(database_windows[safe_labels], dtype=np.float32))
        scores = cp.einsum("bd,bkd->bk", query_device, candidates, optimize=True)
        scores = cp.asnumpy(scores).astype(np.float32, copy=False)
        scores[~valid] = -np.inf
        block_width = min(output_k, stop - start)
        block_top_labels, block_top_scores = _top_k(block_labels, scores, block_width)
        if block_width < output_k:
            padded_labels = np.full((batch_count, output_k), INVALID_LABEL, dtype=np.int64)
            padded_scores = np.full((batch_count, output_k), -np.inf, dtype=np.float32)
            padded_labels[:, :block_width] = block_top_labels
            padded_scores[:, :block_width] = block_top_scores
            block_top_labels, block_top_scores = padded_labels, padded_scores
        final_labels, final_scores = _merge_top_k(
            final_labels, final_scores, block_top_labels, block_top_scores, output_k
        )
    cp.cuda.Stream.null.synchronize()
    return final_labels, final_scores


def rerank_label_matrix(
    database_windows: Any,
    query_windows: np.ndarray,
    labels: np.ndarray,
    output_k: int,
    *,
    batch_size: int = 32,
    candidate_batch_size: int = 2048,
    workers: int = 1,
    device: str = "cpu",
) -> tuple[np.ndarray, np.ndarray, float]:
    """Exact-rerank candidate labels using normalized original window vectors.

    ``database_windows`` may be a NumPy array, memory map, or an object exposing
    ``shape`` and NumPy-style integer indexing.  The latter lets the pipeline
    materialize only the candidates selected by an ANN index instead of storing
    another full concatenated-window matrix.
    """
    queries = np.asarray(query_windows, dtype=np.float32)
    candidate_labels = np.asarray(labels, dtype=np.int64)
    if queries.ndim != 2 or candidate_labels.ndim != 2:
        raise ValueError("query windows and candidate labels must be 2-D")
    if candidate_labels.shape[0] != queries.shape[0]:
        raise ValueError("candidate labels and query windows have incompatible row counts")
    if output_k < 1 or output_k > candidate_labels.shape[1]:
        raise ValueError("output_k must be between 1 and the candidate count")
    if batch_size < 1 or candidate_batch_size < 1:
        raise ValueError("batch_size and candidate_batch_size must be positive")
    normalized_device = str(device).lower()
    if normalized_device == "auto":
        try:
            import cupy as cp

            normalized_device = "cuda" if int(cp.cuda.runtime.getDeviceCount()) > 0 else "cpu"
        except (ImportError, RuntimeError):
            normalized_device = "cpu"
    if normalized_device not in {"cpu", "cuda"}:
        raise ValueError("device must be cpu, cuda, or auto")
    if normalized_device == "cuda":
        workers = 1
    elif workers == 0:
        workers = max(1, int(os.cpu_count() or 1))
    elif workers < 1:
        raise ValueError("workers must be zero or positive")

    ranges = [
        (start, min(start + batch_size, queries.shape[0]))
        for start in range(0, queries.shape[0], batch_size)
    ]
    started = time.perf_counter()

    def run_one(bounds: tuple[int, int]) -> tuple[int, np.ndarray, np.ndarray]:
        start, stop = bounds
        if normalized_device == "cuda":
            result_labels, result_scores = _score_gpu_batch(
                database_windows,
                queries[start:stop],
                candidate_labels[start:stop],
                output_k,
                candidate_batch_size,
            )
        else:
            result_labels, result_scores = _score_cpu_batch(
                database_windows,
                queries[start:stop],
                candidate_labels[start:stop],
                output_k,
                candidate_batch_size,
            )
        return start, result_labels, result_scores

    final_labels = np.full((queries.shape[0], output_k), INVALID_LABEL, dtype=np.int64)
    final_scores = np.full((queries.shape[0], output_k), -np.inf, dtype=np.float32)
    if normalized_device == "cuda" or workers == 1 or len(ranges) < 2:
        completed = [run_one(bounds) for bounds in ranges]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            completed = list(executor.map(run_one, ranges))
    for start, batch_labels, batch_scores in completed:
        stop = start + batch_labels.shape[0]
        final_labels[start:stop] = batch_labels
        final_scores[start:stop] = batch_scores
    return final_labels, final_scores, time.perf_counter() - started
