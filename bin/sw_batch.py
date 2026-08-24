#!/usr/bin/env python3
"""Batch Smith–Waterman over independent score matrices.

Uses GINFINITY-SW ``align_scores_many`` when present (1.2+). The 1.1
container falls back to the same Numba ``prange`` kernel against ``_core``.
"""
from __future__ import annotations

import os
from collections.abc import Iterable

import numpy as np
from ginfinity_sw import HAVE_NUMBA as SW_HAVE_NUMBA
from ginfinity_sw import ScoringParameters, align_scores
from ginfinity_sw.core import _core

try:
    from ginfinity_sw import align_scores_many as _library_align_scores_many
except ImportError:  # ginfinity-sw < 1.2
    _library_align_scores_many = None

try:
    from numba import njit, prange, set_num_threads
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return lambda function: function

    def prange(*args):
        return range(*args)

    def set_num_threads(_n: int) -> None:
        return None


_SCORE_BATCH_BYTES = 64 * 1024 * 1024


@njit(parallel=True, cache=True)
def _core_many(scores, q_len, t_len, gap_open, gap_extend, traceback):
    n = scores.shape[0]
    bests = np.zeros(n)
    counts = np.zeros(n, dtype=np.int64)
    max_path = scores.shape[1] + scores.shape[2] + 2
    buffers = np.empty((n, max_path, 2), dtype=np.int64)
    for index in prange(n):
        query_length = int(q_len[index])
        target_length = int(t_len[index])
        if query_length == 0 or target_length == 0:
            continue
        blocked_query = np.zeros(query_length, dtype=np.bool_)
        blocked_target = np.zeros(target_length, dtype=np.bool_)
        best, count, buffer = _core(
            scores[index, :query_length, :target_length],
            gap_open, gap_extend, traceback, blocked_query, blocked_target)
        bests[index] = best
        counts[index] = count
        if traceback and count > 0:
            for step in range(count):
                buffers[index, step, 0] = buffer[step, 0]
                buffers[index, step, 1] = buffer[step, 1]
    return bests, counts, buffers


def _alignment_from_core(best, count, buffer, rows, traceback):
    from ginfinity_sw import Alignment

    if not traceback or best <= 0.0:
        return Alignment(float(best), (0, 0), (0, 0), (), rows)
    columns = tuple(
        (int(buffer[index, 0]), int(buffer[index, 1]))
        for index in range(int(count) - 1, -1, -1))
    query_indices = [query for query, _ in columns if query >= 0]
    target_indices = [target for _, target in columns if target >= 0]
    query_span = ((min(query_indices), max(query_indices) + 1)
                  if query_indices else (0, 0))
    target_span = ((min(target_indices), max(target_indices) + 1)
                   if target_indices else (0, 0))
    return Alignment(float(best), query_span, target_span, columns, rows)


def _score_batch_slices(shapes: list[tuple[int, int]]):
    start = 0
    n = len(shapes)
    while start < n:
        max_query, max_target = shapes[start]
        stop = start + 1
        while stop < n:
            query = max(max_query, shapes[stop][0])
            target = max(max_target, shapes[stop][1])
            cells = query * target * (stop - start + 1)
            if cells * 8 > _SCORE_BATCH_BYTES and stop > start:
                break
            max_query, max_target = query, target
            stop += 1
        yield start, stop, max_query, max_target
        start = stop


def _fallback_align_scores_many(
    score_matrices: list[np.ndarray],
    parameters: ScoringParameters,
    *,
    traceback: bool,
    max_cells: int,
    workers: int,
):
    if workers <= 1 or len(score_matrices) <= 1 or not HAVE_NUMBA or not SW_HAVE_NUMBA:
        return [
            align_scores(matrix, parameters, traceback=traceback, max_cells=max_cells)
            for matrix in score_matrices
        ]
    thread_count = max(1, workers)
    try:
        from numba.config import NUMBA_NUM_THREADS
        thread_count = min(thread_count, int(NUMBA_NUM_THREADS))
    except Exception:
        pass
    set_num_threads(thread_count)
    results = []
    gap_open = float(parameters.gap_open)
    gap_extend = float(parameters.gap_extend)
    do_traceback = bool(traceback)
    shapes = [matrix.shape for matrix in score_matrices]
    for start, stop, max_query, max_target in _score_batch_slices(shapes):
        batch = stop - start
        padded = np.zeros((batch, max_query, max_target), dtype=np.float64)
        q_len = np.empty(batch, dtype=np.int64)
        t_len = np.empty(batch, dtype=np.int64)
        for offset, matrix in enumerate(score_matrices[start:stop]):
            query_length, target_length = matrix.shape
            padded[offset, :query_length, :target_length] = matrix
            q_len[offset] = query_length
            t_len[offset] = target_length
        bests, counts, buffers = _core_many(
            padded, q_len, t_len, gap_open, gap_extend, do_traceback)
        for offset in range(batch):
            results.append(_alignment_from_core(
                bests[offset], counts[offset], buffers[offset],
                int(q_len[offset]), do_traceback))
    return results


def align_score_matrices(
    score_matrices: Iterable[np.ndarray],
    parameters: ScoringParameters,
    *,
    traceback: bool = True,
    max_cells: int = 16_777_216,
    workers: int = 1,
):
    prepared = [np.ascontiguousarray(np.asarray(matrix, dtype=np.float64))
                for matrix in score_matrices]
    if not prepared:
        return []
    thread_count = max(1, int(workers or 1))
    if _library_align_scores_many is not None:
        return _library_align_scores_many(
            prepared,
            parameters,
            traceback=traceback,
            max_cells=max_cells,
            workers=thread_count,
        )
    return _fallback_align_scores_many(
        prepared,
        parameters,
        traceback=traceback,
        max_cells=max_cells,
        workers=thread_count,
    )


def pin_blas_threads() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(name, "1")
