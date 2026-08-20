#!/usr/bin/env python3
"""Evaluate labels emitted by the cuVS CAGRA-to-CPU-HNSW research drivers.

The conversion and CPU search are implemented in the companion C++ programs
``cagra_to_hnsw_cpu.cpp`` and ``search_cuvs_hnsw_cpu.cpp``. This script keeps
the reference comparison and exact original-window rerank in Python, where it
can reuse NumPy memory maps without loading the full window matrix in RAM.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def parse_int_list(value: str) -> list[int]:
    values = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not values or any(item < 1 for item in values):
        raise ValueError("integer lists must contain positive values")
    return values


def quantize_windows(values: np.ndarray, scale: float) -> np.ndarray:
    """Quantize normalized windows for the int8 cuVS driver input."""
    rows = np.asarray(values, dtype=np.float32)
    if rows.ndim != 2:
        raise ValueError(f"windows must be a 2-D array, got {rows.shape}")
    if not np.isfinite(rows).all():
        raise ValueError("windows contain non-finite values")
    scaled = np.rint(rows * np.float32(scale))
    if scaled.size and (scaled.min() < -127 or scaled.max() > 127):
        raise ValueError("scaled windows exceed the signed int8 range")
    return np.ascontiguousarray(scaled, dtype=np.int8)


def recall_at(reference: np.ndarray, approximate: np.ndarray, k: int) -> float | None:
    if k < 1 or reference.shape[1] < k or approximate.shape[1] < k:
        return None
    values = [
        len(set(expected.tolist()) & set(actual.tolist())) / float(k)
        for expected, actual in zip(reference[:, :k], approximate[:, :k])
    ]
    return float(np.mean(values)) if values else None


def rerank_candidates(
    database_windows: np.ndarray,
    query_windows: np.ndarray,
    labels: np.ndarray,
    output_k: int,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Score candidates with original windows and return exact top-k labels."""
    if labels.ndim != 2 or labels.shape[0] != query_windows.shape[0]:
        raise ValueError("candidate labels and queries have incompatible shapes")
    if output_k < 1 or output_k > labels.shape[1]:
        raise ValueError("output_k must be within the candidate-label width")
    if batch_size < 1:
        raise ValueError("rerank batch size must be positive")

    final_labels = np.empty((labels.shape[0], output_k), dtype=np.int64)
    final_scores = np.empty((labels.shape[0], output_k), dtype=np.float32)
    started = time.perf_counter()
    database = np.asarray(database_windows)
    queries = np.asarray(query_windows, dtype=np.float32)
    for start in range(0, labels.shape[0], batch_size):
        stop = min(start + batch_size, labels.shape[0])
        batch_labels = np.asarray(labels[start:stop], dtype=np.int64)
        if (batch_labels < 0).any() or (batch_labels >= database.shape[0]).any():
            raise ValueError("candidate labels contain an invalid database ID")
        candidates = np.asarray(database[batch_labels], dtype=np.float32)
        scores = np.einsum("bd,bkd->bk", queries[start:stop], candidates, optimize=True)
        order = np.argsort(-scores, axis=1, kind="stable")[:, :output_k]
        row_ids = np.arange(stop - start)[:, None]
        final_labels[start:stop] = batch_labels[row_ids, order]
        final_scores[start:stop] = scores[row_ids, order]
    return final_labels, final_scores, time.perf_counter() - started


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    efs = parse_int_list(args.ef_search_values)
    if args.candidate_k < 1 or args.output_k < 1 or args.output_k > args.candidate_k:
        raise ValueError("candidate-k/output-k values are incompatible")

    queries = np.asarray(np.load(args.queries), dtype=np.float32)
    database = np.load(args.database_windows, mmap_mode="r")
    if queries.ndim != 2 or database.ndim != 2 or queries.shape[1] != database.shape[1]:
        raise ValueError("database windows and queries must be compatible 2-D arrays")
    reference = None
    if args.reference_labels is not None:
        reference = np.asarray(np.load(args.reference_labels), dtype=np.int64)
        if reference.ndim != 2 or reference.shape[0] != queries.shape[0]:
            raise ValueError("reference labels and queries have incompatible shapes")

    search_metrics: dict[int, dict[str, Any]] = {}
    if args.search_metrics is not None:
        raw_metrics = json.loads(args.search_metrics.read_text())
        search_metrics = {
            int(row["ef_search"]): row for row in raw_metrics.get("results", [])
        }

    rows: list[dict[str, Any]] = []
    for ef in efs:
        label_path = Path(f"{args.labels_prefix}{ef}.u64")
        raw_labels = np.fromfile(label_path, dtype=np.uint64)
        expected = queries.shape[0] * args.candidate_k
        if raw_labels.size != expected:
            raise ValueError(f"{label_path} has {raw_labels.size} labels; expected {expected}")
        labels = raw_labels.reshape(queries.shape[0], args.candidate_k).astype(np.int64)
        final_labels, final_scores, rerank_seconds = rerank_candidates(
            database, queries, labels, args.output_k, args.rerank_batch_size
        )
        row: dict[str, Any] = {
            "ef_search": ef,
            "candidate_k": args.candidate_k,
            "output_k": args.output_k,
            "exact_rerank_seconds": rerank_seconds,
            "mean_final_score": float(np.mean(final_scores[:, 0])),
        }
        if ef in search_metrics:
            row.update(
                {
                    "search_seconds": search_metrics[ef].get("search_seconds"),
                    "search_plus_rerank_seconds": (
                        float(search_metrics[ef]["search_seconds"]) + rerank_seconds
                    ),
                }
            )
        if reference is not None:
            for k in (1, 5, 10, 50):
                row[f"candidate_recall_at_{k}"] = recall_at(reference, labels, k)
                row[f"final_recall_at_{k}"] = recall_at(reference, final_labels, k)
        rows.append(row)

    result: dict[str, Any] = {
        "backend": "cuvs_cagra_to_cpu_hnsw_evaluation",
        "database_windows": str(args.database_windows),
        "queries": str(args.queries),
        "reference_labels": str(args.reference_labels) if args.reference_labels else None,
        "n_windows": int(database.shape[0]),
        "n_queries": int(queries.shape[0]),
        "window_dim": int(database.shape[1]),
        "results": rows,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-windows", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--reference-labels", type=Path)
    parser.add_argument("--labels-prefix", required=True)
    parser.add_argument("--search-metrics", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidate-k", type=int, required=True)
    parser.add_argument("--output-k", type=int, default=50)
    parser.add_argument("--ef-search-values", default="50")
    parser.add_argument("--rerank-batch-size", type=int, default=32)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = evaluate(args)
        print(json.dumps(result, indent=2))
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
