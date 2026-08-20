#!/usr/bin/env python3
"""Benchmark FAISS HNSW over fixed-scale int8 window vectors.

The CAGRA experiment uses the original normalized window vectors, rounded to
int8 with one global scale, for candidate selection. FAISS's regular
``IndexHNSWFlat`` accepts float32, so this harness measures both a true
``IndexHNSWSQ`` direct-signed-int8 index and a float-backed control containing
the same integer-valued coordinates. Exact reranking always uses the original
window matrix.
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def quantize(values: np.ndarray, scale: float) -> np.ndarray:
    scaled = np.rint(np.asarray(values, dtype=np.float32) * np.float32(scale))
    if scaled.size and (scaled.min() < -127 or scaled.max() > 127):
        raise ValueError("scaled values exceed the signed int8 range")
    return np.ascontiguousarray(scaled, dtype=np.float32)


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
    if labels.ndim != 2 or labels.shape[0] != query_windows.shape[0]:
        raise ValueError("candidate labels and queries have incompatible shapes")
    if output_k < 1 or output_k > labels.shape[1]:
        raise ValueError("output_k must be within the candidate-label width")
    final_labels = np.empty((labels.shape[0], output_k), dtype=np.int64)
    final_scores = np.empty((labels.shape[0], output_k), dtype=np.float32)
    started = time.perf_counter()
    for start in range(0, labels.shape[0], batch_size):
        stop = min(start + batch_size, labels.shape[0])
        batch_labels = np.asarray(labels[start:stop], dtype=np.int64)
        if (batch_labels < 0).any() or (batch_labels >= database_windows.shape[0]).any():
            raise ValueError("candidate labels contain an invalid database ID")
        candidates = np.asarray(database_windows[batch_labels], dtype=np.float32)
        scores = np.einsum(
            "bd,bkd->bk", query_windows[start:stop], candidates, optimize=True
        )
        order = np.argsort(-scores, axis=1, kind="stable")[:, :output_k]
        rows = np.arange(stop - start)[:, None]
        final_labels[start:stop] = batch_labels[rows, order]
        final_scores[start:stop] = scores[rows, order]
    return final_labels, final_scores, time.perf_counter() - started


def max_rss_bytes() -> int:
    # Linux reports ru_maxrss in KiB; the benchmark image is Linux-only.
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def make_index(faiss: Any, dim: int, variant: str, m: int, metric: int) -> Any:
    if variant == "sq-int8":
        qtype = faiss.ScalarQuantizer.QT_8bit_direct_signed
        index = faiss.IndexHNSWSQ(dim, qtype, m, metric)
    elif variant in {"flat-int8", "flat-fp32"}:
        index = faiss.IndexHNSWFlat(dim, m, metric)
    else:
        raise ValueError(f"unknown FAISS HNSW variant: {variant}")
    return index


def build_index(
    faiss: Any,
    source: np.ndarray,
    variant: str,
    scale: float,
    m: int,
    ef_construction: int,
    add_batch_size: int,
    metric: int,
) -> tuple[Any, float, float]:
    index = make_index(faiss, source.shape[1], variant, m, metric)
    index.hnsw.efConstruction = ef_construction
    started = time.perf_counter()
    if not index.is_trained:
        train_source = np.asarray(source[: min(add_batch_size, source.shape[0])], dtype=np.float32)
        if variant == "sq-int8":
            train_source = quantize(train_source, scale)
        index.train(np.ascontiguousarray(train_source, dtype=np.float32))
    train_seconds = time.perf_counter() - started

    add_started = time.perf_counter()
    for start in range(0, source.shape[0], add_batch_size):
        stop = min(start + add_batch_size, source.shape[0])
        batch = np.asarray(source[start:stop], dtype=np.float32)
        if variant in {"sq-int8", "flat-int8"}:
            batch = quantize(batch, scale)
        index.add(np.ascontiguousarray(batch, dtype=np.float32))
    add_seconds = time.perf_counter() - add_started
    return index, train_seconds, add_seconds


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    if args.m < 2 or args.ef_construction < 1 or args.add_batch_size < 1:
        raise ValueError("M, ef-construction, and add-batch-size must be positive")
    if args.num_threads < 1:
        raise ValueError("num-threads must be positive")
    efs = [int(value) for value in args.ef_search_values.split(",") if value.strip()]
    if not efs or any(ef < args.candidate_k for ef in efs):
        raise ValueError("every ef-search value must be >= candidate-k")
    if args.output_k < 1 or args.output_k > args.candidate_k:
        raise ValueError("output-k must be between 1 and candidate-k")

    import faiss

    faiss.omp_set_num_threads(args.num_threads)
    source = np.load(args.database_windows, mmap_mode="r")
    if source.ndim != 2 or source.shape[0] < 1:
        raise ValueError("database windows must be a non-empty 2-D array")
    if args.max_vectors is not None:
        source = source[: args.max_vectors]
    queries = np.asarray(np.load(args.queries), dtype=np.float32)
    if queries.ndim != 2 or queries.shape[1] != source.shape[1]:
        raise ValueError("query and database dimensions do not match")
    reference = None
    if args.reference_labels is not None:
        reference = np.asarray(np.load(args.reference_labels), dtype=np.int64)
        if reference.shape[0] != queries.shape[0]:
            raise ValueError("reference and query row counts do not match")

    metric_name = args.metric.lower()
    metric = faiss.METRIC_INNER_PRODUCT if metric_name == "ip" else faiss.METRIC_L2
    args.outdir.mkdir(parents=True, exist_ok=True)
    loaded_existing_index = args.load_index is not None
    if loaded_existing_index:
        index_path = args.load_index
        load_started = time.perf_counter()
        index = faiss.read_index(str(index_path))
        load_seconds = time.perf_counter() - load_started
        if index.ntotal != source.shape[0] or index.d != source.shape[1]:
            raise ValueError("loaded index and database windows are incompatible")
        train_seconds = None
        add_seconds = None
        build_seconds = None
        serialize_seconds = None
    else:
        build_started = time.perf_counter()
        index, train_seconds, add_seconds = build_index(
            faiss,
            source,
            args.variant,
            args.int8_scale,
            args.m,
            args.ef_construction,
            args.add_batch_size,
            metric,
        )
        build_seconds = time.perf_counter() - build_started

        index_path = args.outdir / "index.faiss"
        serialize_started = time.perf_counter()
        faiss.write_index(index, str(index_path))
        serialize_seconds = time.perf_counter() - serialize_started
        load_seconds = None

    query_started = time.perf_counter()
    query_vectors = quantize(queries, args.int8_scale)
    if args.variant == "flat-fp32":
        query_vectors = np.ascontiguousarray(queries, dtype=np.float32)
    query_transform_seconds = time.perf_counter() - query_started

    rows: list[dict[str, Any]] = []
    for ef in efs:
        index.hnsw.efSearch = ef
        index.search(query_vectors[:1], args.candidate_k)
        search_started = time.perf_counter()
        distances, labels = index.search(query_vectors, args.candidate_k)
        search_seconds = time.perf_counter() - search_started
        labels = np.asarray(labels, dtype=np.int64)
        distances = np.asarray(distances, dtype=np.float32)
        np.save(args.outdir / f"labels-{ef}.npy", labels)
        np.save(args.outdir / f"distances-{ef}.npy", distances)
        final_labels, final_scores, rerank_seconds = rerank_candidates(
            source, queries, labels, args.output_k, args.rerank_batch_size
        )
        row: dict[str, Any] = {
            "ef_search": ef,
            "candidate_k": args.candidate_k,
            "output_k": args.output_k,
            "search_seconds": search_seconds,
            "exact_rerank_seconds": rerank_seconds,
            "search_plus_rerank_seconds": search_seconds + rerank_seconds,
            "mean_candidate_distance_or_score": float(np.mean(distances[:, 0])),
            "mean_final_score": float(np.mean(final_scores[:, 0])),
        }
        if reference is not None:
            for k in (1, 5, 10, 50):
                row[f"candidate_recall_at_{k}"] = recall_at(reference, labels, k)
                row[f"final_recall_at_{k}"] = recall_at(reference, final_labels, k)
        rows.append(row)

    result: dict[str, Any] = {
        "backend": "faiss_hnsw",
        "faiss_version": faiss.__version__,
        "variant": args.variant,
        "metric": metric_name,
        "candidate_representation": (
            "original_normalized_window_int8" if args.variant != "flat-fp32" else "original_normalized_window_float32"
        ),
        "int8_scale": args.int8_scale if args.variant != "flat-fp32" else None,
        "database_windows": str(args.database_windows),
        "queries": str(args.queries),
        "reference_labels": str(args.reference_labels) if args.reference_labels else None,
        "n_windows": int(source.shape[0]),
        "n_queries": int(queries.shape[0]),
        "window_dim": int(source.shape[1]),
        "m": args.m,
        "ef_construction": args.ef_construction,
        "num_threads": args.num_threads,
        "add_batch_size": args.add_batch_size,
        "train_seconds": train_seconds,
        "add_seconds": add_seconds,
        "build_seconds": build_seconds,
        "serialize_seconds": serialize_seconds,
        "query_transform_seconds": query_transform_seconds,
        "loaded_existing_index": loaded_existing_index,
        "load_seconds": load_seconds,
        "index_bytes": index_path.stat().st_size,
        "peak_rss_bytes": max_rss_bytes(),
        "peak_rss_mib": max_rss_bytes() / (1024 * 1024),
        "index_path": str(index_path),
        "results": rows,
    }
    (args.outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-windows", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--reference-labels", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--load-index",
        type=Path,
        help="reuse an existing FAISS index instead of rebuilding it",
    )
    parser.add_argument("--variant", choices=("sq-int8", "flat-int8", "flat-fp32"), default="sq-int8")
    parser.add_argument("--metric", choices=("ip", "l2"), default="ip")
    parser.add_argument("--int8-scale", type=float, default=850.0)
    parser.add_argument("--m", type=int, default=32)
    parser.add_argument("--ef-construction", type=int, default=200)
    parser.add_argument("--ef-search-values", default="50,100,200,400,800")
    parser.add_argument("--candidate-k", type=int, default=50)
    parser.add_argument("--output-k", type=int, default=50)
    parser.add_argument("--num-threads", type=int, default=16)
    parser.add_argument("--add-batch-size", type=int, default=32768)
    parser.add_argument("--rerank-batch-size", type=int, default=32)
    parser.add_argument("--max-vectors", type=int)
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
