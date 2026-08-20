#!/usr/bin/env python3
"""Benchmark compact C++ custom-distance hnswlib against exact FlatIP."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_hnswlib import (
    build_codes,
    build_flat_index,
    build_flat_windows,
    file_size,
    fit_centroids,
    load_queries,
    load_structures,
    load_or_fit_centroids,
    normalize_rows,
    parse_int_list,
    positive_rank,
    recall_at,
    self_hit_rate,
    source_embeddings,
)


def log(message: str) -> None:
    print(message, flush=True)


def write_raw_window_codes(
    codes: np.ndarray,
    records,
    window_size: int,
    stride: int,
    output_path: Path,
) -> int:
    n_windows = records[-1].window_start + records[-1].n_windows
    expected_shape = (n_windows, window_size)
    if output_path.exists() and output_path.stat().st_size == n_windows * window_size * np.dtype(np.uint16).itemsize:
        return n_windows
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output = np.memmap(output_path, mode="w+", dtype=np.uint16, shape=expected_shape)
    for record in records:
        if record.n_windows == 0:
            continue
        values = np.asarray(codes[record.node_start : record.node_start + record.length])
        views = np.lib.stride_tricks.sliding_window_view(values, window_size)[::stride]
        output[record.window_start : record.window_start + record.n_windows] = np.asarray(views, dtype=np.uint16)
    output.flush()
    del output
    return n_windows


def query_code_windows(codes: np.ndarray, queries, records, window_size: int) -> np.ndarray:
    by_id = {record.identifier: record for record in records}
    output = np.empty((len(queries), window_size), dtype=np.uint16)
    for row, query in enumerate(queries):
        record = by_id[query.transcript_id]
        start = record.node_start + query.window_offset
        output[row] = np.asarray(codes[start : start + window_size], dtype=np.uint16)
    return output


def parse_float_list(raw: str) -> list[float]:
    values = [float(value.strip()) for value in str(raw).split(",") if value.strip()]
    if not values:
        raise ValueError("threshold list must contain at least one value")
    return values


def write_similarity(centroids: np.ndarray, output_path: Path) -> None:
    values = np.ascontiguousarray(centroids, dtype=np.float32)
    similarity = np.ascontiguousarray(values @ values.T, dtype=np.float32)
    similarity.tofile(output_path)


def run_compact_search(
    executable: Path,
    index_path: Path,
    similarity_path: Path,
    query_codes_path: Path,
    query_count: int,
    window_size: int,
    n_centroids: int,
    top_k: int,
    ef_search: int,
    threads: int,
    workdir: Path,
) -> tuple[np.ndarray, np.ndarray, float]:
    result_prefix = index_path.stem
    labels_path = workdir / f"{result_prefix}_ef{ef_search}_k{top_k}_labels.bin"
    distances_path = workdir / f"{result_prefix}_ef{ef_search}_k{top_k}_distances.bin"
    expected_labels = query_count * top_k * np.dtype("<u8").itemsize
    expected_distances = query_count * top_k * np.dtype("<f4").itemsize
    started = time.perf_counter()
    if (
        not labels_path.exists()
        or not distances_path.exists()
        or labels_path.stat().st_size != expected_labels
        or distances_path.stat().st_size != expected_distances
    ):
        command = [
            str(executable),
            "search",
            "--codes", str(query_codes_path),
            "--similarity", str(similarity_path),
            "--query-count", str(query_count),
            "--window-size", str(window_size),
            "--n-centroids", str(n_centroids),
            "--index", str(index_path),
            "--k", str(top_k),
            "--ef-search", str(ef_search),
            "--threads", str(threads),
            "--labels-out", str(labels_path),
            "--distances-out", str(distances_path),
        ]
        subprocess.run(command, check=True)
    elapsed = time.perf_counter() - started
    labels = np.fromfile(labels_path, dtype="<u8").reshape(query_count, top_k)
    distances = np.fromfile(distances_path, dtype="<f4").reshape(query_count, top_k)
    return labels, distances, elapsed


def rerank_candidates(
    windows: np.ndarray,
    query_ids: np.ndarray,
    candidate_labels: np.ndarray,
    output_k: int,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Rerank a compact-HNSW candidate pool with original window vectors."""
    result = np.empty((candidate_labels.shape[0], output_k), dtype=np.int64)
    result_scores = np.empty((candidate_labels.shape[0], output_k), dtype=np.float32)
    query_vectors = np.ascontiguousarray(windows[query_ids], dtype=np.float32)
    for start in range(0, candidate_labels.shape[0], batch_size):
        stop = min(start + batch_size, candidate_labels.shape[0])
        labels = candidate_labels[start:stop]
        values = np.ascontiguousarray(windows[labels], dtype=np.float32)
        scores = np.einsum("qd,qkd->qk", query_vectors[start:stop], values, optimize=True)
        order = np.argpartition(-scores, output_k - 1, axis=1)[:, :output_k]
        order_scores = np.take_along_axis(scores, order, axis=1)
        order = np.take_along_axis(order, np.argsort(-order_scores, axis=1), axis=1)
        result[start:stop] = np.take_along_axis(labels, order, axis=1)
        result_scores[start:stop] = np.take_along_axis(scores, order, axis=1)
    return result, result_scores


def threshold_stats(
    exact_labels: np.ndarray,
    exact_scores: np.ndarray,
    approximate_labels: np.ndarray,
    approximate_scores: np.ndarray,
    thresholds: list[float],
) -> list[dict[str, object]]:
    """Compare thresholded top-k sets in the same per-position score space."""
    if exact_labels.shape != approximate_labels.shape:
        raise ValueError("exact and approximate label arrays must have the same shape")
    if exact_scores.shape != approximate_scores.shape or exact_scores.shape != exact_labels.shape:
        raise ValueError("threshold score arrays must match the label arrays")
    rows: list[dict[str, object]] = []
    for threshold in thresholds:
        exact_mask = exact_scores >= threshold
        approximate_mask = approximate_scores >= threshold
        exact_count = int(exact_mask.sum())
        approximate_count = int(approximate_mask.sum())
        overlap_count = 0
        query_recalls: list[float] = []
        query_precisions: list[float] = []
        for row in range(exact_labels.shape[0]):
            exact_set = set(exact_labels[row, exact_mask[row]].tolist())
            approximate_set = set(approximate_labels[row, approximate_mask[row]].tolist())
            overlap = len(exact_set & approximate_set)
            overlap_count += overlap
            if exact_set:
                query_recalls.append(overlap / float(len(exact_set)))
            if approximate_set:
                query_precisions.append(overlap / float(len(approximate_set)))
        rows.append(
            {
                "threshold": float(threshold),
                "exact_hits": exact_count,
                "approximate_hits": approximate_count,
                "overlap_hits": overlap_count,
                "exact_hits_per_query": exact_count / float(exact_labels.shape[0]),
                "approximate_hits_per_query": approximate_count / float(approximate_labels.shape[0]),
                "micro_recall": overlap_count / float(exact_count) if exact_count else None,
                "micro_precision": overlap_count / float(approximate_count) if approximate_count else None,
                "mean_query_recall": float(np.mean(query_recalls)) if query_recalls else None,
                "mean_query_precision": float(np.mean(query_precisions)) if query_precisions else None,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> int:
    ks = parse_int_list(args.k_values)
    ms = parse_int_list(args.m_values)
    efs = parse_int_list(args.ef_search_values)
    thresholds = parse_float_list(args.thresholds)
    if any(value < -1.0 or value > 1.0 for value in thresholds):
        raise ValueError("thresholds must be in [-1, 1] per-position cosine score units")
    records, node_count, window_count = load_structures(args.structures, args.window_size, args.stride)
    queries, query_stats = load_queries(args.queries, records, args.max_queries)
    embeddings = source_embeddings(args.embeddings, node_count, args.original_dtype)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "benchmark.json").write_text(
        json.dumps(
            {
                "backend": "hnswlib_cpp_compact",
                "structures": str(args.structures),
                "queries": str(args.queries),
                "embeddings": str(args.embeddings),
                "original_dtype": args.original_dtype,
                "window_size": args.window_size,
                "stride": args.stride,
                "record_count": len(records),
                "node_count": node_count,
                "window_count": window_count,
                "query_stats": query_stats,
                "thresholds": thresholds,
                "candidate_k_values": (
                    parse_int_list(args.candidate_k_values)
                    if args.candidate_k_values
                    else [args.candidate_k or args.top_k]
                ),
            },
            indent=2,
        )
        + "\n"
    )

    flat_windows_path = args.outdir / f"flat_windows_{args.original_dtype}_w{args.window_size}_s{args.stride}.npy"
    windows = build_flat_windows(
        embeddings, records, args.window_size, args.stride, args.original_dtype, flat_windows_path
    )
    flat_index, flat_build_seconds = build_flat_index(windows, args.outdir / "flatip.index", args.batch_size)
    query_ids = np.asarray([query.window_id for query in queries], dtype=np.int64)
    query_vectors = np.ascontiguousarray(windows[query_ids], dtype=np.float32)
    exact_started = time.perf_counter()
    exact_distances, exact_labels = flat_index.search(query_vectors, args.top_k)
    flat_query_seconds = time.perf_counter() - exact_started
    (args.outdir / "flatip.json").write_text(
        json.dumps(
            {
                "backend": "FlatIP",
                "n_queries": len(queries),
                "n_windows": window_count,
                "top_k": args.top_k,
                "build_seconds": flat_build_seconds,
                "query_seconds": flat_query_seconds,
                "index_bytes": file_size(args.outdir / "flatip.index"),
                "self_hit_at_1": self_hit_rate(exact_labels, query_ids, 1),
                "self_hit_at_5": self_hit_rate(exact_labels, query_ids, min(5, args.top_k)),
                "self_hit_at_10": self_hit_rate(exact_labels, query_ids, min(10, args.top_k)),
                "self_hit_at_50": self_hit_rate(exact_labels, query_ids, min(50, args.top_k)),
            },
            indent=2,
        )
        + "\n"
    )

    rows: list[dict[str, object]] = []
    threshold_rows: list[dict[str, object]] = []
    candidate_ks = (
        parse_int_list(args.candidate_k_values)
        if args.candidate_k_values
        else [args.candidate_k or args.top_k]
    )
    if any(value < args.top_k for value in candidate_ks):
        raise ValueError("every candidate-k must be >= top-k")
    for k in ks:
        centroids = load_or_fit_centroids(
            embeddings,
            k,
            args.centroids_dir,
            args.sample_size,
            args.niter,
            args.seed,
            args.outdir / f"centroids_k{k}.float16.npy",
        )
        codes = build_codes(embeddings, centroids, args.outdir / f"codes_k{k}.npy", args.batch_size)
        code_path = args.outdir / f"window_codes_k{k}.bin"
        write_raw_window_codes(codes, records, args.window_size, args.stride, code_path)
        query_code_path = args.outdir / f"query_codes_k{k}.bin"
        query_code_windows(codes, queries, records, args.window_size).tofile(query_code_path)
        similarity_path = args.outdir / f"similarity_k{k}.f32.bin"
        if not similarity_path.exists() or similarity_path.stat().st_size != k * k * np.dtype(np.float32).itemsize:
            write_similarity(centroids, similarity_path)
        for m in ms:
            index_path = args.outdir / f"compact_hnsw_k{k}_m{m}_efc{args.ef_construction}.bin"
            build_started = time.perf_counter()
            if not index_path.exists():
                command = [
                    str(args.executable),
                    "build",
                    "--codes", str(code_path),
                    "--similarity", str(similarity_path),
                    "--count", str(window_count),
                    "--window-size", str(args.window_size),
                    "--n-centroids", str(k),
                    "--index", str(index_path),
                    "--m", str(m),
                    "--ef-construction", str(args.ef_construction),
                    "--ef-search", str(max(efs)),
                    "--random-seed", str(args.seed),
                    "--threads", str(args.num_threads),
                ]
                subprocess.run(command, check=True)
            build_seconds = time.perf_counter() - build_started
            for candidate_k in candidate_ks:
                for ef_search in efs:
                    labels, distances, query_seconds = run_compact_search(
                        args.executable,
                        index_path,
                        similarity_path,
                        query_code_path,
                        len(queries),
                        args.window_size,
                        k,
                        candidate_k,
                        ef_search,
                        args.num_threads,
                        args.outdir,
                    )
                    rerank_started = time.perf_counter()
                    if args.rerank:
                        final_labels, final_scores = rerank_candidates(
                            windows,
                            query_ids,
                            labels,
                            args.top_k,
                            args.rerank_batch_size,
                        )
                        threshold_score_space = "original_window_cosine"
                    else:
                        final_labels = labels[:, : args.top_k]
                        final_scores = (-distances[:, : args.top_k]) / float(args.window_size)
                        threshold_score_space = "centroid_score_per_position"
                    rerank_seconds = time.perf_counter() - rerank_started if args.rerank else 0.0
                    threshold_rows.append(
                        {
                            "k": k,
                            "m": m,
                            "ef_construction": args.ef_construction,
                            "ef_search": ef_search,
                            "candidate_k": candidate_k,
                            "rerank": bool(args.rerank),
                            "score_space": threshold_score_space,
                            "rows": threshold_stats(
                                exact_labels,
                                exact_distances[:, : args.top_k],
                                final_labels,
                                final_scores,
                                thresholds,
                            ),
                        }
                    )
                    fields: dict[str, object] = {
                        "backend": "hnswlib_cpp_compact",
                        "k": k,
                        "m": m,
                        "ef_construction": args.ef_construction,
                        "ef_search": ef_search,
                        "candidate_k": candidate_k,
                        "rerank": bool(args.rerank),
                        "n_queries": len(queries),
                        "n_windows": window_count,
                        "window_size": args.window_size,
                        "code_bytes_per_window": args.window_size * np.dtype(np.uint16).itemsize,
                        "flat_build_seconds": flat_build_seconds,
                        "flat_query_seconds": flat_query_seconds,
                        "hnsw_build_seconds": build_seconds,
                        "hnsw_query_seconds": query_seconds,
                        "rerank_seconds": rerank_seconds,
                        "total_approx_query_seconds": query_seconds + rerank_seconds,
                        "hnsw_index_bytes": file_size(index_path),
                        "candidate_recall_at_1": recall_at(exact_labels, labels, 1),
                        "candidate_recall_at_5": recall_at(exact_labels, labels, min(5, args.top_k)),
                        "candidate_recall_at_10": recall_at(exact_labels, labels, min(10, args.top_k)),
                        "candidate_recall_at_50": recall_at(exact_labels, labels, min(50, args.top_k)),
                        "recall_at_1": recall_at(exact_labels, final_labels, 1),
                        "recall_at_5": recall_at(exact_labels, final_labels, min(5, args.top_k)),
                        "recall_at_10": recall_at(exact_labels, final_labels, min(10, args.top_k)),
                        "recall_at_50": recall_at(exact_labels, final_labels, min(50, args.top_k)),
                        "candidate_self_hit_at_1": self_hit_rate(labels, query_ids, 1),
                        "self_hit_at_1": self_hit_rate(final_labels, query_ids, 1),
                        "self_hit_at_5": self_hit_rate(final_labels, query_ids, min(5, args.top_k)),
                        "self_hit_at_10": self_hit_rate(final_labels, query_ids, min(10, args.top_k)),
                        "self_hit_at_50": self_hit_rate(final_labels, query_ids, min(50, args.top_k)),
                        "mean_self_rank_at_50": positive_rank(final_labels, query_ids, min(50, args.top_k)),
                        "mean_exact_top1_score": float(np.mean(exact_distances[:, 0])),
                        "mean_compact_top1_score": float(np.mean(-distances[:, 0])),
                        "threshold_score_space": threshold_score_space,
                    }
                    rows.append(fields)
                    log(
                        f"w={args.window_size} k={k} M={m} ef={ef_search} candidate={candidate_k}: "
                        f"candidate-R@50={fields['candidate_recall_at_50']:.4f} "
                        f"R@50={fields['recall_at_50']:.4f} "
                        f"self@1={fields['self_hit_at_1']:.4f} query={query_seconds:.3f}s "
                        f"index={fields['hnsw_index_bytes']}"
                    )
    write_csv(args.outdir / "results.csv", rows)
    (args.outdir / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
    (args.outdir / "thresholds.json").write_text(json.dumps(threshold_rows, indent=2) + "\n")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--centroids-dir", type=Path)
    parser.add_argument("--k-values", default="2048")
    parser.add_argument("--m-values", default="32")
    parser.add_argument("--ef-search-values", default="100")
    parser.add_argument("--ef-construction", type=int, default=200)
    parser.add_argument("--num-threads", type=int, default=0)
    parser.add_argument("--sample-size", type=int, default=500_000)
    parser.add_argument("--niter", type=int, default=25)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32_768)
    parser.add_argument("--top-k", type=int, default=50)
    parser.add_argument("--candidate-k", type=int)
    parser.add_argument("--candidate-k-values")
    parser.add_argument("--rerank", action="store_true")
    parser.add_argument("--rerank-batch-size", type=int, default=32)
    parser.add_argument(
        "--thresholds",
        default="0.7,0.75,0.8,0.85,0.9,0.95",
        help="Comma-separated per-position score thresholds for hit-count and overlap analysis.",
    )
    parser.add_argument("--max-queries", type=int)
    parser.add_argument("--original-dtype", choices=("float16", "float32"), default="float16")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.window_size < 1 or args.stride < 1 or args.batch_size < 1 or args.top_k < 1 or args.rerank_batch_size < 1:
            raise ValueError("window size, stride, batch size, top-k, and rerank batch size must be >= 1")
        if args.ef_construction < 1 or args.num_threads < 0:
            raise ValueError("ef-construction must be >= 1 and num-threads must be >= 0")
        return run(args)
    except (OSError, KeyError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
