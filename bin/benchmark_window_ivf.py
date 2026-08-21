#!/usr/bin/env python3
"""Benchmark the standalone custom-distance WindowIVF prototype.

The node embeddings are trained/encoded with product quantization first.  The
window index stores only packed node-PQ codes and uses the positional lookup
table in ``bin/window_ivf.cpp`` for coarse routing and list scanning.  Exact
reranking always reads the preserved original window vectors.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from benchmark_faiss_hnsw_int8 import recall_at
from benchmark_hnswlib import load_queries, load_structures, source_embeddings
from benchmark_pq_hnswlib import (
    OriginalNodeWindowStore,
    pack_window_codes,
    query_window_codes,
    run_timed,
    train_and_encode_nodes,
)
from rerank_core import rerank_label_matrix


def parse_int_list(value: str) -> list[int]:
    values = [int(item.strip()) for item in str(value).split(",") if item.strip()]
    if not values or any(item < 1 for item in values):
        raise ValueError(f"expected positive comma-separated integers, got {value!r}")
    return values


def run_ivf_build(
    executable: Path,
    codes_path: Path,
    similarity_path: Path,
    index_path: Path,
    count: int,
    window_size: int,
    pq_m: int,
    nbits: int,
    nlist: int,
    niter: int,
    threads: int,
    output_dir: Path,
) -> tuple[float, float]:
    signature = {
        "count": count,
        "window_size": window_size,
        "pq_m": pq_m,
        "nbits": nbits,
        "nlist": nlist,
        "niter": niter,
        "threads": threads,
        "codes_bytes": codes_path.stat().st_size,
        "similarity_bytes": similarity_path.stat().st_size,
    }
    metadata_path = index_path.with_suffix(".build.json")
    metrics_path = index_path.with_suffix(".build.metrics.json")
    if index_path.exists() and metadata_path.exists() and metrics_path.exists():
        try:
            if json.loads(metadata_path.read_text()) == signature:
                metrics = json.loads(metrics_path.read_text())
                return float(metrics["elapsed_seconds"]), float(metrics["peak_rss_mib"])
        except (OSError, TypeError, ValueError, KeyError):
            pass
    command = [
        str(executable),
        "build",
        "--codes", str(codes_path),
        "--similarity", str(similarity_path),
        "--count", str(count),
        "--window-size", str(window_size),
        "--pq-m", str(pq_m),
        "--nbits", str(nbits),
        "--nlist", str(nlist),
        "--niter", str(niter),
        "--threads", str(threads),
        "--index", str(index_path),
    ]
    elapsed, peak_rss_mib = run_timed(command, output_dir, f"ivf-build-nlist{nlist}")
    metadata_path.write_text(json.dumps(signature, indent=2) + "\n")
    metrics_path.write_text(json.dumps({
        "elapsed_seconds": elapsed,
        "peak_rss_mib": peak_rss_mib,
    }, indent=2) + "\n")
    return elapsed, peak_rss_mib


def run_ivf_search(
    executable: Path,
    index_path: Path,
    query_codes_path: Path,
    similarity_path: Path,
    query_count: int,
    window_size: int,
    pq_m: int,
    nbits: int,
    nprobe: int,
    candidate_k: int,
    threads: int,
    output_dir: Path,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    stem = f"ivf-search-nlist{json.loads(index_path.with_suffix('.build.json').read_text())['nlist']}-probe{nprobe}-k{candidate_k}"
    labels_path = output_dir / f"{stem}.labels.u64"
    distances_path = output_dir / f"{stem}.distances.f32"
    metadata_path = output_dir / f"{stem}.json"
    metrics_path = output_dir / f"{stem}.metrics.json"
    signature = {
        "query_count": query_count,
        "window_size": window_size,
        "pq_m": pq_m,
        "nbits": nbits,
        "nprobe": nprobe,
        "candidate_k": candidate_k,
        "threads": threads,
        "index_bytes": index_path.stat().st_size,
    }
    expected_labels = query_count * candidate_k * np.dtype("<u8").itemsize
    expected_distances = query_count * candidate_k * np.dtype("<f4").itemsize
    if labels_path.exists() and distances_path.exists() and metadata_path.exists() and metrics_path.exists():
        try:
            if (
                json.loads(metadata_path.read_text()) == signature
                and labels_path.stat().st_size == expected_labels
                and distances_path.stat().st_size == expected_distances
            ):
                metrics = json.loads(metrics_path.read_text())
                return (
                    np.fromfile(labels_path, dtype="<u8").reshape(query_count, candidate_k),
                    np.fromfile(distances_path, dtype="<f4").reshape(query_count, candidate_k),
                    float(metrics["elapsed_seconds"]),
                    float(metrics["peak_rss_mib"]),
                )
        except (OSError, TypeError, ValueError, KeyError):
            pass
    command = [
        str(executable),
        "search",
        "--codes", str(query_codes_path),
        "--similarity", str(similarity_path),
        "--query-count", str(query_count),
        "--window-size", str(window_size),
        "--pq-m", str(pq_m),
        "--nbits", str(nbits),
        "--index", str(index_path),
        "--nprobe", str(nprobe),
        "--k", str(candidate_k),
        "--threads", str(threads),
        "--labels-out", str(labels_path),
        "--distances-out", str(distances_path),
    ]
    elapsed, peak_rss_mib = run_timed(command, output_dir, stem)
    metadata_path.write_text(json.dumps(signature, indent=2) + "\n")
    metrics_path.write_text(json.dumps({
        "elapsed_seconds": elapsed,
        "peak_rss_mib": peak_rss_mib,
    }, indent=2) + "\n")
    return (
        np.fromfile(labels_path, dtype="<u8").reshape(query_count, candidate_k),
        np.fromfile(distances_path, dtype="<f4").reshape(query_count, candidate_k),
        elapsed,
        peak_rss_mib,
    )


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    pq_m_values = parse_int_list(args.pq_m_values)
    nbits_values = parse_int_list(args.nbits_values)
    nlist_values = parse_int_list(args.nlist_values)
    nprobe_values = parse_int_list(args.nprobe_values)
    candidate_values = parse_int_list(args.candidate_k_values)
    recall_values = parse_int_list(args.recall_k_values)
    if any(bits > 8 for bits in nbits_values):
        raise ValueError("nbits must be <= 8")
    if any(probe > nlist for probe in nprobe_values for nlist in nlist_values):
        raise ValueError("every nprobe value must be <= every nlist value")
    if args.output_k > min(candidate_values):
        raise ValueError("output-k must not exceed the smallest candidate-k")

    records, node_count, window_count = load_structures(args.structures, args.window_size, args.stride)
    queries, query_stats = load_queries(args.query_selections, records, args.max_queries)
    embeddings = source_embeddings(args.embeddings, node_count, args.original_dtype)
    database_windows = np.load(args.database_windows, mmap_mode="r")
    query_windows = np.asarray(np.load(args.queries), dtype=np.float32)[: len(queries)]
    reference = np.asarray(np.load(args.reference_labels), dtype=np.int64)[: len(queries)]
    if database_windows.shape[0] != window_count:
        raise ValueError("database windows and structures have different window counts")
    if query_windows.shape != (len(queries), database_windows.shape[1]):
        raise ValueError("query windows and query selections are incompatible")
    original_windows = OriginalNodeWindowStore(embeddings, records, args.window_size, args.stride)
    if reference.shape[0] != len(queries):
        raise ValueError("reference labels and query selections are incompatible")

    args.outdir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for pq_m in pq_m_values:
        if embeddings.shape[1] % pq_m:
            continue
        for nbits in nbits_values:
            config_dir = args.outdir / f"pq_m{pq_m}_b{nbits}"
            config_dir.mkdir(parents=True, exist_ok=True)
            model = train_and_encode_nodes(
                embeddings,
                pq_m,
                nbits,
                args.sample_size,
                args.pq_niter,
                args.seed,
                args.batch_size,
                args.original_dtype,
                config_dir,
            )
            node_codes = np.load(model["node_code_path"], mmap_mode="r")
            database_codes_path = config_dir / "database_windows.codes.bin"
            query_codes_path = config_dir / "query_windows.codes.bin"
            pack_started = time.perf_counter()
            actual_windows = pack_window_codes(
                node_codes,
                records,
                args.window_size,
                args.stride,
                pq_m,
                nbits,
                database_codes_path,
            )
            query_window_codes(
                node_codes,
                queries,
                records,
                args.window_size,
                args.stride,
                pq_m,
                nbits,
            ).tofile(query_codes_path)
            pack_seconds = time.perf_counter() - pack_started
            if actual_windows != window_count:
                raise ValueError("packed window count does not match structures")
            for nlist in nlist_values:
                index_path = config_dir / f"window_ivf_nlist{nlist}.bin"
                build_seconds, build_peak_rss_mib = run_ivf_build(
                    args.executable,
                    database_codes_path,
                    Path(model["similarity_path"]),
                    index_path,
                    window_count,
                    args.window_size,
                    pq_m,
                    nbits,
                    nlist,
                    args.niter,
                    args.threads,
                    config_dir,
                )
                for nprobe in nprobe_values:
                    for candidate_k in candidate_values:
                        labels, distances, search_seconds, search_peak_rss_mib = run_ivf_search(
                            args.executable,
                            index_path,
                            query_codes_path,
                            Path(model["similarity_path"]),
                            len(queries),
                            args.window_size,
                            pq_m,
                            nbits,
                            nprobe,
                            candidate_k,
                            args.threads,
                            config_dir,
                        )
                        final_labels, final_scores, rerank_seconds = rerank_label_matrix(
                            original_windows,
                            query_windows,
                            labels,
                            args.output_k,
                            batch_size=args.rerank_batch_size,
                            candidate_batch_size=args.candidate_batch_size,
                            workers=args.rerank_workers,
                            device=args.rerank_device,
                        )
                        row: dict[str, Any] = {
                            "backend": "window_ivf_cpp_node_pq",
                            "pq_m": pq_m,
                            "nbits": nbits,
                            "nlist": nlist,
                            "nprobe": nprobe,
                            "candidate_k": candidate_k,
                            "output_k": args.output_k,
                            "n_windows": window_count,
                            "n_queries": len(queries),
                            "node_code_bytes": int(model["node_code_bytes"]),
                            "window_code_bytes": int(args.window_size * model["node_code_bytes"]),
                            "node_code_storage_bytes": int(model["node_code_storage_bytes"]),
                            "codebook_bytes": int(model["codebook_bytes"]),
                            "similarity_bytes": int(model["similarity_bytes"]),
                            "pq_train_seconds": model["train_seconds"],
                            "pq_encode_seconds": model["encode_seconds"],
                            "window_pack_seconds": pack_seconds,
                            "build_seconds": build_seconds,
                            "build_peak_rss_mib": build_peak_rss_mib,
                            "search_seconds": search_seconds,
                            "search_peak_rss_mib": search_peak_rss_mib,
                            "exact_rerank_seconds": rerank_seconds,
                            "search_plus_rerank_seconds": search_seconds + rerank_seconds,
                            "rerank_workers": args.rerank_workers,
                            "rerank_device": args.rerank_device,
                            "index_bytes": index_path.stat().st_size,
                            "candidate_mean_score": float(np.mean(-distances[:, 0])),
                            "final_mean_score": float(np.mean(final_scores[:, 0])),
                        }
                        for recall_k in recall_values:
                            if recall_k <= labels.shape[1] and recall_k <= reference.shape[1]:
                                row[f"candidate_recall_at_{recall_k}"] = recall_at(reference, labels, recall_k)
                            if recall_k <= final_labels.shape[1] and recall_k <= reference.shape[1]:
                                row[f"final_recall_at_{recall_k}"] = recall_at(reference, final_labels, recall_k)
                        results.append(row)
                        print(json.dumps(row), flush=True)
    benchmark = {
        "backend": "window_ivf_cpp_node_pq",
        "embeddings": str(args.embeddings),
        "structures": str(args.structures),
        "query_selections": str(args.query_selections),
        "database_windows": str(args.database_windows),
        "queries": str(args.queries),
        "reference_labels": str(args.reference_labels),
        "original_dtype": args.original_dtype,
        "window_size": args.window_size,
        "stride": args.stride,
        "node_count": node_count,
        "n_windows": window_count,
        "n_queries": len(queries),
        "query_stats": query_stats,
        "results": results,
    }
    (args.outdir / args.results_name).write_text(json.dumps(benchmark, indent=2) + "\n")
    return benchmark


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--query-selections", type=Path, required=True)
    parser.add_argument("--database-windows", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--reference-labels", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--results-name", default="results.json")
    parser.add_argument("--pq-m-values", default="16")
    parser.add_argument("--nbits-values", default="4,8")
    parser.add_argument("--nlist-values", default="256,1024")
    parser.add_argument("--nprobe-values", default="8,32,128")
    parser.add_argument("--candidate-k-values", default="1000,5000,10000")
    parser.add_argument("--output-k", type=int, default=500)
    parser.add_argument("--recall-k-values", default="100,500")
    parser.add_argument("--sample-size", type=int, default=500_000)
    parser.add_argument("--pq-niter", type=int, default=25)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32768)
    parser.add_argument("--rerank-batch-size", type=int, default=8)
    parser.add_argument("--candidate-batch-size", type=int, default=1024)
    parser.add_argument("--rerank-workers", type=int, default=8)
    parser.add_argument("--rerank-device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--niter", type=int, default=1)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--original-dtype", choices=("float16", "float32"), default="float16")
    parser.add_argument("--max-queries", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(evaluate(parse_args()), indent=2))
