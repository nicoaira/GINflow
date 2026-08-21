#!/usr/bin/env python3
"""Research benchmark for a GPU CAGRA companion index.

The script is intentionally independent of Nextflow. It uploads a NumPy
window matrix in batches, builds or loads a CAGRA graph, searches a query
matrix, and optionally compares the labels with an exact reference saved as a
NumPy array. It remains a lightweight research harness for reproducing GPU
parameter sweeps; production HNSWLIB search is implemented in
``hnswlib_gpu.py`` and performs the final original-embedding rerank.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from rerank_core import rerank_label_matrix


def convert(values: np.ndarray, dtype: str, int8_scale: float) -> np.ndarray:
    if dtype == "int8":
        return np.clip(np.rint(np.asarray(values, dtype=np.float32) * int8_scale), -128, 127).astype(
            np.int8
        )
    return np.asarray(values, dtype=np.dtype(dtype))


def parse_int_list(value: str) -> list[int]:
    values = [int(item.strip()) for item in str(value).split(",") if item.strip()]
    if not values or any(item < 1 for item in values):
        raise ValueError(f"expected positive comma-separated integers, got {value!r}")
    return values


def recall_at(reference: np.ndarray, approximate: np.ndarray, k: int) -> float | None:
    if k < 1 or reference.shape[1] < k or approximate.shape[1] < k:
        return None
    overlaps = [
        len(set(expected.tolist()) & set(actual.tolist())) / float(k)
        for expected, actual in zip(reference[:, :k], approximate[:, :k])
    ]
    return float(np.mean(overlaps)) if overlaps else None


def upload(path: Path, dtype: str, batch_size: int, int8_scale: float):
    import cupy as cp

    source = np.load(path, mmap_mode="r")
    if source.ndim != 2:
        raise ValueError(f"{path} must be a 2-D NumPy array, got {source.shape}")
    target_dtype = np.dtype(dtype)
    device = cp.empty(source.shape, dtype=target_dtype)
    started = time.perf_counter()
    for start in range(0, source.shape[0], batch_size):
        stop = min(start + batch_size, source.shape[0])
        device[start:stop] = cp.asarray(convert(source[start:stop], dtype, int8_scale))
    cp.cuda.Stream.null.synchronize()
    return device, time.perf_counter() - started


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--index")
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--dtype", choices=("float16", "float32", "int8"), default="float16")
    parser.add_argument("--int8-scale", type=float, default=127.0)
    parser.add_argument("--batch-size", type=int, default=32768)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--itopk-size", type=int, default=64)
    parser.add_argument("--intermediate-graph-degree", type=int, default=128)
    parser.add_argument("--graph-degree", type=int, default=64)
    parser.add_argument("--build-algo", default="nn_descent")
    parser.add_argument("--metric", choices=("sqeuclidean", "inner_product"), default="sqeuclidean")
    parser.add_argument("--max-queries", type=int)
    parser.add_argument("--search-batch-size", type=int, default=512)
    parser.add_argument("--output-k", type=int)
    parser.add_argument("--rerank-batch-size", type=int, default=8)
    parser.add_argument("--candidate-batch-size", type=int, default=2048)
    parser.add_argument("--rerank-workers", type=int, default=1)
    parser.add_argument("--rerank-device", choices=("cpu", "cuda", "auto"), default="cuda")
    parser.add_argument("--recall-k-values", default="1,5,10,50,100,500")
    parser.add_argument("--skip-save", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1 or args.search_batch_size < 1 or args.k < 1 or args.itopk_size < args.k:
        raise ValueError("batch sizes and k must be positive; itopk-size must be >= k")

    import cupy as cp
    from cuvs.neighbors import cagra

    args.outdir.mkdir(parents=True, exist_ok=True)
    source = np.load(args.data, mmap_mode="r")
    if source.ndim != 2:
        raise ValueError(f"{args.data} must be a 2-D NumPy array, got {source.shape}")
    index_path = Path(args.index) if args.index else args.outdir / "cagra.index"
    upload_seconds = 0.0
    data = None
    if index_path.exists():
        # A serialized CAGRA index already owns its dataset.  Uploading a
        # second copy can consume the memory needed by a large search beam.
        started = time.perf_counter()
        index = cagra.load(str(index_path))
        cp.cuda.Stream.null.synchronize()
        load_seconds = time.perf_counter() - started
        dimension = int(index.dim)
    else:
        data, upload_seconds = upload(args.data, args.dtype, args.batch_size, args.int8_scale)
        dimension = int(data.shape[1])
        index = None
    query_source = np.load(args.queries, mmap_mode="r")
    if query_source.ndim != 2 or query_source.shape[1] != dimension:
        raise ValueError(f"query shape {query_source.shape} is incompatible with data dimension {dimension}")
    if args.max_queries is not None:
        query_source = query_source[: args.max_queries]
    # CAGRA requires query and dataset dtypes to match for this release.
    queries = cp.asarray(convert(query_source, args.dtype, args.int8_scale))

    build_seconds = 0.0
    if index is None:
        started = time.perf_counter()
        params = cagra.IndexParams(
            metric=args.metric,
            intermediate_graph_degree=args.intermediate_graph_degree,
            graph_degree=args.graph_degree,
            build_algo=args.build_algo,
        )
        index = cagra.build(params, data)
        cp.cuda.Stream.null.synchronize()
        build_seconds = time.perf_counter() - started
        if not args.skip_save:
            cagra.save(str(index_path), index, include_dataset=True)
        load_seconds = 0.0

    search_params = cagra.SearchParams(
        itopk_size=args.itopk_size,
        max_queries=min(int(queries.shape[0]), args.search_batch_size),
    )
    # Warm up kernels and measure synchronized device time.
    cagra.search(search_params, index, queries[:1], args.k)
    cp.cuda.Stream.null.synchronize()
    started = time.perf_counter()
    distance_parts = []
    label_parts = []
    for start in range(0, int(queries.shape[0]), args.search_batch_size):
        stop = min(start + args.search_batch_size, int(queries.shape[0]))
        distances, labels = cagra.search(search_params, index, queries[start:stop], args.k)
        cp.cuda.Stream.null.synchronize()
        distance_parts.append(cp.asnumpy(cp.asarray(distances)))
        label_parts.append(cp.asnumpy(cp.asarray(labels)))
    distances = np.concatenate(distance_parts, axis=0)
    labels = np.concatenate(label_parts, axis=0)
    search_seconds = time.perf_counter() - started
    labels = np.asarray(labels, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.float32)
    np.save(args.outdir / "labels.npy", labels)
    np.save(args.outdir / "distances.npy", distances)

    exact_queries = np.ascontiguousarray(query_source, dtype=np.float32)
    exact_rerank_seconds = 0.0
    if args.output_k is not None:
        if args.output_k < 1 or args.output_k > labels.shape[1]:
            raise ValueError("output-k must be between 1 and the CAGRA candidate k")
        final_labels, final_scores, exact_rerank_seconds = rerank_label_matrix(
            source,
            exact_queries,
            labels,
            args.output_k,
            batch_size=args.rerank_batch_size,
            candidate_batch_size=args.candidate_batch_size,
            workers=args.rerank_workers,
            device=args.rerank_device,
        )
    else:
        final_labels = labels
        final_scores = -distances
    np.save(args.outdir / "final_labels.npy", final_labels)
    np.save(args.outdir / "final_scores.npy", final_scores)

    result = {
        "backend": "cuvs_cagra",
        "data": str(args.data),
        "queries": str(args.queries),
        "data_shape": list(source.shape),
        "data_dtype": args.dtype,
        "int8_scale": args.int8_scale if args.dtype == "int8" else None,
        "query_count": int(queries.shape[0]),
        "k": int(args.k),
        "itopk_size": int(args.itopk_size),
        "search_batch_size": int(args.search_batch_size),
        "intermediate_graph_degree": int(args.intermediate_graph_degree),
        "graph_degree": int(args.graph_degree),
        "build_algo": args.build_algo,
        "metric": args.metric,
        "upload_seconds": upload_seconds,
        "build_seconds": build_seconds,
        "load_seconds": load_seconds,
        "search_seconds": search_seconds,
        "exact_rerank_seconds": exact_rerank_seconds,
        "search_plus_rerank_seconds": search_seconds + exact_rerank_seconds,
        "output_k": int(final_labels.shape[1]),
        "rerank_device": args.rerank_device if args.output_k is not None else None,
    }
    if args.reference:
        reference = np.asarray(np.load(args.reference), dtype=np.int64)
        if reference.shape[0] != labels.shape[0]:
            raise ValueError(f"reference rows {reference.shape[0]} != result rows {labels.shape[0]}")
        for k in parse_int_list(args.recall_k_values):
            if k > labels.shape[1] or k > reference.shape[1]:
                continue
            result[f"candidate_recall_at_{k}"] = recall_at(reference, labels, k)
            result[f"final_recall_at_{k}"] = recall_at(reference, final_labels, k)
    (args.outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
