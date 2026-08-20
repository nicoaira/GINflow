#!/usr/bin/env python3
"""Benchmark quantized-window HNSWLIB candidate recall against exact FlatIP.

The benchmark is intentionally independent of Nextflow.  It consumes a flat
node-embedding cache and a structures TSV, builds the exact normalized window
vectors used by the existing FlatIP backend, and then evaluates one or more
centroid codebooks with HNSWLIB.  Query rows are resolved by ``transcript_id``
and ``window_offset``; the ``vector_id`` and ``record_ordinal`` columns are
retained as provenance because those fields can refer to a different source
ordering.

All generated arrays and indexes are placed below ``--outdir``.  This makes a
large run resumable and keeps research artifacts separate from published
pipeline output.  The intended location for the Rouskin run is
``/mnt/ssd_samsung/ginflow-hnsw-research``.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np


@dataclass(frozen=True)
class Record:
    identifier: str
    length: int
    node_start: int
    window_start: int
    n_windows: int


@dataclass(frozen=True)
class Query:
    query_ordinal: int
    vector_id: int
    record_ordinal: int
    transcript_id: str
    window_offset: int
    selection_cycle: int
    window_id: int


def log(message: str) -> None:
    print(message, flush=True)


def parse_int_list(value: str) -> list[int]:
    values = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not values or any(item < 1 for item in values):
        raise ValueError(f"expected a comma-separated list of positive integers, got {value!r}")
    return values


def load_structures(path: Path, window_size: int, stride: int) -> tuple[list[Record], int, int]:
    records: list[Record] = []
    seen: set[str] = set()
    node_start = 0
    window_start = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "sequence"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must have columns {sorted(required)}")
        for row_number, row in enumerate(reader, start=2):
            identifier = str(row["transcript_id"])
            sequence = str(row["sequence"])
            if not identifier or identifier in seen:
                raise ValueError(f"duplicate or empty transcript_id at {path}:{row_number}: {identifier!r}")
            if not sequence:
                raise ValueError(f"empty sequence at {path}:{row_number}")
            seen.add(identifier)
            length = len(sequence)
            n_windows = max(0, (length - window_size) // stride + 1)
            records.append(Record(identifier, length, node_start, window_start, n_windows))
            node_start += length
            window_start += n_windows
    if not records or node_start < 1 or window_start < 1:
        raise ValueError(f"{path} contains no usable records/windows")
    return records, node_start, window_start


def load_queries(path: Path, records: list[Record], max_queries: int | None) -> tuple[list[Query], dict[str, int]]:
    by_id = {record.identifier: record for record in records}
    queries: list[Query] = []
    stats = {
        "rows": 0,
        "resolved": 0,
        "record_ordinal_matches": 0,
        "vector_id_matches": 0,
    }
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "query_ordinal",
            "vector_id",
            "record_ordinal",
            "transcript_id",
            "window_offset",
            "selection_cycle",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must have columns {sorted(required)}")
        for row in reader:
            stats["rows"] += 1
            identifier = str(row["transcript_id"])
            record = by_id.get(identifier)
            if record is None:
                continue
            offset = int(row["window_offset"])
            if offset < 0 or offset >= record.n_windows:
                continue
            actual_window_id = record.window_start + offset
            stats["resolved"] += 1
            if int(row["record_ordinal"]) == records.index(record):
                stats["record_ordinal_matches"] += 1
            if int(row["vector_id"]) == actual_window_id:
                stats["vector_id_matches"] += 1
            queries.append(
                Query(
                    query_ordinal=int(row["query_ordinal"]),
                    vector_id=int(row["vector_id"]),
                    record_ordinal=int(row["record_ordinal"]),
                    transcript_id=identifier,
                    window_offset=offset,
                    selection_cycle=int(row["selection_cycle"]),
                    window_id=actual_window_id,
                )
            )
            if max_queries is not None and len(queries) >= max_queries:
                break
    if not queries:
        raise ValueError(f"no query rows from {path} could be resolved against {len(records)} records")
    return queries, stats


def source_embeddings(path: Path, expected_nodes: int, original_dtype: str) -> np.ndarray:
    values = np.load(path, mmap_mode="r")
    if values.ndim != 2 or values.shape[0] != expected_nodes:
        raise ValueError(
            f"{path} has shape {values.shape}; expected ({expected_nodes}, embedding_dim)"
        )
    if original_dtype == "float16" and values.dtype != np.float16:
        log(f"Casting source node embeddings {values.dtype} -> float16 for pipeline parity")
    return values


def normalize_rows(values: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype=np.float32)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    return array / np.maximum(norms, np.float32(1e-12))


def cache_matches(path: Path, shape: tuple[int, ...], dtype: np.dtype) -> bool:
    if not path.exists():
        return False
    try:
        values = np.load(path, mmap_mode="r")
    except (OSError, ValueError):
        return False
    return values.shape == shape and values.dtype == dtype


def build_flat_windows(
    embeddings: np.ndarray,
    records: list[Record],
    window_size: int,
    stride: int,
    original_dtype: str,
    output_path: Path,
) -> np.memmap:
    embedding_dim = int(embeddings.shape[1])
    shape = (records[-1].window_start + records[-1].n_windows, window_size * embedding_dim)
    dtype = np.dtype(np.float32)
    if not cache_matches(output_path, shape, dtype):
        output_path.parent.mkdir(parents=True, exist_ok=True)
        windows = np.lib.format.open_memmap(output_path, mode="w+", dtype=dtype, shape=shape)
        for record in records:
            if record.n_windows == 0:
                continue
            stop = record.node_start + record.length
            node_values = np.asarray(embeddings[record.node_start:stop])
            if original_dtype == "float16":
                node_values = np.asarray(node_values, dtype=np.float16)
            views = np.lib.stride_tricks.sliding_window_view(
                node_values, window_size, axis=0
            )[::stride]
            flat = np.ascontiguousarray(
                views.transpose(0, 2, 1).reshape(record.n_windows, window_size * embedding_dim),
                dtype=np.float32,
            )
            norms = np.linalg.norm(flat, axis=1, keepdims=True)
            windows[record.window_start : record.window_start + record.n_windows] = (
                flat / np.maximum(norms, np.float32(1e-12))
            )
        windows.flush()
        del windows
    return np.load(output_path, mmap_mode="r+")


def build_flat_index(windows: np.ndarray, output_path: Path, batch_size: int):
    import faiss

    if output_path.exists():
        log(f"Using cached FlatIP index: {output_path}")
        return faiss.read_index(str(output_path)), 0.0
    started = time.perf_counter()
    index = faiss.IndexFlatIP(int(windows.shape[1]))
    for start in range(0, windows.shape[0], batch_size):
        stop = min(start + batch_size, windows.shape[0])
        index.add(np.ascontiguousarray(windows[start:stop], dtype=np.float32))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    faiss.write_index(index, str(output_path))
    return index, time.perf_counter() - started


def find_centroid_file(root: Path, k: int) -> Path | None:
    if root.is_file():
        return root
    candidates = [
        root / f"centroids_k{k}.float16.npy",
        root / f"centroids_k{k}.float32.npy",
        root / f"centroids_k{k}.npy",
        root / f"spherical_kmeans_k{k}.npz",
        root / f"k{k}.npy",
        root / f"k{k}.npz",
    ]
    return next((path for path in candidates if path.exists()), None)


def load_centroids(path: Path, embedding_dim: int) -> np.ndarray:
    if path.suffix == ".npz":
        with np.load(path) as archive:
            if "centers" not in archive.files:
                raise ValueError(f"{path} must contain a 'centers' array")
            values = np.asarray(archive["centers"])
    else:
        values = np.asarray(np.load(path))
    if values.ndim != 2 or values.shape[1] != embedding_dim:
        raise ValueError(f"{path} has shape {values.shape}; expected (*, {embedding_dim})")
    if values.dtype != np.float16:
        values = normalize_rows(values)
        values = values.astype(np.float16)
    if values.shape[0] < 1:
        raise ValueError(f"{path} contains no centroids")
    return np.ascontiguousarray(values)


def fit_centroids(embeddings: np.ndarray, k: int, sample_size: int, niter: int, seed: int) -> np.ndarray:
    import faiss

    rng = np.random.default_rng(seed)
    sample_count = min(sample_size, embeddings.shape[0])
    indices = rng.choice(embeddings.shape[0], size=sample_count, replace=False)
    training = normalize_rows(np.asarray(embeddings[indices]))
    effective_k = min(k, sample_count)
    kmeans = faiss.Kmeans(
        int(training.shape[1]),
        int(effective_k),
        niter=int(niter),
        nredo=1,
        verbose=False,
        spherical=True,
        seed=int(seed),
    )
    kmeans.train(np.ascontiguousarray(training, dtype=np.float32))
    return np.ascontiguousarray(normalize_rows(kmeans.centroids), dtype=np.float16)


def load_or_fit_centroids(
    embeddings: np.ndarray,
    k: int,
    centroid_root: Path | None,
    sample_size: int,
    niter: int,
    seed: int,
    cache_path: Path,
) -> np.ndarray:
    if cache_path.exists():
        values = np.load(cache_path)
        if values.shape == (k, embeddings.shape[1]) and values.dtype == np.float16:
            return np.ascontiguousarray(values)
    source = find_centroid_file(centroid_root, k) if centroid_root else None
    if source is not None:
        log(f"Loading centroid cache for k={k}: {source}")
        centroids = load_centroids(source, int(embeddings.shape[1]))
        if centroids.shape[0] != k:
            raise ValueError(f"{source} contains k={centroids.shape[0]}, requested k={k}")
    else:
        log(f"Fitting spherical k-means codebook for k={k}")
        centroids = fit_centroids(embeddings, k, sample_size, niter, seed)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache_path, centroids)
    return centroids


def build_codes(
    embeddings: np.ndarray,
    centroids: np.ndarray,
    output_path: Path,
    batch_size: int,
) -> np.memmap:
    code_dtype = np.uint16 if centroids.shape[0] <= np.iinfo(np.uint16).max else np.uint32
    shape = (embeddings.shape[0],)
    if not cache_matches(output_path, shape, np.dtype(code_dtype)):
        output_path.parent.mkdir(parents=True, exist_ok=True)
        codes = np.lib.format.open_memmap(output_path, mode="w+", dtype=code_dtype, shape=shape)
        codebook = np.ascontiguousarray(centroids, dtype=np.float32)
        for start in range(0, embeddings.shape[0], batch_size):
            stop = min(start + batch_size, embeddings.shape[0])
            values = normalize_rows(np.asarray(embeddings[start:stop]))
            codes[start:stop] = np.argmax(values @ codebook.T, axis=1).astype(code_dtype)
        codes.flush()
        del codes
    return np.load(output_path, mmap_mode="r")


def iter_encoded_windows(
    codes: np.ndarray,
    records: list[Record],
    centroids: np.ndarray,
    window_size: int,
    stride: int,
) -> Iterator[tuple[int, np.ndarray]]:
    codebook = np.ascontiguousarray(centroids, dtype=np.float32)
    for record in records:
        if record.n_windows == 0:
            continue
        node_values = np.asarray(codes[record.node_start : record.node_start + record.length])
        views = np.lib.stride_tricks.sliding_window_view(node_values, window_size)[::stride]
        windows = np.ascontiguousarray(
            codebook[views].reshape(record.n_windows, window_size * codebook.shape[1]),
            dtype=np.float32,
        )
        yield record.window_start, windows


def query_encoded_windows(
    codes: np.ndarray,
    queries: list[Query],
    records: list[Record],
    centroids: np.ndarray,
    window_size: int,
) -> np.ndarray:
    by_id = {record.identifier: record for record in records}
    codebook = np.ascontiguousarray(centroids, dtype=np.float32)
    output = np.empty((len(queries), window_size * codebook.shape[1]), dtype=np.float32)
    for row, query in enumerate(queries):
        record = by_id[query.transcript_id]
        start = record.node_start + query.window_offset
        values = np.asarray(codes[start : start + window_size], dtype=np.int64)
        output[row] = codebook[values].reshape(-1)
    return output


def build_hnsw(
    codes: np.ndarray,
    records: list[Record],
    centroids: np.ndarray,
    window_size: int,
    stride: int,
    m: int,
    ef_construction: int,
    num_threads: int,
    output_path: Path,
):
    import hnswlib

    dimension = window_size * int(centroids.shape[1])
    if output_path.exists():
        index = hnswlib.Index(space="ip", dim=dimension)
        index.load_index(str(output_path), max_elements=int(records[-1].window_start + records[-1].n_windows))
        if num_threads > 0:
            index.set_num_threads(int(num_threads))
        return index, 0.0
    started = time.perf_counter()
    index = hnswlib.Index(space="ip", dim=dimension)
    index.init_index(
        max_elements=int(records[-1].window_start + records[-1].n_windows),
        M=int(m),
        ef_construction=int(ef_construction),
        random_seed=1,
    )
    if num_threads > 0:
        index.set_num_threads(int(num_threads))
    for base, vectors in iter_encoded_windows(codes, records, centroids, window_size, stride):
        labels = np.arange(base, base + vectors.shape[0], dtype=np.int64)
        index.add_items(vectors, labels)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    index.save_index(str(output_path))
    return index, time.perf_counter() - started


def recall_at(exact: np.ndarray, approximate: np.ndarray, k: int) -> float:
    exact = exact[:, :k]
    approximate = approximate[:, :k]
    values = [len(set(row_a.tolist()).intersection(row_b.tolist())) / float(k) for row_a, row_b in zip(exact, approximate)]
    return float(np.mean(values))


def self_hit_rate(labels: np.ndarray, query_ids: np.ndarray, k: int) -> float:
    return float(np.mean(np.any(labels[:, :k] == query_ids[:, None], axis=1)))


def positive_rank(labels: np.ndarray, query_ids: np.ndarray, k: int) -> float:
    ranks = []
    for row, target in zip(labels[:, :k], query_ids):
        matches = np.flatnonzero(row == target)
        ranks.append(float(matches[0] + 1) if len(matches) else float(k + 1))
    return float(np.mean(ranks))


def file_size(path: Path) -> int:
    return path.stat().st_size if path.exists() else 0


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> int:
    ks = parse_int_list(args.k_values)
    ms = parse_int_list(args.m_values)
    efs = parse_int_list(args.ef_search_values)
    records, node_count, window_count = load_structures(args.structures, args.window_size, args.stride)
    queries, query_stats = load_queries(args.queries, records, args.max_queries)
    embeddings = source_embeddings(args.embeddings, node_count, args.original_dtype)
    if embeddings.shape[1] != 128:
        log(f"Note: benchmark embedding dimension is {embeddings.shape[1]}, not 128")
    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata = {
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
    }
    (args.outdir / "benchmark.json").write_text(json.dumps(metadata, indent=2) + "\n")

    flat_windows_path = args.outdir / f"flat_windows_{args.original_dtype}_w{args.window_size}_s{args.stride}.npy"
    windows = build_flat_windows(
        embeddings,
        records,
        args.window_size,
        args.stride,
        args.original_dtype,
        flat_windows_path,
    )
    flat_index, flat_build_seconds = build_flat_index(
        windows, args.outdir / "flatip.index", args.batch_size
    )
    query_ids = np.asarray([query.window_id for query in queries], dtype=np.int64)
    query_vectors = np.ascontiguousarray(windows[query_ids], dtype=np.float32)
    exact_started = time.perf_counter()
    exact_distances, exact_labels = flat_index.search(query_vectors, args.top_k)
    flat_query_seconds = time.perf_counter() - exact_started
    baseline = {
        "backend": "FlatIP",
        "index_type": "IndexFlatIP",
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
    }
    (args.outdir / "flatip.json").write_text(json.dumps(baseline, indent=2) + "\n")

    rows: list[dict[str, object]] = []
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
        codes = build_codes(
            embeddings,
            centroids,
            args.outdir / f"codes_k{k}.npy",
            args.batch_size,
        )
        query_vectors_quantized = query_encoded_windows(
            codes, queries, records, centroids, args.window_size
        )
        for m in ms:
            index_path = args.outdir / f"hnsw_k{k}_m{m}_efc{args.ef_construction}.bin"
            index, hnsw_build_seconds = build_hnsw(
                codes,
                records,
                centroids,
                args.window_size,
                args.stride,
                m,
                args.ef_construction,
                args.num_threads,
                index_path,
            )
            for ef_search in efs:
                index.set_ef(int(ef_search))
                if args.num_threads > 0:
                    index.set_num_threads(int(args.num_threads))
                started = time.perf_counter()
                labels, distances = index.knn_query(query_vectors_quantized, k=args.top_k)
                query_seconds = time.perf_counter() - started
                fields: dict[str, object] = {
                    "backend": "hnswlib",
                    "k": k,
                    "m": m,
                    "ef_construction": args.ef_construction,
                    "ef_search": ef_search,
                    "n_queries": len(queries),
                    "n_windows": window_count,
                    "window_size": args.window_size,
                    "flat_build_seconds": flat_build_seconds,
                    "flat_query_seconds": flat_query_seconds,
                    "hnsw_build_seconds": hnsw_build_seconds,
                    "hnsw_query_seconds": query_seconds,
                    "hnsw_index_bytes": file_size(index_path),
                    "recall_at_1": recall_at(exact_labels, labels, 1),
                    "recall_at_5": recall_at(exact_labels, labels, min(5, args.top_k)),
                    "recall_at_10": recall_at(exact_labels, labels, min(10, args.top_k)),
                    "recall_at_50": recall_at(exact_labels, labels, min(50, args.top_k)),
                    "self_hit_at_1": self_hit_rate(labels, query_ids, 1),
                    "self_hit_at_5": self_hit_rate(labels, query_ids, min(5, args.top_k)),
                    "self_hit_at_10": self_hit_rate(labels, query_ids, min(10, args.top_k)),
                    "self_hit_at_50": self_hit_rate(labels, query_ids, min(50, args.top_k)),
                    "mean_self_rank_at_50": positive_rank(labels, query_ids, min(50, args.top_k)),
                    "mean_exact_top1_score": float(np.mean(exact_distances[:, 0])),
                    "mean_hnsw_top1_score": float(np.mean(1.0 - distances[:, 0])),
                }
                rows.append(fields)
                log(
                    f"k={k} M={m} ef={ef_search}: "
                    f"R@10={fields['recall_at_10']:.4f} "
                    f"R@50={fields['recall_at_50']:.4f} "
                    f"self@1={fields['self_hit_at_1']:.4f} "
                    f"query={query_seconds:.3f}s"
                )
    write_csv(args.outdir / "results.csv", rows)
    (args.outdir / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
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
    parser.add_argument("--max-queries", type=int)
    parser.add_argument("--original-dtype", choices=("float16", "float32"), default="float16")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.window_size < 1 or args.stride < 1 or args.batch_size < 1 or args.top_k < 1:
            raise ValueError("window size, stride, batch size, and top-k must be >= 1")
        if args.ef_construction < 1 or args.num_threads < 0:
            raise ValueError("ef-construction must be >= 1 and num-threads must be >= 0")
        return run(args)
    except (OSError, KeyError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
