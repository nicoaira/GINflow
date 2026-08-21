#!/usr/bin/env python3
"""Pack published GINflow artifacts into resumable ANN benchmark inputs."""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_hnswlib import (  # noqa: E402
    build_flat_windows,
    load_queries,
    load_structures,
    source_embeddings,
)
from rerank_candidates import load_query_windows  # noqa: E402


def load_records(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "transcript_id" not in reader.fieldnames:
            raise ValueError(f"{path} must contain transcript_id")
        return [str(row["transcript_id"]) for row in reader]


def pack_nodes(
    database: Path,
    identifiers: list[str],
    output: Path,
    dtype: np.dtype,
) -> np.memmap:
    with np.load(database / "embeddings.npz") as archive:
        arrays = {identifier: np.asarray(archive[identifier]) for identifier in archive.files}
    missing = [identifier for identifier in identifiers if identifier not in arrays]
    if missing:
        raise KeyError(f"database embeddings are missing {missing[:5]}")
    total_nodes = sum(int(arrays[identifier].shape[0]) for identifier in identifiers)
    dimension = int(arrays[identifiers[0]].shape[1])
    output.parent.mkdir(parents=True, exist_ok=True)
    nodes = np.lib.format.open_memmap(output, mode="w+", dtype=dtype, shape=(total_nodes, dimension))
    start = 0
    for identifier in identifiers:
        values = np.asarray(arrays[identifier], dtype=dtype)
        stop = start + values.shape[0]
        nodes[start:stop] = values
        start = stop
    nodes.flush()
    return np.load(output, mmap_mode="r")


def pack_query_vectors(
    query_windows_dir: Path,
    query_manifests_dir: Path,
    query_structures: Path,
    query_selections: Path,
    window_size: int,
    stride: int,
    output: Path,
) -> tuple[np.ndarray, dict[str, int]]:
    records, _nodes, _windows = load_structures(query_structures, window_size, stride)
    queries, stats = load_queries(query_selections, records, None)
    paths = sorted(query_windows_dir.rglob("*.windows.npz"))
    manifests = sorted(query_manifests_dir.rglob("*.windows.manifest.json"))
    vectors = load_query_windows(paths, manifests)
    selected = []
    for query in queries:
        key = (
            query.transcript_id,
            query.window_offset * stride,
            query.window_offset * stride + window_size,
        )
        if key not in vectors:
            raise KeyError(f"query window {key} is missing from published windows")
        selected.append(vectors[key])
    result = np.ascontiguousarray(np.stack(selected), dtype=np.float32)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.save(output, result)
    return result, stats


def build_reference(
    database_windows: np.ndarray,
    queries: np.ndarray,
    output: Path,
    k: int,
    batch_size: int,
) -> dict:
    if k < 1:
        raise ValueError("reference k must be positive")
    if batch_size < 1:
        raise ValueError("reference batch size must be positive")
    if database_windows.ndim != 2 or queries.ndim != 2:
        raise ValueError("reference arrays must be 2-D")
    if database_windows.shape[1] != queries.shape[1]:
        raise ValueError("database and query dimensions differ")
    k = min(k, database_windows.shape[0])
    started = time.perf_counter()
    best_scores = np.full((queries.shape[0], k), -np.inf, dtype=np.float32)
    best_labels = np.full((queries.shape[0], k), -1, dtype=np.int64)
    query_values = np.ascontiguousarray(queries, dtype=np.float32)
    for start in range(0, database_windows.shape[0], batch_size):
        stop = min(start + batch_size, database_windows.shape[0])
        block = np.ascontiguousarray(database_windows[start:stop], dtype=np.float32)
        scores = query_values @ block.T
        local_k = min(k, scores.shape[1])
        local_indices = np.argpartition(scores, scores.shape[1] - local_k, axis=1)[:, -local_k:]
        local_scores = np.take_along_axis(scores, local_indices, axis=1)
        local_order = np.argsort(-local_scores, axis=1, kind="stable")
        local_scores = np.take_along_axis(local_scores, local_order, axis=1)
        local_labels = start + np.take_along_axis(local_indices, local_order, axis=1)
        merged_scores = np.concatenate((best_scores, local_scores), axis=1)
        merged_labels = np.concatenate((best_labels, local_labels), axis=1)
        chosen = np.argpartition(merged_scores, merged_scores.shape[1] - k, axis=1)[:, -k:]
        best_scores = np.take_along_axis(merged_scores, chosen, axis=1)
        best_labels = np.take_along_axis(merged_labels, chosen, axis=1)
        order = np.argsort(-best_scores, axis=1, kind="stable")
        best_scores = np.take_along_axis(best_scores, order, axis=1)
        best_labels = np.take_along_axis(best_labels, order, axis=1)
    np.save(output, best_labels)
    elapsed = time.perf_counter() - started
    return {
        "backend": "FlatIP-equivalent",
        "index_type": "blockwise-IndexFlatIP",
        "k": k,
        "build_seconds": 0.0,
        "search_seconds": elapsed,
        "batch_size": batch_size,
        "index_path": None,
    }


def run(args: argparse.Namespace) -> dict:
    db_records = load_records(args.database / "records.tsv")
    source_records, node_count, window_count = load_structures(
        args.structures, args.window_size, args.stride
    )
    source_ids = [record.identifier for record in source_records]
    db_id_set = set(db_records)
    missing = [identifier for identifier in source_ids if identifier not in db_id_set]
    extra = [identifier for identifier in db_records if identifier not in set(source_ids)]
    if missing or extra or len(source_records) != len(db_records):
        raise ValueError(
            "--structures and database/records.tsv contain different identifiers "
            f"(missing={missing[:3]}, extra={extra[:3]})"
        )
    nodes = pack_nodes(
        args.database,
        source_ids,
        args.outdir / "node_embeddings.float16.npy",
        np.dtype(args.original_dtype),
    )
    records = source_records
    database_windows_path = (
        args.outdir
        / f"database_windows_{args.original_dtype}_w{args.window_size}_s{args.stride}_source_order.npy"
    )
    database_windows = build_flat_windows(
        nodes,
        records,
        args.window_size,
        args.stride,
        args.original_dtype,
        database_windows_path,
    )
    queries, query_stats = pack_query_vectors(
        args.query_windows_dir,
        args.query_manifests_dir,
        args.query_structures,
        args.query_selections,
        args.window_size,
        args.stride,
        args.outdir / "query_windows.float32.npy",
    )
    reference = None
    if not args.skip_reference:
        reference = build_reference(
            database_windows,
            queries,
            args.outdir / "reference_labels.npy",
            args.reference_k,
            args.reference_batch_size,
        )
    metadata = {
        "database": str(args.database),
        "structures": str(args.structures),
        "query_structures": str(args.query_structures),
        "query_selections": str(args.query_selections),
        "query_windows_dir": str(args.query_windows_dir),
        "query_manifests_dir": str(args.query_manifests_dir),
        "original_dtype": args.original_dtype,
        "node_embeddings": str(args.outdir / "node_embeddings.float16.npy"),
        "database_windows": str(database_windows_path),
        "window_size": args.window_size,
        "stride": args.stride,
        "node_count": node_count,
        "window_count": window_count,
        "n_queries": int(queries.shape[0]),
        "embedding_dim": int(nodes.shape[1]),
        "query_stats": query_stats,
        "reference": reference,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    return metadata


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--query-structures", type=Path, required=True)
    parser.add_argument("--query-selections", type=Path, required=True)
    parser.add_argument("--query-windows-dir", type=Path, required=True)
    parser.add_argument("--query-manifests-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--original-dtype", choices=("float16", "float32"), default="float16")
    parser.add_argument("--reference-k", type=int, default=500)
    parser.add_argument("--reference-batch-size", type=int, default=32768)
    parser.add_argument("--skip-reference", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.window_size < 1 or args.stride < 1:
            raise ValueError("window-size and stride must be positive")
        run(args)
    except (OSError, KeyError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
