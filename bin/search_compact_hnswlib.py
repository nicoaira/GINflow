#!/usr/bin/env python3
"""Search a compact custom-distance hnswlib index over quantized windows."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compact_hnswlib_io import pack_code_windows, read_compact_results, write_similarity_binary
from hnswlib_index import load_quantization
from search_hnswlib import check_compatible, load_targets, write_seeds


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def load_embeddings(paths: list[Path]) -> dict[str, np.ndarray]:
    embeddings: dict[str, np.ndarray] = {}
    for path in paths:
        with np.load(path) as archive:
            for identifier in archive.files:
                if identifier in embeddings:
                    raise ValueError(f"duplicate embedding identifier {identifier!r}")
                embeddings[identifier] = np.asarray(archive[identifier])
    return embeddings


def normalized_window(
    embeddings: dict[str, np.ndarray], identifier: str, start: int, end: int
) -> np.ndarray:
    if identifier not in embeddings:
        raise KeyError(f"missing original embedding for {identifier}")
    values = np.asarray(embeddings[identifier][start:end], dtype=np.float32)
    if values.shape[0] != end - start:
        raise ValueError(f"window [{start}, {end}) is outside the embedding for {identifier}")
    flat = np.ascontiguousarray(values.reshape(-1), dtype=np.float32)
    norm = np.linalg.norm(flat)
    return flat / max(float(norm), 1e-12)


def rerank_candidates(
    labels: np.ndarray,
    query_mapping: list[tuple[str, int, int]],
    targets: list[tuple[str, int, int]],
    query_embeddings: dict[str, np.ndarray],
    database_embeddings: dict[str, np.ndarray],
    output_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    if labels.ndim != 2 or labels.shape[0] != len(query_mapping):
        raise ValueError("candidate labels must be a 2-D array with one row per query window")
    if not query_mapping:
        return (
            np.empty((0, output_k), dtype=np.uint64),
            np.empty((0, output_k), dtype=np.float32),
        )
    window_size = int(query_mapping[0][2] - query_mapping[0][1])
    if window_size < 1 or any(int(end - start) != window_size for _id, start, end in query_mapping):
        raise ValueError("all query windows must have the same positive length")
    if any(int(end - start) != window_size for _id, start, end in targets):
        raise ValueError("all target windows must have the same length as query windows")

    query_vectors = np.ascontiguousarray(
        np.stack(
            [normalized_window(query_embeddings, identifier, start, end) for identifier, start, end in query_mapping]
        ),
        dtype=np.float32,
    )
    query_count, candidate_count = labels.shape
    if output_k < 1 or candidate_count < output_k:
        raise ValueError("candidate count must be >= the requested output k")
    flat_labels = np.asarray(labels, dtype=np.uint64).reshape(-1)
    invalid = np.iinfo(np.uint64).max
    valid = flat_labels != invalid
    if np.any(flat_labels[valid] >= len(targets)):
        raise ValueError("candidate labels contain an invalid database element")

    target_identifier_to_group: dict[str, int] = {}
    target_group_identifiers: list[str] = []
    target_group = np.empty(len(targets), dtype=np.int32)
    target_starts = np.empty(len(targets), dtype=np.int64)
    for target_index, (target_id, target_start, _target_end) in enumerate(targets):
        group = target_identifier_to_group.get(target_id)
        if group is None:
            group = len(target_group_identifiers)
            target_identifier_to_group[target_id] = group
            target_group_identifiers.append(target_id)
        target_group[target_index] = group
        target_starts[target_index] = target_start

    flat_scores = np.full(flat_labels.shape[0], -np.inf, dtype=np.float32)
    valid_positions = np.flatnonzero(valid)
    if valid_positions.size:
        candidate_groups = target_group[flat_labels[valid_positions].astype(np.int64)]
        sorted_valid_order = np.argsort(candidate_groups, kind="stable")
        sorted_positions = valid_positions[sorted_valid_order]
        sorted_groups = candidate_groups[sorted_valid_order]
        group_starts = np.r_[0, np.flatnonzero(np.diff(sorted_groups)) + 1, sorted_positions.size]
        offsets = np.arange(window_size, dtype=np.int64)
        flat_query_rows = np.arange(query_count, dtype=np.int64)[:, None]
        flat_query_rows = np.broadcast_to(flat_query_rows, labels.shape).reshape(-1)
        for start, stop in zip(group_starts[:-1], group_starts[1:]):
            positions = sorted_positions[start:stop]
            group = int(sorted_groups[start])
            target_id = target_group_identifiers[group]
            embedding = np.ascontiguousarray(database_embeddings[target_id], dtype=np.float32)
            target_labels = flat_labels[positions].astype(np.int64)
            window_indices = target_starts[target_labels, None] + offsets[None, :]
            values = np.ascontiguousarray(embedding[window_indices], dtype=np.float32).reshape(len(positions), -1)
            values /= np.maximum(np.linalg.norm(values, axis=1, keepdims=True), np.float32(1e-12))
            flat_scores[positions] = np.einsum(
                "nd,nd->n", values, query_vectors[flat_query_rows[positions]], optimize=True
            )

    candidate_scores = flat_scores.reshape(query_count, candidate_count)
    output_labels = np.empty((query_count, output_k), dtype=np.uint64)
    output_scores = np.empty((query_count, output_k), dtype=np.float32)
    order = np.argpartition(-candidate_scores, output_k - 1, axis=1)[:, :output_k]
    order_scores = np.take_along_axis(candidate_scores, order, axis=1)
    order = np.take_along_axis(order, np.argsort(-order_scores, axis=1, kind="stable"), axis=1)
    output_labels[:] = np.take_along_axis(labels, order, axis=1)
    output_scores[:] = np.take_along_axis(candidate_scores, order, axis=1)
    return output_labels, output_scores


def search_windows(
    query_windows_dir: Path,
    database: Path,
    executable: Path,
    k: int,
    candidate_k: int | None,
    rerank: bool,
    query_embeddings_paths: list[Path] | None,
    min_similarity: float,
    ef_search: int | None,
    num_threads: int,
) -> list[tuple]:
    if k < 1:
        raise ValueError("--k must be >= 1")
    candidate_k = int(candidate_k or k)
    if candidate_k < k:
        raise ValueError("--candidate-k must be >= --k")
    if rerank and not query_embeddings_paths:
        raise ValueError("--rerank requires --query-embeddings")
    db_meta = load_json(database / "meta.json")
    if str(db_meta.get("index_type", "")).upper() != "HNSWLIB_COMPACT":
        raise ValueError(f"{database} is not a compact HNSWLIB database")
    centroids, similarity, _metadata = load_quantization(database / "quantization")
    if similarity.shape != (centroids.shape[0], centroids.shape[0]):
        raise ValueError("centroid similarity matrix has an invalid shape")
    targets = load_targets(database / "windows.tsv")
    search_k = min(candidate_k, len(targets))
    query_count = 0
    executable_path = executable if executable.is_absolute() else Path.cwd() / executable
    with tempfile.TemporaryDirectory(prefix="compact-hnsw-query-") as work:
        workdir = Path(work)
        query_codes = workdir / "query_codes.bin"
        query_count, window_size, query_mapping, _manifests = pack_code_windows(query_windows_dir, query_codes)
        for manifest in _manifests:
            check_compatible(manifest, db_meta)
        similarity_path = workdir / "similarity.bin"
        write_similarity_binary(database / "quantization", similarity_path)
        labels_path = workdir / "labels.bin"
        distances_path = workdir / "distances.bin"
        command = [
            str(executable_path),
            "search",
            "--codes", str(query_codes),
            "--similarity", str(similarity_path),
            "--query-count", str(query_count),
            "--window-size", str(window_size),
            "--n-centroids", str(centroids.shape[0]),
            "--index", str(database / "index.bin"),
            "--k", str(search_k),
            "--ef-search", str(ef_search if ef_search is not None else db_meta.get("hnsw_ef_search", 100)),
            "--threads", str(num_threads),
            "--labels-out", str(labels_path),
            "--distances-out", str(distances_path),
        ]
        subprocess.run(command, check=True)
        labels, distances = read_compact_results(labels_path, distances_path, query_count, search_k)

    if rerank:
        database_embeddings = load_embeddings([database / "embeddings.npz"])
        query_embeddings = load_embeddings(query_embeddings_paths or [])
        final_labels, final_scores = rerank_candidates(
            labels,
            query_mapping,
            targets,
            query_embeddings,
            database_embeddings,
            k,
        )
    else:
        final_labels = labels[:, :k]
        final_scores = -distances[:, :k]

    hits: list[tuple] = []
    for query_row, ((identifier, query_start, query_end), row_labels, row_scores) in enumerate(
        zip(query_mapping, final_labels, final_scores)
    ):
        rank = 0
        for target_label, score in zip(row_labels, row_scores):
            if int(target_label) == np.iinfo(np.uint64).max:
                continue
            score = float(score)
            if score < min_similarity:
                continue
            rank += 1
            target_id, target_start, target_end = targets[int(target_label)]
            hits.append(
                (
                    identifier,
                    query_start,
                    query_end,
                    target_id,
                    target_start,
                    target_end,
                    score,
                    rank,
                )
            )
            if rank >= k:
                break
    return hits


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--candidate-k", type=int)
    parser.add_argument("--rerank", action="store_true")
    parser.add_argument("--query-embeddings", type=Path, nargs="*")
    parser.add_argument("--min-similarity", type=float, default=0.8)
    parser.add_argument("--ef-search", type=int)
    parser.add_argument("--num-threads", type=int, default=0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        hits = search_windows(
            args.windows_dir,
            args.database,
            args.executable,
            args.k,
            args.candidate_k,
            args.rerank,
            args.query_embeddings,
            args.min_similarity,
            args.ef_search,
            args.num_threads,
        )
        write_seeds(args.output, hits)
    except (OSError, KeyError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"output": str(args.output), "n_seeds": len(hits)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
