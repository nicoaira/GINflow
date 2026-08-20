#!/usr/bin/env python3
"""Search a quantized-node HNSWLIB database and emit candidate seeds."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hnswlib_index import encode_code_windows, load_index, load_quantization, similarity_from_distance


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def pair_window_shards(windows_dir: Path) -> list[tuple[Path, Path]]:
    window_paths = sorted(windows_dir.glob("*.windows.npz"))
    manifest_paths = sorted(windows_dir.glob("*.windows.manifest.json"))
    window_map = {path.name[: -len(".windows.npz")]: path for path in window_paths}
    manifest_map = {path.name[: -len(".windows.manifest.json")]: path for path in manifest_paths}
    only_windows = sorted(set(window_map) - set(manifest_map))
    only_manifests = sorted(set(manifest_map) - set(window_map))
    if only_windows or only_manifests:
        raise ValueError(
            "window / manifest shard names do not match: "
            f"only-windows={only_windows} only-manifests={only_manifests}"
        )
    if not window_map:
        raise ValueError(f"no quantized query windows found in {windows_dir}")
    return [(window_map[key], manifest_map[key]) for key in sorted(window_map)]


def load_targets(path: Path) -> list[tuple[str, int, int]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"faiss_id", "transcript_id", "start", "end"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must have columns {sorted(required)}")
        rows: list[tuple[str, int, int] | None] = []
        for record in reader:
            faiss_id = int(record["faiss_id"])
            while len(rows) <= faiss_id:
                rows.append(None)
            rows[faiss_id] = (
                record["transcript_id"],
                int(record["start"]),
                int(record["end"]),
            )
    missing = [idx for idx, row in enumerate(rows) if row is None]
    if missing:
        raise ValueError(f"{path} has gaps at faiss_id {missing[:8]}")
    return [(row[0], row[1], row[2]) for row in rows]  # type: ignore[index]


def check_compatible(manifest: dict, db_meta: dict) -> None:
    keys = (
        ("window_size", "window_size"),
        ("stride", "window_stride"),
        ("window_dim", "window_dim"),
        ("checkpoint_sha256", "checkpoint_sha256"),
        ("n_centroids", "n_centroids"),
        ("embedding_dim", "embedding_dim"),
    )
    mismatches = [
        f"{query_key}={manifest.get(query_key)!r} vs {db_key}={db_meta.get(db_key)!r}"
        for query_key, db_key in keys
        if manifest.get(query_key) != db_meta.get(db_key)
    ]
    if mismatches:
        raise ValueError("query windows are incompatible with the HNSWLIB database: " + "; ".join(mismatches))


def search_windows(
    query_windows_dir: Path,
    database: Path,
    k: int,
    min_similarity: float,
    ef_search: int | None,
    num_threads: int,
) -> list[tuple]:
    if k < 1:
        raise ValueError("--k must be >= 1")
    db_meta = load_json(database / "meta.json")
    if str(db_meta.get("backend", "")).lower() != "hnswlib":
        raise ValueError(f"{database} is not an HNSWLIB database")
    centroids, similarity_matrix, _quantization_metadata = load_quantization(database / "quantization")
    if similarity_matrix.shape != (centroids.shape[0], centroids.shape[0]):
        raise ValueError("centroid similarity matrix has an invalid shape")
    search_ef = int(ef_search if ef_search is not None else db_meta.get("hnsw_ef_search", 100))
    index = load_index(database / "index.bin", int(db_meta["window_dim"]), search_ef, num_threads)
    targets = load_targets(database / "windows.tsv")
    if index.get_current_count() != len(targets):
        raise ValueError(
            f"HNSWLIB element_count={index.get_current_count()} does not match windows.tsv rows={len(targets)}"
        )
    search_k = min(k, len(targets))
    hits: list[tuple] = []
    for window_path, manifest_path in pair_window_shards(query_windows_dir):
        manifest = load_json(manifest_path)
        check_compatible(manifest, db_meta)
        window_size = int(manifest["window_size"])
        stride = int(manifest["stride"])
        with np.load(window_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                if identifier not in arrays.files:
                    raise KeyError(f"{identifier} is in {manifest_path} but missing from {window_path}")
                codes = np.asarray(arrays[identifier])
                vectors = encode_code_windows(codes, centroids)
                if vectors.shape[0] == 0:
                    continue
                labels, distances = index.knn_query(vectors, k=search_k, num_threads=num_threads or -1)
                scores = similarity_from_distance(distances)
                for offset, (row_labels, row_scores) in enumerate(zip(labels, scores)):
                    query_start = offset * stride
                    query_end = query_start + window_size
                    rank = 0
                    for target_label, score in zip(row_labels, row_scores):
                        if int(target_label) < 0:
                            continue
                        value = float(score)
                        if value < min_similarity:
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
                                value,
                                rank,
                            )
                        )
                        if rank >= k:
                            break
    return hits


def write_seeds(path: Path, hits: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "query_id",
                "query_start",
                "query_end",
                "target_id",
                "target_start",
                "target_end",
                "score",
                "rank",
            ]
        )
        for row in hits:
            writer.writerow([row[0], row[1], row[2], row[3], row[4], row[5], f"{row[6]:.6f}", row[7]])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=50)
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
            args.k,
            args.min_similarity,
            args.ef_search,
            args.num_threads,
        )
        write_seeds(args.output, hits)
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"output": str(args.output), "n_seeds": len(hits)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
