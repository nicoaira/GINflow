#!/usr/bin/env python3
"""Search a FAISS window index and emit seeds above a cosine threshold."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import faiss
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from faiss_index import (
    IndexOptions,
    distances_to_similarity,
    prepare_search_index,
)


COMPAT_KEYS = (
    ("window_size", "window_size"),
    ("stride", "window_stride"),
    ("window_dim", "window_dim"),
    ("checkpoint_sha256", "checkpoint_sha256"),
)


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


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


def check_compatible(query_manifest: dict, db_meta: dict) -> None:
    mismatches = []
    for query_key, db_key in COMPAT_KEYS:
        if query_manifest.get(query_key) != db_meta.get(db_key):
            mismatches.append(
                f"{query_key}={query_manifest.get(query_key)!r} vs {db_key}={db_meta.get(db_key)!r}"
            )
    if mismatches:
        raise ValueError("query windows are incompatible with the FAISS database: " + "; ".join(mismatches))


def search_shard(
    windows_path: Path,
    manifest: dict,
    index: faiss.Index,
    targets: list[tuple[str, int, int]],
    k: int,
    min_similarity: float,
    metric: str = "inner_product",
    lsh_nbits: int | None = None,
) -> list[tuple]:
    arrays = np.load(windows_path)
    window_size = int(manifest["window_size"])
    stride = int(manifest["stride"])
    hits: list[tuple] = []
    search_k = min(k, index.ntotal)
    if search_k < 1:
        return hits

    for record in manifest.get("records", []):
        identifier = record["identifier"]
        if identifier not in arrays.files:
            raise KeyError(f"{identifier} is in the query manifest but missing from {windows_path}")
        query = np.ascontiguousarray(arrays[identifier], dtype=np.float32)
        if query.shape[0] == 0:
            continue
        distances, labels = index.search(query, search_k)
        similarities = distances_to_similarity(distances, metric, lsh_nbits)
        for offset, (scores, ids) in enumerate(zip(similarities, labels)):
            query_start = offset * stride
            query_end = query_start + window_size
            rank = 0
            for score, faiss_id in zip(scores, ids):
                if faiss_id < 0:
                    continue
                similarity = float(score)
                if similarity < min_similarity:
                    continue
                rank += 1
                target_id, target_start, target_end = targets[int(faiss_id)]
                hits.append((
                    identifier,
                    query_start,
                    query_end,
                    target_id,
                    target_start,
                    target_end,
                    similarity,
                    rank,
                ))
                if rank >= k:
                    break
    return hits


def write_seeds(path: Path, hits: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "query_id",
            "query_start",
            "query_end",
            "target_id",
            "target_start",
            "target_end",
            "score",
            "rank",
        ])
        for row in hits:
            writer.writerow([
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
                f"{row[6]:.6f}",
                row[7],
            ])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True, help="Directory with index.faiss, windows.tsv, meta.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--min-similarity", type=float, default=0.8)
    parser.add_argument("--nprobe", type=int)
    parser.add_argument("--hnsw-ef-search", type=int)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-device", type=int, default=0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.k < 1:
        print("error: --k must be >= 1", file=sys.stderr)
        return 2
    try:
        query_manifest = load_json(args.manifest)
        db_meta = load_json(args.database / "meta.json")
        check_compatible(query_manifest, db_meta)
        index = faiss.read_index(str(args.database / "index.faiss"))
        targets = load_targets(args.database / "windows.tsv")
        if index.ntotal != len(targets):
            raise ValueError(
                f"index ntotal={index.ntotal} does not match windows.tsv rows={len(targets)}"
            )
        search_options = IndexOptions(
            nprobe=args.nprobe,
            hnsw_ef_search=args.hnsw_ef_search,
            gpu=args.gpu,
            gpu_device=args.gpu_device,
        )
        index, metric, lsh_nbits = prepare_search_index(index, db_meta, search_options)
        hits = search_shard(
            args.windows,
            query_manifest,
            index,
            targets,
            args.k,
            args.min_similarity,
            metric,
            lsh_nbits,
        )
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    write_seeds(args.output, hits)
    print(json.dumps({"output": str(args.output), "n_seeds": len(hits)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
