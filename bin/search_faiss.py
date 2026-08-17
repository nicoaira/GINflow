#!/usr/bin/env python3
"""Search a FAISS window index and emit seeds above a cosine threshold."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from faiss_index import (
    IndexOptions,
    distances_to_similarity,
    prepare_search_index,
)
from cuvs_index import is_cuvs_database, load_index as load_cuvs_index
from ngt_index import is_ngt_database, load_index as load_ngt_index
from scann_index import (
    apply_search_params as apply_scann_search_params,
    is_scann_database,
    load_index as load_scann_index,
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
        raise ValueError(
            "query windows are incompatible with the vector database: " + "; ".join(mismatches)
        )


def search_shard(
    windows_path: Path,
    manifest: dict,
    index: Any,
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
    parser.add_argument(
        "--database",
        type=Path,
        required=True,
        help="Directory with windows.tsv, meta.json, and either index.faiss or scann/",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--min-similarity", type=float, default=0.8)
    parser.add_argument("--nprobe", type=int)
    parser.add_argument("--cuvs-n-probes", type=int)
    parser.add_argument("--hnsw-ef-search", type=int)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-device", type=int, default=0)
    parser.add_argument("--scann-reorder", type=int)
    parser.add_argument("--scann-leaves-to-search", type=int)
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
        targets = load_targets(args.database / "windows.tsv")
        search_options = IndexOptions(
            nprobe=args.nprobe,
            hnsw_ef_search=args.hnsw_ef_search,
            gpu=args.gpu,
            gpu_device=args.gpu_device,
            scann_reorder=args.scann_reorder if args.scann_reorder is not None else 100,
        )
        if is_cuvs_database(args.database, db_meta):
            if not args.gpu:
                raise ValueError("cuVS databases require GPU search; use -profile gpu")
            index = load_cuvs_index(args.database / "cuvs", db_meta, args.cuvs_n_probes)
            metric, lsh_nbits = str(db_meta.get("metric") or index.metric), None
        elif is_ngt_database(args.database, db_meta):
            index = load_ngt_index(args.database / "ngt", db_meta)
            metric, lsh_nbits = str(db_meta.get("metric") or index.metric), None
        elif is_scann_database(args.database, db_meta):
            if args.gpu:
                raise ValueError(
                    "--faiss_gpu is not supported for ScaNN. ScaNN is CPU-only (AVX/FMA)."
                )
            ntotal = db_meta.get("n_windows")
            index = load_scann_index(
                args.database / "scann",
                ntotal=int(ntotal) if ntotal is not None else None,
                leaves_to_search=db_meta.get("nprobe"),
                reorder=db_meta.get("scann_reorder"),
            )
            apply_scann_search_params(
                index,
                nprobe=args.scann_leaves_to_search,
                reorder=args.scann_reorder,
            )
            metric, lsh_nbits = "inner_product", None
        else:
            import faiss

            index = faiss.read_index(str(args.database / "index.faiss"))
            index, metric, lsh_nbits = prepare_search_index(index, db_meta, search_options)
        if index.ntotal != len(targets):
            raise ValueError(
                f"index ntotal={index.ntotal} does not match windows.tsv rows={len(targets)}"
            )
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
