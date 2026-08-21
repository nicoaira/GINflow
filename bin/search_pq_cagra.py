#!/usr/bin/env python3
"""Search a PQ-CAGRA index with ADC (GPU or CPU)."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from node_quantization import load_json, normalize_rows


def pair_window_shards(windows_dir: Path) -> list[tuple[Path, Path]]:
    window_paths = sorted(windows_dir.glob("*.windows.npz"))
    manifest_paths = sorted(windows_dir.glob("*.windows.manifest.json"))
    window_map = {path.name[: -len(".windows.npz")]: path for path in window_paths}
    manifest_map = {path.name[: -len(".windows.manifest.json")]: path for path in manifest_paths}
    if set(window_map) != set(manifest_map) or not window_map:
        raise ValueError(f"no query windows in {windows_dir}")
    return [(window_map[key], manifest_map[key]) for key in sorted(window_map)]


def load_targets(path: Path) -> list[tuple[str, int, int]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[tuple[str, int, int] | None] = []
        for record in reader:
            faiss_id = int(record["faiss_id"])
            while len(rows) <= faiss_id:
                rows.append(None)
            rows[faiss_id] = (record["transcript_id"], int(record["start"]), int(record["end"]))
    return [(row[0], row[1], row[2]) for row in rows if row is not None]


def load_float_windows(windows_dir: Path) -> tuple[np.ndarray, list[tuple[str, int, int]]]:
    blocks: list[np.ndarray] = []
    mapping: list[tuple[str, int, int]] = []
    for window_path, manifest_path in pair_window_shards(windows_dir):
        manifest = load_json(manifest_path)
        window_size = int(manifest["window_size"])
        stride = int(manifest["stride"])
        with np.load(window_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                values = np.ascontiguousarray(arrays[identifier], dtype=np.float32)
                blocks.append(values)
                for offset in range(values.shape[0]):
                    start = offset * stride
                    mapping.append((identifier, start, start + window_size))
    return np.concatenate(blocks, axis=0), mapping


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--min-similarity", type=float, default=float("-inf"))
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--itopk-size", type=int, default=0)
    parser.add_argument("--num-threads", type=int, default=0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import pq_cagra_adc as pq
    except ImportError as exc:
        raise SystemExit(
            "pq_cagra_adc is not installed. Use -profile conda or -profile docker "
            "so SEARCH_PQ_CAGRA gets nicolas.aira::pq-cagra-adc."
        ) from exc

    db_meta = load_json(args.database / "meta.json")
    queries, query_map = load_float_windows(args.windows_dir)
    window_size = int(db_meta["window_size"])
    dim = int(db_meta.get("embedding_dim", 128))
    if queries.ndim != 2 or queries.shape[1] != window_size * dim:
        raise SystemExit(f"query width {queries.shape} does not match {window_size}x{dim}")
    queries = queries.reshape(queries.shape[0], window_size, dim)
    queries = np.stack([normalize_rows(window) for window in queries]).astype(np.float32, copy=False)
    rotation = None
    rotation_path = args.database / "quantization" / "rotation.npy"
    if rotation_path.exists():
        rotation = np.load(rotation_path)
    index = pq.load_index(args.database / "cagra.index")
    itopk = args.itopk_size or int(db_meta.get("itopk_size") or max(args.k, 64))
    kwargs: dict = {}
    if args.device == "cpu":
        if args.num_threads:
            kwargs["num_threads"] = args.num_threads
    else:
        kwargs["itopk_size"] = max(itopk, args.k)
    labels, distances = pq.search(
        index,
        queries,
        k=args.k,
        device=args.device,
        rotation=rotation,
        **kwargs,
    )
    targets = load_targets(args.database / "windows.tsv")
    hits = []
    for query_id, (query_name, query_start, query_end) in enumerate(query_map):
        for rank in range(args.k):
            label = int(labels[query_id, rank])
            if label < 0 or label >= len(targets):
                continue
            score = float(-distances[query_id, rank]) / float(window_size)
            if score < args.min_similarity:
                continue
            target_id, target_start, target_end = targets[label]
            hits.append(
                (
                    query_name,
                    query_start,
                    query_end,
                    target_id,
                    target_start,
                    target_end,
                    score,
                    rank + 1,
                )
            )
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["query_id", "query_start", "query_end", "target_id", "target_start", "target_end", "score", "rank"]
        )
        writer.writerows(hits)
    print(json.dumps({"n_seeds": len(hits), "device": args.device}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
