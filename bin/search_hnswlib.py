#!/usr/bin/env python3
"""ADC search of a node-PQ HNSWLIB index using original float query windows."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
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
    only_windows = sorted(set(window_map) - set(manifest_map))
    only_manifests = sorted(set(manifest_map) - set(window_map))
    if only_windows or only_manifests:
        raise ValueError(
            "window / manifest shard names do not match: "
            f"only-windows={only_windows} only-manifests={only_manifests}"
        )
    if not window_map:
        raise ValueError(f"no query windows found in {windows_dir}")
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


def load_float_windows(windows_dir: Path) -> tuple[np.ndarray, list[tuple[str, int, int]], dict]:
    mapping: list[tuple[str, int, int]] = []
    blocks: list[np.ndarray] = []
    reference = None
    for window_path, manifest_path in pair_window_shards(windows_dir):
        manifest = load_json(manifest_path)
        if reference is None:
            reference = manifest
        window_size = int(manifest["window_size"])
        stride = int(manifest["stride"])
        with np.load(window_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                values = np.asarray(arrays[identifier])
                if values.ndim != 2:
                    raise ValueError(f"{identifier} has invalid window shape {values.shape}")
                blocks.append(np.ascontiguousarray(values, dtype=np.float32))
                for offset in range(values.shape[0]):
                    start = offset * stride
                    mapping.append((identifier, start, start + window_size))
    if not blocks or reference is None:
        raise ValueError("no query windows found")
    stacked = np.concatenate(blocks, axis=0)
    return stacked, mapping, reference


def search_windows(
    query_windows_dir: Path,
    database: Path,
    executable: Path,
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
    queries, query_map, query_meta = load_float_windows(query_windows_dir)
    window_size = int(db_meta["window_size"])
    dim = int(db_meta.get("embedding_dim", 128))
    if queries.shape[1] != window_size * dim:
        raise ValueError(
            f"query windows have width {queries.shape[1]}, expected {window_size * dim}"
        )
    queries = queries.reshape(queries.shape[0], window_size, dim)
    queries = np.concatenate(
        [normalize_rows(queries[i])[None, ...] for i in range(queries.shape[0])], axis=0
    ).astype(np.float32, copy=False)
    import tempfile

    tmp = Path(tempfile.mkdtemp(prefix="ginflow-hnsw-search-"))
    query_bin = tmp / "query_windows.bin"
    labels_bin = tmp / "labels.bin"
    distances_bin = tmp / "distances.bin"
    queries.reshape(queries.shape[0], -1).tofile(query_bin)
    codebook = np.load(database / "quantization" / "codebook.npy")
    codebook_bin = database / "quantization" / "codebook.bin"
    if not codebook_bin.exists():
        np.ascontiguousarray(codebook, dtype=np.float32).tofile(codebook_bin)
    rotation = database / "quantization" / "rotation.npy"
    rotation_bin = database / "quantization" / "rotation.bin"
    if rotation.exists() and not rotation_bin.exists():
        np.ascontiguousarray(np.load(rotation), dtype=np.float32).tofile(rotation_bin)
    search_ef = int(ef_search if ef_search is not None else db_meta.get("hnsw_ef_search", 200))
    command = [
        str(executable),
        "search",
        "--queries",
        str(query_bin),
        "--codebook",
        str(codebook_bin),
        "--index",
        str(database / "index.bin"),
        "--query-count",
        str(queries.shape[0]),
        "--window-size",
        str(window_size),
        "--pq-m",
        str(int(db_meta["pq_m"])),
        "--nbits",
        str(int(db_meta["pq_nbits"])),
        "--dim",
        str(dim),
        "--k",
        str(k),
        "--ef-search",
        str(search_ef),
        "--threads",
        str(num_threads),
        "--labels-out",
        str(labels_bin),
        "--distances-out",
        str(distances_bin),
    ]
    if rotation.exists():
        command.extend(["--rotation", str(rotation_bin)])
    try:
        subprocess.run(command, check=True)
        labels = np.fromfile(labels_bin, dtype=np.uint64).reshape(queries.shape[0], k)
        distances = np.fromfile(distances_bin, dtype=np.float32).reshape(queries.shape[0], k)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    targets = load_targets(database / "windows.tsv")
    hits: list[tuple] = []
    for query_id, (query_name, query_start, query_end) in enumerate(query_map):
        for rank in range(k):
            label = int(labels[query_id, rank])
            if label < 0 or label >= len(targets) or label == np.iinfo(np.uint64).max:
                continue
            # ADC distances are negative inner products.
            score = float(-distances[query_id, rank]) / float(window_size)
            if score < min_similarity:
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
    _ = query_meta
    return hits


def write_seeds(path: Path, hits: list[tuple]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["query_id", "query_start", "query_end", "target_id", "target_start", "target_end", "score", "rank"]
        )
        writer.writerows(hits)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--min-similarity", type=float, default=float("-inf"))
    parser.add_argument("--ef-search", type=int, default=None)
    parser.add_argument("--num-threads", type=int, default=0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    hits = search_windows(
        args.windows_dir,
        args.database,
        args.executable,
        args.k,
        args.min_similarity,
        args.ef_search,
        args.num_threads,
    )
    write_seeds(args.output, hits)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
