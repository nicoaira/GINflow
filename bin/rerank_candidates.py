#!/usr/bin/env python3
"""Exact-rerank ANN seed candidates with preserved original embeddings.

The input candidate table contains database window coordinates, not quantized
vectors.  The database residue embeddings remain the source of truth; node-PQ
codes are used only upstream by candidate indexes.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rerank_core import rerank_label_matrix  # noqa: E402


SEED_COLUMNS = (
    "query_id",
    "query_start",
    "query_end",
    "target_id",
    "target_start",
    "target_end",
    "score",
    "rank",
)


@dataclass
class CandidateGroup:
    query_id: str
    query_start: int
    query_end: int
    targets: list[tuple[str, int, int]]


class OriginalWindowStore:
    """Lazy exact-window view over residue embeddings and windows.tsv rows."""

    def __init__(self, targets: list[tuple[str, int, int]], embeddings: dict[str, np.ndarray]):
        if not targets:
            raise ValueError("database windows.tsv contains no targets")
        self.targets = targets
        self.embeddings = embeddings
        self.window_dim = None
        for identifier, start, end in targets:
            if identifier not in embeddings:
                raise KeyError(f"database target {identifier!r} is missing from embeddings.npz")
            values = np.asarray(embeddings[identifier])
            if values.ndim != 2 or end <= start or end > values.shape[0]:
                raise ValueError(f"invalid target window {identifier}:{start}-{end}")
            dimension = int((end - start) * values.shape[1])
            if self.window_dim is None:
                self.window_dim = dimension
            elif dimension != self.window_dim:
                raise ValueError("database windows do not have a consistent dimension")
        self.shape = (len(targets), int(self.window_dim))

    def _window(self, ordinal: int) -> np.ndarray:
        identifier, start, end = self.targets[int(ordinal)]
        values = np.asarray(self.embeddings[identifier][start:end], dtype=np.float32).reshape(-1)
        norm = max(float(np.linalg.norm(values)), 1e-12)
        return np.ascontiguousarray(values / np.float32(norm), dtype=np.float32)

    def __getitem__(self, indices: object) -> np.ndarray:
        values = np.asarray(indices)
        if values.ndim == 0:
            return self._window(int(values))
        flat = values.reshape(-1)
        rows = np.stack([self._window(int(index)) for index in flat], axis=0)
        return rows.reshape(*values.shape, self.shape[1])


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def load_candidates(path: Path) -> list[CandidateGroup]:
    groups: dict[tuple[str, int, int], CandidateGroup] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        missing = [column for column in SEED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        for row in reader:
            key = (row["query_id"], int(row["query_start"]), int(row["query_end"]))
            group = groups.setdefault(key, CandidateGroup(*key, []))
            group.targets.append((row["target_id"], int(row["target_start"]), int(row["target_end"])))
    return list(groups.values())


def load_query_windows(
    window_paths: Iterable[Path],
    manifest_paths: Iterable[Path],
) -> dict[tuple[str, int, int], np.ndarray]:
    manifests = {path.name.replace(".windows.manifest.json", ""): load_json(path) for path in manifest_paths}
    result: dict[tuple[str, int, int], np.ndarray] = {}
    paths = list(window_paths)
    if not paths:
        raise ValueError("no query window NPZ files were supplied")
    for path in paths:
        stem = path.name.replace(".windows.npz", "")
        manifest = manifests.get(stem)
        if manifest is None:
            raise ValueError(f"no query manifest matches {path.name}")
        stride = int(manifest["stride"])
        window_size = int(manifest["window_size"])
        with np.load(path) as archive:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                if identifier not in archive.files:
                    raise KeyError(f"{identifier} is missing from {path}")
                values = np.asarray(archive[identifier], dtype=np.float32)
                for offset in range(values.shape[0]):
                    key = (identifier, offset * stride, offset * stride + window_size)
                    if key in result:
                        raise ValueError(f"duplicate query window {key}")
                    result[key] = np.ascontiguousarray(values[offset], dtype=np.float32)
    return result


def load_targets(path: Path) -> list[tuple[str, int, int]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"faiss_id", "transcript_id", "start", "end"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain {sorted(required)}")
        rows: list[tuple[str, int, int] | None] = []
        for row in reader:
            ordinal = int(row["faiss_id"])
            while len(rows) <= ordinal:
                rows.append(None)
            rows[ordinal] = (row["transcript_id"], int(row["start"]), int(row["end"]))
    if not rows or any(row is None for row in rows):
        raise ValueError(f"{path} has missing or empty faiss_id rows")
    return [row for row in rows if row is not None]


def load_embeddings(path: Path) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    with np.load(path) as archive:
        for identifier in archive.files:
            result[identifier] = np.asarray(archive[identifier])
    return result


def build_label_matrix(
    groups: list[CandidateGroup],
    query_vectors: dict[tuple[str, int, int], np.ndarray],
    target_index: dict[tuple[str, int, int], int],
) -> tuple[np.ndarray, np.ndarray, int]:
    if not groups:
        raise ValueError("candidate table contains no rows")
    query_rows: list[np.ndarray] = []
    max_candidates = max(len(group.targets) for group in groups)
    labels = np.full((len(groups), max_candidates), -1, dtype=np.int64)
    for row, group in enumerate(groups):
        query_key = (group.query_id, group.query_start, group.query_end)
        if query_key not in query_vectors:
            raise KeyError(f"query window {query_key} is missing from the query windows")
        query_rows.append(query_vectors[query_key])
        for column, target in enumerate(group.targets):
            if target not in target_index:
                raise KeyError(f"candidate target window {target} is missing from windows.tsv")
            labels[row, column] = target_index[target]
    queries = np.ascontiguousarray(np.stack(query_rows), dtype=np.float32)
    return queries, labels, max_candidates


def write_seeds(
    path: Path,
    groups: list[CandidateGroup],
    targets: list[tuple[str, int, int]],
    labels: np.ndarray,
    scores: np.ndarray,
    min_similarity: float,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SEED_COLUMNS)
        for group, label_row, score_row in zip(groups, labels, scores):
            rank = 0
            for label, score in zip(label_row, score_row):
                if int(label) < 0 or not np.isfinite(score) or float(score) < min_similarity:
                    continue
                rank += 1
                target_id, target_start, target_end = targets[int(label)]
                writer.writerow(
                    [
                        group.query_id,
                        group.query_start,
                        group.query_end,
                        target_id,
                        target_start,
                        target_end,
                        f"{float(score):.6f}",
                        rank,
                    ]
                )
                rows += 1
    return rows


def rerank(args: argparse.Namespace) -> dict:
    groups = load_candidates(args.candidates)
    targets = load_targets(args.database / "windows.tsv")
    target_index = {target: index for index, target in enumerate(targets)}
    embeddings = load_embeddings(args.database / "embeddings.npz")
    database_windows = OriginalWindowStore(targets, embeddings)
    query_vectors = load_query_windows(args.query_windows, args.query_manifests)
    queries, labels, max_candidates = build_label_matrix(groups, query_vectors, target_index)
    if args.output_k < 1 or args.output_k > max_candidates:
        raise ValueError("output-k must be between 1 and the largest candidate pool")
    final_labels, final_scores, elapsed = rerank_label_matrix(
        database_windows,
        queries,
        labels,
        args.output_k,
        batch_size=args.batch_size,
        candidate_batch_size=args.candidate_batch_size,
        workers=args.workers,
        device=args.device,
    )
    n_seeds = write_seeds(
        args.output,
        groups,
        targets,
        final_labels,
        final_scores,
        args.min_similarity,
    )
    metrics = {
        "backend": "exact_original_window_rerank",
        "device": args.device,
        "workers": int(args.workers),
        "batch_size": int(args.batch_size),
        "candidate_batch_size": int(args.candidate_batch_size),
        "candidate_pool_max": int(max_candidates),
        "output_k": int(args.output_k),
        "n_query_windows": len(groups),
        "n_database_windows": len(targets),
        "n_seeds": n_seeds,
        "window_dim": int(database_windows.shape[1]),
        "exact_rerank_seconds": elapsed,
        "query_windows_per_second": len(groups) / elapsed if elapsed > 0 else None,
    }
    if args.metrics:
        args.metrics.parent.mkdir(parents=True, exist_ok=True)
        args.metrics.write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--query-windows", type=Path, nargs="+", required=True)
    parser.add_argument("--query-manifests", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path)
    parser.add_argument("--output-k", type=int, default=50)
    parser.add_argument("--min-similarity", type=float, default=0.8)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--candidate-batch-size", type=int, default=2048)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="cpu")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        metrics = rerank(args)
    except (OSError, KeyError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(metrics))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
