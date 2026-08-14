#!/usr/bin/env python3
"""Cluster FAISS seeds into HSP-like diagonal groups."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


SEED_COLUMNS = [
    "query_id",
    "query_start",
    "query_end",
    "target_id",
    "target_start",
    "target_end",
    "score",
    "rank",
]

CLUSTER_COLUMNS = [
    "cluster_id",
    "query_id",
    "target_id",
    "seed_count",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "max_score",
    "diagonal",
    "diagonal_min",
    "diagonal_max",
    "diagonal_span",
]


def load_seeds(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        missing = [name for name in SEED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(missing)}")
        return list(reader)


def cluster_pair(
    rows: list[dict],
    span: int,
    min_seeds: int,
    diag_tol: int,
    max_diag_span: int,
    start_id: int,
) -> tuple[list[dict], list[dict], int]:
    n = len(rows)
    q_start = np.fromiter((int(row["query_start"]) for row in rows), dtype=np.int32, count=n)
    q_end = np.fromiter((int(row["query_end"]) for row in rows), dtype=np.int32, count=n)
    t_start = np.fromiter((int(row["target_start"]) for row in rows), dtype=np.int32, count=n)
    t_end = np.fromiter((int(row["target_end"]) for row in rows), dtype=np.int32, count=n)
    score = np.fromiter((float(row["score"]) for row in rows), dtype=np.float64, count=n)
    diag = q_start.astype(np.int32) - t_start.astype(np.int32)
    used = np.zeros(n, dtype=bool)
    clusters: list[dict] = []
    members: list[dict] = []
    cluster_id = start_id

    for idx in np.argsort(-score, kind="stable"):
        seed = int(idx)
        if used[seed]:
            continue
        members_idx = [seed]
        used[seed] = True
        q0, q1 = int(q_start[seed]), int(q_end[seed])
        t0, t1 = int(t_start[seed]), int(t_end[seed])
        d0 = d1 = int(diag[seed])

        while True:
            cand = ~used
            if not cand.any():
                break
            cand &= (diag >= d0 - diag_tol) & (diag <= d1 + diag_tol)
            if max_diag_span > 0:
                cand &= (np.maximum(d1, diag) - np.minimum(d0, diag)) <= max_diag_span
            cand &= (q_end >= q0 - span) & (q_start <= q1 + span)
            cand &= (t_end >= t0 - span) & (t_start <= t1 + span)
            if not cand.any():
                break
            best = int(np.flatnonzero(cand)[np.argmax(score[cand])])
            used[best] = True
            members_idx.append(best)
            q0 = min(q0, int(q_start[best]))
            q1 = max(q1, int(q_end[best]))
            t0 = min(t0, int(t_start[best]))
            t1 = max(t1, int(t_end[best]))
            d0 = min(d0, int(diag[best]))
            d1 = max(d1, int(diag[best]))

        if len(members_idx) < min_seeds:
            continue
        diag_values = diag[members_idx]
        clusters.append({
            "cluster_id": cluster_id,
            "query_id": rows[members_idx[0]]["query_id"],
            "target_id": rows[members_idx[0]]["target_id"],
            "seed_count": len(members_idx),
            "query_start": q0,
            "query_end": q1,
            "target_start": t0,
            "target_end": t1,
            "max_score": float(np.max(score[members_idx])),
            "diagonal": int(np.median(diag_values)),
            "diagonal_min": d0,
            "diagonal_max": d1,
            "diagonal_span": d1 - d0,
        })
        for member in members_idx:
            members.append({"cluster_id": cluster_id, **rows[member], "diagonal": int(diag[member])})
        cluster_id += 1

    return clusters, members, cluster_id


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, required=True)
    parser.add_argument("--output-clusters", type=Path, required=True)
    parser.add_argument("--output-members", type=Path, required=True)
    parser.add_argument("--cluster-span", type=int, default=80)
    parser.add_argument("--min-seeds", type=int, default=2)
    parser.add_argument("--diagonal-tolerance", type=int, default=12)
    parser.add_argument("--max-diagonal-span", type=int, default=96)
    parser.add_argument("--max-seed-rank", type=int, default=10)
    parser.add_argument("--stats-json", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.cluster_span < 0 or args.min_seeds < 1 or args.diagonal_tolerance < 0:
        print("error: cluster-span, min-seeds, and diagonal-tolerance must be valid", file=sys.stderr)
        return 2
    try:
        seeds = load_seeds(args.seeds)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    before = len(seeds)
    if args.max_seed_rank > 0:
        seeds = [row for row in seeds if int(row["rank"]) <= args.max_seed_rank]
    after = len(seeds)

    groups: dict[tuple[str, str], list[dict]] = {}
    for row in seeds:
        groups.setdefault((row["query_id"], row["target_id"]), []).append(row)

    clusters: list[dict] = []
    members: list[dict] = []
    next_id = 0
    for key in sorted(groups):
        cls, mem, next_id = cluster_pair(
            groups[key],
            span=args.cluster_span,
            min_seeds=args.min_seeds,
            diag_tol=args.diagonal_tolerance,
            max_diag_span=args.max_diagonal_span,
            start_id=next_id,
        )
        clusters.extend(cls)
        members.extend(mem)

    write_tsv(args.output_clusters, clusters, CLUSTER_COLUMNS)
    write_tsv(args.output_members, members, ["cluster_id", *SEED_COLUMNS, "diagonal"])
    stats = {
        "clusters": len(clusters),
        "clustered_seeds": len(members),
        "seeds_before_rank_filter": before,
        "seeds_after_rank_filter": after,
    }
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
