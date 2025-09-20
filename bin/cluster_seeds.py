#!/usr/bin/env python3
"""Cluster BLAST-style seeds along shared diagonals."""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List

import pandas as pd


@dataclass
class Seed:
    query_transcript: str
    target_transcript: str
    query_window_start: int
    query_window_end: int
    target_window_start: int
    target_window_end: int
    query_window_index: int
    target_window_index: int
    similarity: float
    rank: int
    diagonal: int


@dataclass
class Cluster:
    cluster_id: int
    query_transcript: str
    target_transcript: str
    diagonal: int
    seed_count: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    max_similarity: float


CLUSTER_COLUMNS = [
    "cluster_id",
    "query_transcript",
    "target_transcript",
    "diagonal",
    "seed_count",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "max_similarity",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Group seeds into diagonal clusters")
    parser.add_argument("--seeds", required=True, help="TSV emitted by query_faiss_index.py")
    parser.add_argument("--output-clusters", required=True, help="TSV with cluster summaries")
    parser.add_argument("--output-members", required=True, help="TSV mapping seeds to cluster IDs")
    parser.add_argument("--cluster-span", type=int, default=80, help="Maximum nucleotide distance between neighboring seeds in a cluster")
    parser.add_argument("--min-seeds", type=int, default=2, help="Minimum number of seeds to keep a cluster")
    parser.add_argument("--stats-json", default=None, help="Optional path to write clustering stats")
    return parser.parse_args()


def load_seeds(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return df
    required = {
        "query_transcript",
        "target_transcript",
        "query_window_start",
        "query_window_end",
        "target_window_start",
        "target_window_end",
        "similarity",
        "diagonal",
    }
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Seeds file missing columns: {', '.join(sorted(missing))}")
    return df


def cluster_group(group: pd.DataFrame, span: int, min_seeds: int, start_cluster: int) -> tuple[List[Cluster], List[dict], int]:
    clusters: List[Cluster] = []
    members: List[dict] = []
    cluster_id = start_cluster

    # Sort by query position to apply sliding window grouping
    ordered = group.sort_values("query_window_start")
    active: List[pd.Series] = []

    for row in ordered.itertuples(index=False):
        if not active:
            active.append(row)
            continue
        last = active[-1]
        if (
            row.query_window_start - last.query_window_start <= span
            and row.target_window_start - last.target_window_start <= span
        ):
            active.append(row)
        else:
            if len(active) >= min_seeds:
                clusters.append(build_cluster(active, cluster_id))
                members.extend(
                    {"cluster_id": cluster_id, **seed_as_dict(s)}
                    for s in active
                )
                cluster_id += 1
            active = [row]

    if active and len(active) >= min_seeds:
        clusters.append(build_cluster(active, cluster_id))
        members.extend({"cluster_id": cluster_id, **seed_as_dict(s)} for s in active)
        cluster_id += 1

    return clusters, members, cluster_id


def seed_as_dict(seed: pd.Series) -> dict:
    return {
        "query_transcript": seed.query_transcript,
        "target_transcript": seed.target_transcript,
        "query_window_start": int(seed.query_window_start),
        "query_window_end": int(seed.query_window_end),
        "target_window_start": int(seed.target_window_start),
        "target_window_end": int(seed.target_window_end),
        "query_window_index": int(seed.query_window_index),
        "target_window_index": int(seed.target_window_index),
        "similarity": float(seed.similarity),
        "rank": int(seed.rank),
        "diagonal": int(seed.diagonal),
    }


def build_cluster(seeds: List[pd.Series], cluster_id: int) -> Cluster:
    query_start = min(seed.query_window_start for seed in seeds)
    query_end = max(seed.query_window_end for seed in seeds)
    target_start = min(seed.target_window_start for seed in seeds)
    target_end = max(seed.target_window_end for seed in seeds)
    max_sim = max(seed.similarity for seed in seeds)
    seed0 = seeds[0]
    return Cluster(
        cluster_id=cluster_id,
        query_transcript=seed0.query_transcript,
        target_transcript=seed0.target_transcript,
        diagonal=seed0.diagonal,
        seed_count=len(seeds),
        query_start=int(query_start),
        query_end=int(query_end),
        target_start=int(target_start),
        target_end=int(target_end),
        max_similarity=float(max_sim),
    )


def main() -> None:
    args = parse_args()
    seeds_df = load_seeds(args.seeds)

    if seeds_df.empty:
        pd.DataFrame(columns=CLUSTER_COLUMNS).to_csv(args.output_clusters, sep="\t", index=False)
        pd.DataFrame(columns=["cluster_id"] + SEED_COLUMNS).to_csv(args.output_members, sep="\t", index=False)
        if args.stats_json:
            Path(args.stats_json).write_text(json.dumps({"clusters": 0, "seeds": 0}, indent=2) + "\n")
        return

    clusters: List[Cluster] = []
    members: List[dict] = []
    next_cluster = 0

    for (_, _, diag), group in seeds_df.groupby([
        "query_transcript",
        "target_transcript",
        "diagonal",
    ]):
        cls, mem, next_cluster = cluster_group(group, span=args.cluster_span, min_seeds=args.min_seeds, start_cluster=next_cluster)
        clusters.extend(cls)
        members.extend(mem)

    clusters_df = pd.DataFrame([asdict(c) for c in clusters])
    members_df = pd.DataFrame(members)

    clusters_df.to_csv(args.output_clusters, sep="\t", index=False)
    members_df.to_csv(args.output_members, sep="\t", index=False)

    stats = {"clusters": len(clusters_df), "seeds": len(members_df)}
    if args.stats_json:
        Path(args.stats_json).write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
