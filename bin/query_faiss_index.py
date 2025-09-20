#!/usr/bin/env python3
"""Query FAISS index with window vectors and emit BLAST-style seeds."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple

import faiss
import numpy as np
import pandas as pd


SEED_COLUMNS = [
    "query_transcript",
    "target_transcript",
    "query_window_start",
    "query_window_end",
    "target_window_start",
    "target_window_end",
    "query_window_index",
    "target_window_index",
    "similarity",
    "rank",
    "diagonal",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Query FAISS index with query windows and record seed hits")
    parser.add_argument("--index", required=True, help="Path to FAISS index built over database windows")
    parser.add_argument("--database-vectors", required=True, help=".npy array with database window vectors")
    parser.add_argument("--database-metadata", required=True, help="TSV with metadata for database windows")
    parser.add_argument("--query-vectors", required=True, help=".npy array with query window vectors")
    parser.add_argument("--query-metadata", required=True, help="TSV with metadata for query windows")
    parser.add_argument("--id-column", default="transcript_id", help="Column identifying transcripts")
    parser.add_argument("--top-k", type=int, default=50, help="Neighbours per query window to evaluate")
    parser.add_argument("--similarity-threshold", type=float, default=0.7, help="Minimum cosine similarity to keep a seed")
    parser.add_argument("--metric", default="ip", choices=["ip", "l2"], help="Metric used by the FAISS index")
    parser.add_argument("--output", required=True, help="Output TSV for accepted seeds")
    parser.add_argument("--stats-json", default=None, help="Optional path to write query statistics JSON")
    parser.add_argument("--nprobe", type=int, default=None, help="Optional nprobe override for IVF indices")
    parser.add_argument("--use-gpu", action="store_true", help="Move index to GPU for querying when available")
    return parser.parse_args()


def load_index(path: str, use_gpu: bool, nprobe: int | None = None) -> faiss.Index:
    index = faiss.read_index(path)
    if nprobe is not None and hasattr(index, "nprobe"):
        index.nprobe = nprobe
    if not use_gpu:
        return index
    if not hasattr(faiss, "StandardGpuResources"):
        print("WARNING: FAISS build lacks GPU support; continuing on CPU")
        return index
    res = faiss.StandardGpuResources()
    return faiss.index_cpu_to_gpu(res, 0, index)


def cosine_from_metric(metric: str, distances: np.ndarray) -> np.ndarray:
    if metric == "ip":
        return np.clip(distances, -1.0, 1.0)
    # FAISS L2 returns squared distance. With unit-normalized vectors:
    # d^2 = ||q||^2 + ||x||^2 - 2 q·x = 2 - 2 cosine
    # ⇒ cosine = 1 - d^2 / 2
    cos = 1.0 - (distances / 2.0)
    return np.clip(cos, -1.0, 1.0)


def main() -> None:
    args = parse_args()
    # Database vectors are only needed to validate metadata counts; memory-map them
    # to avoid pulling multi-GB arrays into RAM when the database is large.
    db_vectors = np.load(args.database_vectors, mmap_mode="r")
    db_vector_count = db_vectors.shape[0]
    del db_vectors

    query_vectors = np.load(args.query_vectors)

    db_meta = pd.read_csv(args.database_metadata, sep="\t")
    query_meta = pd.read_csv(args.query_metadata, sep="\t")
    id_col = args.id_column

    if len(query_meta) != len(query_vectors):
        raise SystemExit("Mismatch between query metadata rows and vector count")

    if len(db_meta) != db_vector_count:
        raise SystemExit("Mismatch between database metadata rows and vector count")

    if query_vectors.size == 0:
        pd.DataFrame(columns=SEED_COLUMNS).to_csv(args.output, sep="\t", index=False)
        return

    index = load_index(args.index, args.use_gpu, args.nprobe)
    distances, indices = index.search(query_vectors, args.top_k)

    rows = []
    kept = 0
    for q_idx, meta_row in enumerate(query_meta.to_dict("records")):
        dists = distances[q_idx]
        neighs = indices[q_idx]
        similarities = cosine_from_metric(args.metric, dists)
        for rank, (db_idx, sim) in enumerate(zip(neighs, similarities), start=1):
            if db_idx < 0:
                continue
            if sim < args.similarity_threshold:
                continue
            db_row = db_meta.iloc[int(db_idx)]
            if db_row[id_col] == meta_row[id_col]:
                continue
            rows.append({
                "query_transcript": str(meta_row[id_col]),
                "target_transcript": str(db_row[id_col]),
                "query_window_start": int(meta_row["window_start"]),
                "query_window_end": int(meta_row["window_end"]),
                "target_window_start": int(db_row["window_start"]),
                "target_window_end": int(db_row["window_end"]),
                "query_window_index": int(meta_row["window_index"]),
                "target_window_index": int(db_row["window_index"]),
                "similarity": float(sim),
                "rank": rank,
                "diagonal": int(db_row["window_start"]) - int(meta_row["window_start"]),
            })
            kept += 1

    out_df = pd.DataFrame(rows, columns=SEED_COLUMNS)
    out_df.to_csv(args.output, sep="\t", index=False)

    stats = {
        "query_windows": len(query_vectors),
        "seeds_kept": kept,
        "top_k": args.top_k,
        "similarity_threshold": args.similarity_threshold,
    }
    if args.stats_json:
        Path(args.stats_json).write_text(json.dumps(stats, indent=2) + "\n")

    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
