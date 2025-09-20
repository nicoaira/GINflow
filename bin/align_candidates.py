#!/usr/bin/env python3
"""Run banded Smith–Waterman alignments over clustered seeds."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

NEG_INF = -1e9


@dataclass
class AlignmentResult:
    query_id: str
    target_id: str
    cluster_id: int
    score: float
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    alignment_length: int
    avg_cosine: float
    gap_opens: int
    gap_bases: int
    seed_count: int
    max_seed_similarity: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align clustered BLAST-style seeds with affine-gap Smith–Waterman")
    parser.add_argument("--cluster-members", required=True, help="TSV mapping seeds to cluster IDs")
    parser.add_argument("--cluster-summaries", required=True, help="TSV summarising clusters")
    parser.add_argument("--embeddings", required=True, help="Node-level embeddings TSV")
    parser.add_argument("--id-column", default="transcript_id", help="Identifier column in embeddings/meta tables")
    parser.add_argument("--node-index-column", default="node_index", help="Column giving node order within each transcript")
    parser.add_argument("--embedding-column", default="embedding_vector", help="Column containing comma-separated embedding vectors")
    parser.add_argument("--meta", required=True, help="TSV with transcript metadata (sequence, structure, etc.)")
    parser.add_argument("--sequence-column", default=None, help="Sequence column to include in outputs (auto-detected if omitted)")
    parser.add_argument("--structure-column", default=None, help="Structure column to include in outputs (auto-detected if omitted)")
    parser.add_argument("--background-samples", type=int, default=10000, help="Number of random node pairs for μ₀/σ₀ estimation")
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed for background sampling")
    parser.add_argument("--gamma", type=float, default=1.5, help="Scaling factor for similarity scores")
    parser.add_argument("--band-width", type=int, default=96, help="Diagonal band width for SW (in nucleotides)")
    parser.add_argument("--xdrop", type=float, default=50.0, help="Drop threshold for terminating extension")
    parser.add_argument("--gap-open", type=float, default=12.0, help="Gap opening penalty")
    parser.add_argument("--gap-extend", type=float, default=2.0, help="Gap extension penalty")
    parser.add_argument("--padding", type=int, default=32, help="Extra nucleotides to extend around seed bounds before alignment")
    parser.add_argument("--score-min", type=float, default=-4.0, help="Minimum clamped per-position score")
    parser.add_argument("--score-max", type=float, default=8.0, help="Maximum clamped per-position score")
    parser.add_argument("--top-n", type=int, default=50, help="Report at most this many alignments globally")
    parser.add_argument("--output", required=True, help="Output TSV for alignment summaries")
    parser.add_argument("--stats-json", default=None, help="Optional JSON stats path")
    return parser.parse_args()


def load_embeddings(path: str, id_col: str, node_idx_col: str, emb_col: str) -> Dict[str, np.ndarray]:
    df = pd.read_csv(path, sep="\t", dtype={id_col: str})
    required = {id_col, node_idx_col, emb_col}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Embeddings file missing columns: {', '.join(sorted(missing))}")
    df = df.rename(columns={node_idx_col: "node_index", emb_col: "embedding_vector"})
    embeddings: Dict[str, np.ndarray] = {}
    for tid, group in df.groupby(id_col):
        ordered = group.sort_values("node_index")
        vectors = np.stack(ordered["embedding_vector"].astype(str).str.split(',').map(lambda xs: np.asarray(xs, dtype=np.float32)))
        embeddings[tid] = vectors
    return embeddings


def sample_background(embeddings: Dict[str, np.ndarray], sample_count: int, seed: int) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    ids = [tid for tid, vec in embeddings.items() if len(vec) > 0]
    if not ids:
        raise SystemExit("No embeddings available for background sampling")
    sims: List[float] = []
    for _ in range(sample_count):
        qid = rng.choice(ids)
        tid = rng.choice(ids)
        qvec = embeddings[qid][rng.integers(len(embeddings[qid]))]
        tvec = embeddings[tid][rng.integers(len(embeddings[tid]))]
        sims.append(float(np.dot(qvec, tvec)))
    mu = float(np.mean(sims))
    sigma = float(np.std(sims))
    if sigma < 1e-6:
        sigma = 1.0
    return mu, sigma


def compute_score(dot: float, mu: float, sigma: float, gamma: float, smin: float, smax: float) -> float:
    z = gamma * ((dot - mu) / sigma)
    return float(np.clip(z, smin, smax))


def smith_waterman(
    query: np.ndarray,
    target: np.ndarray,
    mu: float,
    sigma: float,
    gamma: float,
    band: int,
    diag_offset: int,
    gap_open: float,
    gap_extend: float,
    xdrop: float,
    clamp_min: float,
    clamp_max: float,
) -> Tuple[float, List[Tuple[str, int, int]]]:
    m, n = len(query), len(target)
    if m == 0 or n == 0:
        return 0.0, []

    band_half = band // 2

    M_prev = np.zeros(n + 1, dtype=np.float32)
    X_prev = np.full(n + 1, NEG_INF, dtype=np.float32)
    Y_prev = np.full(n + 1, NEG_INF, dtype=np.float32)
    M_curr = np.zeros(n + 1, dtype=np.float32)
    X_curr = np.full(n + 1, NEG_INF, dtype=np.float32)
    Y_curr = np.full(n + 1, NEG_INF, dtype=np.float32)

    trace_M: Dict[Tuple[int, int], Optional[Tuple[str, int, int]]] = {}
    trace_X: Dict[Tuple[int, int], Optional[Tuple[str, int, int]]] = {}
    trace_Y: Dict[Tuple[int, int], Optional[Tuple[str, int, int]]] = {}

    best_score = 0.0
    best_state = ("M", 0, 0)

    for i in range(1, m + 1):
        M_curr.fill(NEG_INF)
        X_curr.fill(NEG_INF)
        Y_curr.fill(NEG_INF)
        M_curr[0] = 0.0

        row_best = NEG_INF
        center = i + diag_offset
        j_start = max(1, center - band_half)
        j_end = min(n, center + band_half)
        if j_start > j_end:
            continue

        for j in range(j_start, j_end + 1):
            i_global = i - 1
            j_global = j - 1
            dot = float(np.dot(query[i_global], target[j_global]))
            score = compute_score(dot, mu, sigma, gamma, clamp_min, clamp_max)

            # Match/mismatch state
            candidates = [
                (0.0, None),
                (float(M_prev[j - 1]), ("M", i - 1, j - 1)),
                (float(X_prev[j - 1]), ("X", i - 1, j - 1)),
                (float(Y_prev[j - 1]), ("Y", i - 1, j - 1)),
            ]
            base_val, base_ptr = max(candidates, key=lambda item: item[0])
            match_val = base_val + score
            M_curr[j] = match_val
            trace_M[(i, j)] = base_ptr

            # Gap in query (insertion in target)
            left_m = M_curr[j - 1] - gap_open
            left_x = X_curr[j - 1] - gap_extend
            if left_m >= left_x:
                X_curr[j] = left_m
                trace_X[(i, j)] = ("M", i, j - 1)
            else:
                X_curr[j] = left_x
                trace_X[(i, j)] = ("X", i, j - 1)

            # Gap in target (deletion)
            up_m = M_prev[j] - gap_open
            up_y = Y_prev[j] - gap_extend
            if up_m >= up_y:
                Y_curr[j] = up_m
                trace_Y[(i, j)] = ("M", i - 1, j)
            else:
                Y_curr[j] = up_y
                trace_Y[(i, j)] = ("Y", i - 1, j)

            cell_best = max(match_val, X_curr[j], Y_curr[j])
            if cell_best > best_score:
                if cell_best == match_val:
                    best_state = ("M", i, j)
                elif cell_best == X_curr[j]:
                    best_state = ("X", i, j)
                else:
                    best_state = ("Y", i, j)
                best_score = cell_best

            row_best = max(row_best, cell_best)

        if best_score - row_best > xdrop:
            break

        M_prev, M_curr = M_curr, M_prev
        X_prev, X_curr = X_curr, X_prev
        Y_prev, Y_curr = Y_curr, Y_prev
        M_prev[0] = 0.0

    path = traceback(best_state, trace_M, trace_X, trace_Y)
    return best_score, path


def traceback(
    best_state: Tuple[str, int, int],
    trace_M: Dict[Tuple[int, int], Optional[Tuple[str, int, int]]],
    trace_X: Dict[Tuple[int, int], Optional[Tuple[str, int, int]]],
    trace_Y: Dict[Tuple[int, int], Optional[Tuple[str, int, int]]],
) -> List[Tuple[str, int, int]]:
    state, i, j = best_state
    if i <= 0 or j <= 0:
        return []
    path: List[Tuple[str, int, int]] = []
    while state in {"M", "X", "Y"} and i >= 0 and j >= 0:
        path.append((state, i, j))
        if state == "M":
            ptr = trace_M.get((i, j))
            if ptr is None:
                break
            state, i, j = ptr
        elif state == "X":
            ptr = trace_X.get((i, j))
            if ptr is None:
                break
            state, i, j = ptr
        else:  # Y
            ptr = trace_Y.get((i, j))
            if ptr is None:
                break
            state, i, j = ptr
    path.reverse()
    return path


def attach_metadata(
    result: AlignmentResult,
    meta_map: Dict[str, dict],
    seq_col: Optional[str],
    struct_col: Optional[str],
) -> dict:
    record = result.__dict__.copy()
    q_meta = meta_map.get(result.query_id, {})
    t_meta = meta_map.get(result.target_id, {})

    def _slice(text: Optional[str], start: int, end: int) -> str:
        if text is None or pd.isna(text):
            return ""
        return str(text)[start:end]

    if seq_col:
        record["query_sequence"] = q_meta.get(seq_col, "")
        record["target_sequence"] = t_meta.get(seq_col, "")
        record["query_sequence_region"] = _slice(record["query_sequence"], result.query_start, result.query_end)
        record["target_sequence_region"] = _slice(record["target_sequence"], result.target_start, result.target_end)
    if struct_col:
        record["query_structure"] = q_meta.get(struct_col, "")
        record["target_structure"] = t_meta.get(struct_col, "")
        record["query_structure_region"] = _slice(record.get("query_structure"), result.query_start, result.query_end)
        record["target_structure_region"] = _slice(record.get("target_structure"), result.target_start, result.target_end)
    return record


def infer_columns(meta: pd.DataFrame, seq_col: Optional[str], struct_col: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    if seq_col and seq_col in meta.columns:
        seq = seq_col
    else:
        candidates = [c for c in meta.columns if "sequence" in c.lower()]
        seq = candidates[0] if candidates else None
    if struct_col and struct_col in meta.columns:
        struct = struct_col
    else:
        candidates = [c for c in meta.columns if "struct" in c.lower()]
        struct = candidates[0] if candidates else None
    return seq, struct


def main() -> None:
    args = parse_args()

    clusters_df = pd.read_csv(args.cluster_summaries, sep="\t")
    members_df = pd.read_csv(args.cluster_members, sep="\t")
    if clusters_df.empty or members_df.empty:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        if args.stats_json:
            Path(args.stats_json).write_text(json.dumps({"alignments": 0}, indent=2) + "\n")
        return

    embeddings = load_embeddings(args.embeddings, args.id_column, args.node_index_column, args.embedding_column)
    mu0, sigma0 = sample_background(embeddings, args.background_samples, args.random_seed)

    meta_df = pd.read_csv(args.meta, sep="\t", dtype={args.id_column: str})
    meta_df = meta_df.drop_duplicates(subset=[args.id_column])
    meta_map = meta_df.set_index(args.id_column).to_dict(orient="index")
    seq_col, struct_col = infer_columns(meta_df, args.sequence_column, args.structure_column)

    results: List[dict] = []

    grouped_members = members_df.groupby("cluster_id")
    for cluster in clusters_df.to_dict("records"):
        cid = cluster["cluster_id"]
        if cid not in grouped_members.groups:
            continue
        seeds = grouped_members.get_group(cid)
        query_id = str(seeds["query_transcript"].iloc[0])
        target_id = str(seeds["target_transcript"].iloc[0])
        diag = int(round(seeds["diagonal"].median()))

        query_vecs = embeddings.get(query_id)
        target_vecs = embeddings.get(target_id)
        if query_vecs is None or target_vecs is None:
            continue

        q_start = max(0, int(seeds["query_window_start"].min()) - args.padding)
        q_end = min(len(query_vecs), int(seeds["query_window_end"].max()) + args.padding)
        t_start = max(0, int(seeds["target_window_start"].min()) - args.padding)
        t_end = min(len(target_vecs), int(seeds["target_window_end"].max()) + args.padding)

        query_slice = query_vecs[q_start:q_end]
        target_slice = target_vecs[t_start:t_end]
        diag_offset = diag + q_start - t_start

        score, path = smith_waterman(
            query_slice,
            target_slice,
            mu0,
            sigma0,
            args.gamma,
            args.band_width,
            diag_offset,
            args.gap_open,
            args.gap_extend,
            args.xdrop,
            args.score_min,
            args.score_max,
        )
        if not path or score <= 0:
            continue

        alignment_info = summarise_alignment(
            path,
            query_slice,
            target_slice,
            q_start,
            t_start,
        )

        result = AlignmentResult(
            query_id=query_id,
            target_id=target_id,
            cluster_id=cid,
            score=float(score),
            query_start=alignment_info["query_start"],
            query_end=alignment_info["query_end"],
            target_start=alignment_info["target_start"],
            target_end=alignment_info["target_end"],
            alignment_length=alignment_info["length"],
            avg_cosine=alignment_info["avg_cosine"],
            gap_opens=alignment_info["gap_opens"],
            gap_bases=alignment_info["gap_bases"],
            seed_count=int(cluster["seed_count"]),
            max_seed_similarity=float(cluster["max_similarity"]),
        )
        results.append(attach_metadata(result, meta_map, seq_col, struct_col))

    if not results:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        if args.stats_json:
            Path(args.stats_json).write_text(json.dumps({"alignments": 0}, indent=2) + "\n")
        return

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("score", ascending=False).head(args.top_n)
    results_df.to_csv(args.output, sep="\t", index=False)

    stats = {"alignments": int(len(results_df)), "mu0": mu0, "sigma0": sigma0}
    if args.stats_json:
        Path(args.stats_json).write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


def summarise_alignment(
    path: List[Tuple[str, int, int]],
    query: np.ndarray,
    target: np.ndarray,
    q_offset: int,
    t_offset: int,
) -> dict:
    q_cursor = q_offset
    t_cursor = t_offset
    query_positions: List[int] = []
    target_positions: List[int] = []
    cosines: List[float] = []
    gap_opens = 0
    gap_bases = 0
    in_gap_q = False
    in_gap_t = False

    for state, i, j in path:
        if state == "M":
            query_positions.append(q_cursor)
            target_positions.append(t_cursor)
            cosines.append(float(np.dot(query[q_cursor - q_offset], target[t_cursor - t_offset])))
            q_cursor += 1
            t_cursor += 1
            in_gap_q = False
            in_gap_t = False
        elif state == "X":  # gap in query (target advances)
            if not in_gap_q:
                gap_opens += 1
                in_gap_q = True
            gap_bases += 1
            t_cursor += 1
            in_gap_t = False
        elif state == "Y":  # gap in target (query advances)
            if not in_gap_t:
                gap_opens += 1
                in_gap_t = True
            gap_bases += 1
            q_cursor += 1
            in_gap_q = False

    if not query_positions or not target_positions:
        return {
            "query_start": q_offset,
            "query_end": q_offset,
            "target_start": t_offset,
            "target_end": t_offset,
            "length": 0,
            "avg_cosine": 0.0,
            "gap_opens": gap_opens,
            "gap_bases": gap_bases,
        }

    return {
        "query_start": int(min(query_positions)),
        "query_end": int(max(query_positions) + 1),
        "target_start": int(min(target_positions)),
        "target_end": int(max(target_positions) + 1),
        "length": len(path),
        "avg_cosine": float(np.mean(cosines) if cosines else 0.0),
        "gap_opens": gap_opens,
        "gap_bases": gap_bases,
    }


if __name__ == "__main__":
    main()
