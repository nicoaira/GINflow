#!/usr/bin/env python3
"""Run banded Smith–Waterman alignments over clustered seeds."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    parser.add_argument(
        "--band-buffer",
        type=int,
        default=32,
        help="Additional allowance (nt) added to the observed diagonal span when adapting the band width",
    )
    parser.add_argument(
        "--band-max-width",
        type=int,
        default=0,
        help="Maximum permitted band width after adaptation; 0 leaves it unbounded",
    )
    parser.add_argument("--xdrop", type=float, default=50.0, help="Drop threshold for terminating extension")
    parser.add_argument("--gap-open", type=float, default=12.0, help="Gap opening penalty")
    parser.add_argument("--gap-extend", type=float, default=2.0, help="Gap extension penalty")
    parser.add_argument("--padding", type=int, default=32, help="Extra nucleotides to extend around seed bounds before alignment")
    parser.add_argument("--score-min", type=float, default=-4.0, help="Minimum clamped per-position score")
    parser.add_argument("--score-max", type=float, default=8.0, help="Maximum clamped per-position score")
    parser.add_argument("--top-n", type=int, default=50, help="Report at most this many alignments globally")
    parser.add_argument("--output", required=True, help="Output TSV for alignment summaries")
    parser.add_argument("--stats-json", default=None, help="Optional JSON stats path")
    parser.add_argument("--dp-output", default=None, help="Optional JSONL file capturing DP traces for reported alignments")
    parser.add_argument("--alignment-text", default=None, help="Optional plain-text alignment dump (BLAST-style)")
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
    dp_records: List[dict] = []
    collect_trace = args.dp_output is not None

    for cluster in clusters_df.to_dict("records"):
        cid = cluster["cluster_id"]
        if cid not in grouped_members.groups:
            continue
        seeds = grouped_members.get_group(cid)
        query_id = str(seeds["query_transcript"].iloc[0])
        target_id = str(seeds["target_transcript"].iloc[0])
        diag_values = seeds["diagonal"].astype(int)
        diag_median = int(np.median(diag_values))
        diag_min = int(diag_values.min())
        diag_max = int(diag_values.max())
        diag_span = diag_max - diag_min

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
        diag_offset = diag_median + q_start - t_start

        effective_band = args.band_width
        adaptive_target = diag_span + args.band_buffer
        if adaptive_target % 2 != 0:
            adaptive_target += 1
        if adaptive_target > effective_band:
            effective_band = adaptive_target
        if args.band_max_width and args.band_max_width > 0:
            effective_band = min(effective_band, args.band_max_width)
        effective_band = max(1, effective_band)

        score, path = smith_waterman(
            query_slice,
            target_slice,
            mu0,
            sigma0,
            args.gamma,
            effective_band,
            diag_offset,
            args.gap_open,
            args.gap_extend,
            args.xdrop,
            args.score_min,
            args.score_max,
        )
        if not path or score <= 0:
            continue

        alignment_info, trace = summarise_alignment(
            path,
            query_slice,
            target_slice,
            q_start,
            t_start,
            collect_trace=collect_trace,
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
        enriched = attach_metadata(result, meta_map, seq_col, struct_col)

        if enriched.get("query_structure") and enriched.get("target_structure"):
            q_alignment, t_alignment, alignment_mask = build_alignment_strings(
                path,
                enriched.get("query_structure"),
                enriched.get("target_structure"),
                q_start,
                t_start,
            )
            if q_alignment:
                enriched["aligned_query_structure"] = q_alignment
            if t_alignment:
                enriched["aligned_target_structure"] = t_alignment
            if alignment_mask:
                enriched["alignment_mask"] = alignment_mask
            enriched["alignment_content"] = "structure"
        elif enriched.get("query_sequence") and enriched.get("target_sequence"):
            q_alignment, t_alignment, alignment_mask = build_alignment_strings(
                path,
                enriched.get("query_sequence"),
                enriched.get("target_sequence"),
                q_start,
                t_start,
            )
            if q_alignment:
                enriched["aligned_query_sequence"] = q_alignment
            if t_alignment:
                enriched["aligned_target_sequence"] = t_alignment
            if alignment_mask:
                enriched["alignment_mask"] = alignment_mask
            enriched["alignment_content"] = "sequence"
        enriched.update(
            {
                "diagonal_median": diag_median,
                "diagonal_min": diag_min,
                "diagonal_max": diag_max,
                "diagonal_span": diag_span,
                "band_width": effective_band,
            }
        )
        results.append(enriched)

        if collect_trace and trace is not None:
            dp_records.append(
                {
                    "query_id": query_id,
                    "target_id": target_id,
                    "cluster_id": cid,
                    "score": float(score),
                    "alignment_start": {
                        "query": alignment_info["query_start"],
                        "target": alignment_info["target_start"],
                    },
                    "alignment_end": {
                        "query": alignment_info["query_end"],
                        "target": alignment_info["target_end"],
                    },
                    "diagonal": {
                        "median": diag_median,
                        "min": diag_min,
                        "max": diag_max,
                        "span": diag_span,
                        "offset": diag_offset,
                    },
                    "band_width": effective_band,
                    "trace": trace,
                }
            )

    if not results:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        if args.dp_output:
            Path(args.dp_output).write_text("")
        if args.alignment_text:
            Path(args.alignment_text).write_text("")
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

    if args.dp_output:
        with Path(args.dp_output).open("w", encoding="utf-8") as handle:
            for record in dp_records:
                handle.write(json.dumps(record) + "\n")
        if not dp_records:
            Path(args.dp_output).touch()

    if args.alignment_text:
        write_alignment_text(Path(args.alignment_text), results_df)


def summarise_alignment(
    path: List[Tuple[str, int, int]],
    query: np.ndarray,
    target: np.ndarray,
    q_offset: int,
    t_offset: int,
    collect_trace: bool = False,
) -> Tuple[dict, Optional[List[dict]]]:
    q_cursor = q_offset
    t_cursor = t_offset
    query_positions: List[int] = []
    target_positions: List[int] = []
    cosines: List[float] = []
    gap_opens = 0
    gap_bases = 0
    in_gap_q = False
    in_gap_t = False

    trace: Optional[List[dict]] = [] if collect_trace else None

    for state, i, j in path:
        if state == "M":
            query_positions.append(q_cursor)
            target_positions.append(t_cursor)
            cosines.append(float(np.dot(query[q_cursor - q_offset], target[t_cursor - t_offset])))
            if trace is not None:
                trace.append({
                    "state": "M",
                    "query_pos": q_cursor,
                    "target_pos": t_cursor,
                })
            q_cursor += 1
            t_cursor += 1
            in_gap_q = False
            in_gap_t = False
        elif state == "X":  # gap in query (target advances)
            if not in_gap_q:
                gap_opens += 1
                in_gap_q = True
            gap_bases += 1
            if trace is not None:
                trace.append({
                    "state": "X",
                    "query_pos": q_cursor,
                    "target_pos": t_cursor,
                })
            t_cursor += 1
            in_gap_t = False
        elif state == "Y":  # gap in target (query advances)
            if not in_gap_t:
                gap_opens += 1
                in_gap_t = True
            gap_bases += 1
            if trace is not None:
                trace.append({
                    "state": "Y",
                    "query_pos": q_cursor,
                    "target_pos": t_cursor,
                })
            q_cursor += 1
            in_gap_q = False

    if not query_positions or not target_positions:
        summary = {
            "query_start": q_offset,
            "query_end": q_offset,
            "target_start": t_offset,
            "target_end": t_offset,
            "length": 0,
            "avg_cosine": 0.0,
            "gap_opens": gap_opens,
            "gap_bases": gap_bases,
        }
        return summary, trace

    summary = {
        "query_start": int(min(query_positions)),
        "query_end": int(max(query_positions) + 1),
        "target_start": int(min(target_positions)),
        "target_end": int(max(target_positions) + 1),
        "length": len(path),
        "avg_cosine": float(np.mean(cosines) if cosines else 0.0),
        "gap_opens": gap_opens,
        "gap_bases": gap_bases,
    }
    return summary, trace


def build_alignment_strings(
    path: List[Tuple[str, int, int]],
    query_sequence: Optional[str],
    target_sequence: Optional[str],
    q_offset: int,
    t_offset: int,
) -> Tuple[str, str, str]:
    if not path:
        return "", "", ""

    q_chars: List[str] = []
    t_chars: List[str] = []
    mask_chars: List[str] = []

    q_idx = q_offset
    t_idx = t_offset

    q_len = len(query_sequence) if query_sequence else 0
    t_len = len(target_sequence) if target_sequence else 0

    def safe_get(seq: Optional[str], idx: int, length: int) -> str:
        if seq is None or idx < 0 or idx >= length:
            return "N"
        return seq[idx]

    for state, _, _ in path:
        if state == "M":
            q_chars.append(safe_get(query_sequence, q_idx, q_len))
            t_chars.append(safe_get(target_sequence, t_idx, t_len))
            mask_chars.append("|")
            q_idx += 1
            t_idx += 1
        elif state == "X":  # gap in query
            q_chars.append("-")
            t_chars.append(safe_get(target_sequence, t_idx, t_len))
            mask_chars.append(" ")
            t_idx += 1
        elif state == "Y":  # gap in target
            q_chars.append(safe_get(query_sequence, q_idx, q_len))
            t_chars.append("-")
            mask_chars.append(" ")
            q_idx += 1

    return "".join(q_chars), "".join(t_chars), "".join(mask_chars)


def write_alignment_text(output_path: Path, df: pd.DataFrame) -> None:
    lines: List[str] = []
    for idx, row in enumerate(df.to_dict("records"), start=1):
        content = row.get("alignment_content")
        if content == "structure":
            q_aln = row.get("aligned_query_structure")
            t_aln = row.get("aligned_target_structure")
        elif content == "sequence":
            q_aln = row.get("aligned_query_sequence")
            t_aln = row.get("aligned_target_sequence")
        else:
            q_aln = row.get("aligned_query_structure") or row.get("aligned_query_sequence")
            t_aln = row.get("aligned_target_structure") or row.get("aligned_target_sequence")
            if row.get("aligned_query_structure") and row.get("aligned_target_structure"):
                content = "structure"
            elif row.get("aligned_query_sequence") and row.get("aligned_target_sequence"):
                content = "sequence"
        mask = row.get("alignment_mask")
        if not q_aln or not t_aln:
            continue
        lines.append(f"Alignment {idx}: {row['query_id']} vs {row['target_id']} (cluster {row['cluster_id']})")
        lines.append(
            f"  Score: {row['score']:.3f}\tAvgCosine: {row['avg_cosine']:.3f}\tLen: {row['alignment_length']}\tGaps: {row['gap_opens']}/{row['gap_bases']}"
        )
        if content:
            lines.append(f"  Content: {content.capitalize()}")
        lines.append(
            f"  Query  {row['query_start']:>6}  {q_aln}  {row['query_end']:>6}"
        )
        if mask:
            lines.append(f"         {'':>6}  {mask}")
        lines.append(
            f"  Target {row['target_start']:>6}  {t_aln}  {row['target_end']:>6}"
        )
        lines.append("")
    output_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


if __name__ == "__main__":
    main()
