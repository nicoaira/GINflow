#!/usr/bin/env python3
"""Run banded Smith–Waterman alignments over clustered seeds."""
from __future__ import annotations

import argparse
import heapq
import json
from collections import defaultdict
from dataclasses import dataclass
from itertools import count
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

NEG_INF = -1e9

TRACE_STOP = 0
TRACE_FROM_M = 1
TRACE_FROM_X = 2
TRACE_FROM_Y = 3
TRACE_GAP_FROM_M = 0
TRACE_GAP_FROM_SELF = 1

EMBEDDING_CHUNK_SIZE = 200_000


def safe_filename(text: str) -> str:
    """Return a filesystem-safe representation of the given text."""
    return "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in text)


def prepare_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


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
    bit_score: Optional[float] = None
    evalue: Optional[float] = None


@dataclass
class AlignmentBundle:
    score: float
    record: dict
    plot_info: Optional[dict]
    dp_record: Optional[dict]


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
    parser.add_argument(
        "--calculate-evalue",
        dest="calculate_evalue",
        action="store_true",
        default=True,
        help="Calculate E-values for alignments (default: enabled)",
    )
    parser.add_argument(
        "--no-calculate-evalue",
        dest="calculate_evalue",
        action="store_false",
        help="Disable E-value calculation",
    )
    parser.add_argument(
        "--evd-samples",
        type=int,
        default=1000,
        help="Number of random alignment samples for EVD parameter estimation (default: 1000)",
    )
    parser.add_argument(
        "--evd-lambda",
        type=float,
        default=None,
        help="Pre-computed lambda parameter for E-value calculation (if provided, skips EVD estimation)",
    )
    parser.add_argument(
        "--evd-K",
        type=float,
        default=None,
        help="Pre-computed K parameter for E-value calculation (if provided, skips EVD estimation)",
    )
    parser.add_argument("--output", required=True, help="Output TSV for alignment summaries")
    parser.add_argument("--stats-json", default=None, help="Optional JSON stats path")
    parser.add_argument("--dp-output", default=None, help="Optional JSONL file capturing DP traces for reported alignments")
    parser.add_argument("--alignment-text", default=None, help="Optional plain-text alignment dump (BLAST-style)")
    parser.add_argument(
        "--plot-scoring-matrices",
        action="store_true",
        help="Plot the banded scoring matrices for the reported alignments",
    )
    parser.add_argument("--plot-dir", default=None, help="Directory for scoring matrix plots")
    return parser.parse_args()


def load_embeddings(
    path: str,
    id_col: str,
    node_idx_col: str,
    emb_col: str,
    required_ids: Optional[Set[str]] = None,
    chunk_size: int = EMBEDDING_CHUNK_SIZE,
) -> Dict[str, np.ndarray]:
    usecols = [id_col, node_idx_col, emb_col]
    read_kwargs = {
        "sep": "\t",
        "dtype": {id_col: str, node_idx_col: np.int64, emb_col: str},
        "usecols": usecols,
    }
    try:
        iterator = pd.read_csv(path, chunksize=chunk_size, **read_kwargs)
    except ValueError as exc:
        raise SystemExit(f"Failed to read embeddings from {path}: {exc}") from exc

    collected: Dict[str, List[Tuple[int, np.ndarray]]] = defaultdict(list)
    for chunk in iterator:
        if required_ids:
            mask = chunk[id_col].isin(required_ids)
            if not mask.any():
                continue
            chunk = chunk.loc[mask]
        ids = chunk[id_col].astype(str)
        node_indices = chunk[node_idx_col].astype(np.int64)
        vectors = chunk[emb_col].astype(str)
        for tid, node_idx, vec_text in zip(ids, node_indices, vectors):
            if not vec_text:
                continue
            vec = np.fromstring(vec_text, sep=",", dtype=np.float32)
            if vec.size == 0:
                continue
            collected[tid].append((int(node_idx), vec))

    if not collected and required_ids:
        raise SystemExit("No embeddings found for the requested transcript IDs")

    embeddings: Dict[str, np.ndarray] = {}
    for tid, entries in collected.items():
        if not entries:
            continue
        entries.sort(key=lambda item: item[0])
        stacked = np.vstack([vec for _, vec in entries]).astype(np.float32, copy=False)
        embeddings[tid] = stacked

    if required_ids:
        missing = required_ids - embeddings.keys()
        if missing:
            preview = ", ".join(sorted(missing)[:10])
            raise SystemExit(f"Embeddings file missing vectors for {len(missing)} transcripts (e.g. {preview})")

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


def shuffle_embeddings(embeddings: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle the order of embedding vectors to create a random sequence.

    This preserves the composition (same vectors) but destroys sequential information.

    Args:
        embeddings: Array of shape (n_nodes, embedding_dim)
        rng: Random number generator

    Returns:
        Shuffled copy of embeddings
    """
    shuffled = embeddings.copy()
    rng.shuffle(shuffled)
    return shuffled


def estimate_evd_parameters(
    embeddings: Dict[str, np.ndarray],
    mu: float,
    sigma: float,
    gamma: float,
    band: int,
    gap_open: float,
    gap_extend: float,
    xdrop: float,
    clamp_min: float,
    clamp_max: float,
    num_samples: int = 1000,
    seed: int = 42,
) -> Tuple[float, float]:
    """Estimate lambda and K parameters for BLAST-like E-value calculation.

    Generates random sequence pairs by shuffling, aligns them to get a null
    distribution of scores, and fits an Extreme Value Distribution (Gumbel).

    Args:
        embeddings: Dictionary mapping sequence IDs to embedding arrays
        mu: Background mean for scoring
        sigma: Background std for scoring
        gamma: Scaling factor for scores
        band: Band width for alignment
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
        xdrop: X-drop threshold
        clamp_min: Minimum score clamp
        clamp_max: Maximum score clamp
        num_samples: Number of random pairs to generate
        seed: Random seed

    Returns:
        Tuple of (lambda, K) parameters
    """
    from scipy.stats import gumbel_r

    rng = np.random.default_rng(seed)

    # Get list of sequences that have sufficient length
    seq_ids = [tid for tid, vec in embeddings.items() if len(vec) >= 20]
    if len(seq_ids) < 2:
        # Not enough data, return default values
        return 0.1, 0.01

    null_scores: List[float] = []

    print(f"Estimating EVD parameters from {num_samples} random alignments...", flush=True)

    for i in range(num_samples):
        # Sample two random sequences
        id1, id2 = rng.choice(seq_ids, size=2, replace=True)

        # Get embeddings and shuffle them
        emb1 = shuffle_embeddings(embeddings[id1], rng)
        emb2 = shuffle_embeddings(embeddings[id2], rng)

        # Align with Smith-Waterman
        score, path, _ = smith_waterman(
            emb1,
            emb2,
            mu,
            sigma,
            gamma,
            band,
            diag_offset=0,  # Random offset for null model
            gap_open=gap_open,
            gap_extend=gap_extend,
            xdrop=xdrop,
            clamp_min=clamp_min,
            clamp_max=clamp_max,
            store_matrix=False,
        )

        if score > 0:
            null_scores.append(float(score))

        if (i + 1) % 100 == 0:
            print(f"  Completed {i + 1}/{num_samples} random alignments...", flush=True)

    if len(null_scores) < 10:
        print("Warning: Very few random alignments produced positive scores. Using default parameters.", flush=True)
        return 0.1, 0.01

    # Fit Gumbel distribution (right-skewed EVD)
    # The Gumbel distribution parameters are (loc=mu, scale=beta)
    fit_mu, fit_beta = gumbel_r.fit(null_scores)

    # Convert to Karlin-Altschul parameters
    lambda_param = 1.0 / fit_beta
    K = np.exp(-lambda_param * fit_mu)

    print(f"EVD fitting complete:", flush=True)
    print(f"  Gumbel location (μ) = {fit_mu:.4f}", flush=True)
    print(f"  Gumbel scale (β) = {fit_beta:.4f}", flush=True)
    print(f"  Lambda (λ) = {lambda_param:.4f}", flush=True)
    print(f"  K = {K:.6f}", flush=True)
    print(f"  Mean null score = {np.mean(null_scores):.4f}", flush=True)
    print(f"  Std null score = {np.std(null_scores):.4f}", flush=True)

    return float(lambda_param), float(K)


def calculate_bit_score(raw_score: float, lambda_param: float, K: float) -> float:
    """Convert raw alignment score to bit score.

    Args:
        raw_score: Raw alignment score from dynamic programming
        lambda_param: Lambda parameter from EVD fitting
        K: K parameter from EVD fitting

    Returns:
        Bit score (normalized, statistically meaningful score)
    """
    if lambda_param <= 0 or K <= 0:
        return 0.0

    bit_score = (lambda_param * raw_score - np.log(K)) / np.log(2)
    return float(bit_score)


def calculate_evalue(
    raw_score: float,
    query_length: int,
    target_length: int,
    lambda_param: float,
    K: float,
) -> Tuple[float, float]:
    """Calculate BLAST-like E-value and bit score for an alignment.

    Uses the Karlin-Altschul statistics based on Extreme Value Distribution.

    The E-value represents the expected number of alignments with score >= S
    that would occur by chance in a search space of this size.

    Args:
        raw_score: Raw alignment score from dynamic programming
        query_length: Length of the query sequence (in nodes)
        target_length: Length of the target sequence (in nodes)
        lambda_param: Lambda parameter from EVD fitting
        K: K parameter from EVD fitting

    Returns:
        Tuple of (E-value, bit_score)
    """
    if query_length == 0 or target_length == 0:
        return float('inf'), 0.0

    if lambda_param <= 0 or K <= 0:
        return float('inf'), 0.0

    # Calculate bit score
    bit_score = calculate_bit_score(raw_score, lambda_param, K)

    # Calculate E-value: E = m × n × 2^(-S')
    # where m and n are sequence lengths and S' is the bit score
    search_space = float(query_length * target_length)
    evalue = search_space * np.power(2.0, -bit_score)

    return float(evalue), float(bit_score)


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
    store_matrix: bool = False,
) -> Tuple[float, List[Tuple[str, int, int]], Optional[np.ndarray]]:
    m, n = len(query), len(target)
    if m == 0 or n == 0:
        return 0.0, [], None

    band = max(1, int(band))
    band_half = band // 2

    M_prev = np.zeros(n + 1, dtype=np.float32)
    X_prev = np.full(n + 1, NEG_INF, dtype=np.float32)
    Y_prev = np.full(n + 1, NEG_INF, dtype=np.float32)
    M_curr = np.zeros(n + 1, dtype=np.float32)
    X_curr = np.full(n + 1, NEG_INF, dtype=np.float32)
    Y_curr = np.full(n + 1, NEG_INF, dtype=np.float32)

    band_capacity = min(n, band + 1) + 2
    trace_M = np.zeros((m + 1, band_capacity), dtype=np.uint8)
    trace_X = np.zeros((m + 1, band_capacity), dtype=np.uint8)
    trace_Y = np.zeros((m + 1, band_capacity), dtype=np.uint8)
    row_offsets = np.zeros(m + 1, dtype=np.int32)
    row_lengths = np.zeros(m + 1, dtype=np.int32)

    matrix_scores: Optional[np.ndarray] = None
    if store_matrix:
        matrix_scores = np.full((m + 1, n + 1), np.nan, dtype=np.float32)
        matrix_scores[0, :] = 0.0
        matrix_scores[:, 0] = 0.0

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
            row_offsets[i] = 1
            row_lengths[i] = 0
            continue

        row_offsets[i] = j_start
        span = j_end - j_start + 1
        row_lengths[i] = span

        for j in range(j_start, j_end + 1):
            i_global = i - 1
            j_global = j - 1
            local_idx = j - j_start
            if local_idx < 0 or local_idx >= band_capacity:
                continue
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
            if base_ptr is None:
                trace_M[i, local_idx] = TRACE_STOP
            else:
                prev_state = base_ptr[0]
                if prev_state == "M":
                    trace_M[i, local_idx] = TRACE_FROM_M
                elif prev_state == "X":
                    trace_M[i, local_idx] = TRACE_FROM_X
                else:
                    trace_M[i, local_idx] = TRACE_FROM_Y

            # Gap in query (insertion in target)
            left_m = M_curr[j - 1] - gap_open
            left_x = X_curr[j - 1] - gap_extend
            if left_m >= left_x:
                X_curr[j] = left_m
                trace_X[i, local_idx] = TRACE_GAP_FROM_M
            else:
                X_curr[j] = left_x
                trace_X[i, local_idx] = TRACE_GAP_FROM_SELF

            # Gap in target (deletion)
            up_m = M_prev[j] - gap_open
            up_y = Y_prev[j] - gap_extend
            if up_m >= up_y:
                Y_curr[j] = up_m
                trace_Y[i, local_idx] = TRACE_GAP_FROM_M
            else:
                Y_curr[j] = up_y
                trace_Y[i, local_idx] = TRACE_GAP_FROM_SELF

            cell_best = max(match_val, X_curr[j], Y_curr[j])
            if matrix_scores is not None:
                matrix_scores[i, j] = cell_best
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

    path = traceback(best_state, trace_M, trace_X, trace_Y, row_offsets, row_lengths)
    return best_score, path, matrix_scores


def traceback(
    best_state: Tuple[str, int, int],
    trace_M: np.ndarray,
    trace_X: np.ndarray,
    trace_Y: np.ndarray,
    row_offsets: np.ndarray,
    row_lengths: np.ndarray,
) -> List[Tuple[str, int, int]]:
    state, i, j = best_state
    if i <= 0 or j <= 0:
        return []
    path: List[Tuple[str, int, int]] = []
    while state in {"M", "X", "Y"} and i > 0 and j > 0:
        local_offset = row_offsets[i]
        local_len = row_lengths[i]
        local_idx = j - local_offset
        if local_idx < 0 or local_idx >= local_len or local_idx >= trace_M.shape[1]:
            break
        path.append((state, i, j))
        if state == "M":
            code = trace_M[i, local_idx]
            if code == TRACE_STOP:
                break
            if code == TRACE_FROM_M:
                state = "M"
            elif code == TRACE_FROM_X:
                state = "X"
            elif code == TRACE_FROM_Y:
                state = "Y"
            else:
                break
            i -= 1
            j -= 1
        elif state == "X":
            code = trace_X[i, local_idx]
            if j <= 0:
                break
            j -= 1
            if code == TRACE_GAP_FROM_M:
                state = "M"
            elif code == TRACE_GAP_FROM_SELF:
                state = "X"
            else:
                break
        else:  # Y
            code = trace_Y[i, local_idx]
            if i <= 0:
                break
            i -= 1
            if code == TRACE_GAP_FROM_M:
                state = "M"
            elif code == TRACE_GAP_FROM_SELF:
                state = "Y"
            else:
                break
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
        if args.dp_output:
            Path(args.dp_output).write_text("")
        if args.alignment_text:
            Path(args.alignment_text).write_text("")
        if args.stats_json:
            Path(args.stats_json).write_text(json.dumps({"alignments": 0}, indent=2) + "\n")
        return

    required_cluster_cols = {"cluster_id", "seed_count", "max_similarity"}
    missing_cluster_cols = required_cluster_cols - set(clusters_df.columns)
    if missing_cluster_cols:
        raise SystemExit(
            f"Cluster summaries file missing columns: {', '.join(sorted(missing_cluster_cols))}"
        )

    member_cols_order = [
        "cluster_id",
        "query_transcript",
        "target_transcript",
        "diagonal",
        "query_window_start",
        "query_window_end",
        "target_window_start",
        "target_window_end",
    ]
    missing_member_cols = set(member_cols_order) - set(members_df.columns)
    if missing_member_cols:
        raise SystemExit(
            f"Cluster members file missing columns: {', '.join(sorted(missing_member_cols))}"
        )
    members_df = members_df.loc[:, member_cols_order]

    members_df["query_transcript"] = members_df["query_transcript"].astype(str)
    members_df["target_transcript"] = members_df["target_transcript"].astype(str)

    query_ids = {tid for tid in members_df["query_transcript"] if tid and tid.lower() != "nan"}
    target_ids = {tid for tid in members_df["target_transcript"] if tid and tid.lower() != "nan"}
    required_ids = query_ids | target_ids

    embeddings = load_embeddings(
        args.embeddings,
        args.id_column,
        args.node_index_column,
        args.embedding_column,
        required_ids=required_ids if required_ids else None,
    )
    mu0, sigma0 = sample_background(embeddings, args.background_samples, args.random_seed)

    # Estimate EVD parameters (lambda and K) for E-value calculation
    lambda_param: Optional[float] = None
    K_param: Optional[float] = None
    if args.calculate_evalue:
        if args.evd_lambda is not None and args.evd_K is not None:
            # Use pre-computed parameters
            lambda_param = args.evd_lambda
            K_param = args.evd_K
            print(f"Using pre-computed EVD parameters: λ={lambda_param:.4f}, K={K_param:.6f}", flush=True)
        else:
            # Estimate parameters from random alignments
            lambda_param, K_param = estimate_evd_parameters(
                embeddings,
                mu0,
                sigma0,
                args.gamma,
                args.band_width,
                args.gap_open,
                args.gap_extend,
                args.xdrop,
                args.score_min,
                args.score_max,
                num_samples=args.evd_samples,
                seed=args.random_seed,
            )

    meta_df = pd.read_csv(args.meta, sep="\t", dtype={args.id_column: str})
    meta_df = meta_df.drop_duplicates(subset=[args.id_column])
    meta_map = meta_df.set_index(args.id_column).to_dict(orient="index")
    seq_col, struct_col = infer_columns(meta_df, args.sequence_column, args.structure_column)

    # Store sequence lengths for E-value computation
    sequence_lengths = {tid: len(vec) for tid, vec in embeddings.items()}

    plot_enabled = args.plot_scoring_matrices
    collect_trace = args.dp_output is not None
    top_n = args.top_n
    store_results = top_n != 0
    use_heap = store_results and top_n > 0

    result_heap: List[Tuple[float, int, AlignmentBundle]] = []
    unlimited_results: List[AlignmentBundle] = []
    heap_counter = count()

    grouped_members = members_df.groupby("cluster_id", sort=False)

    for cluster in clusters_df.itertuples(index=False):
        cid = getattr(cluster, "cluster_id")
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

        score, path, _ = smith_waterman(
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

        # Calculate E-value and bit score if enabled
        evalue = None
        bit_score = None
        if args.calculate_evalue and lambda_param is not None and K_param is not None:
            query_length = sequence_lengths.get(query_id, 0)
            target_length = sequence_lengths.get(target_id, 0)
            if query_length > 0 and target_length > 0:
                evalue, bit_score = calculate_evalue(
                    raw_score=float(score),
                    query_length=query_length,
                    target_length=target_length,
                    lambda_param=lambda_param,
                    K=K_param,
                )

        result = AlignmentResult(
            query_id=query_id,
            target_id=target_id,
            cluster_id=int(cid),
            score=float(score),
            query_start=alignment_info["query_start"],
            query_end=alignment_info["query_end"],
            target_start=alignment_info["target_start"],
            target_end=alignment_info["target_end"],
            alignment_length=alignment_info["length"],
            avg_cosine=alignment_info["avg_cosine"],
            gap_opens=alignment_info["gap_opens"],
            gap_bases=alignment_info["gap_bases"],
            seed_count=int(getattr(cluster, "seed_count")),
            max_seed_similarity=float(getattr(cluster, "max_similarity")),
            bit_score=bit_score,
            evalue=evalue,
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
        plot_info = None
        if plot_enabled:
            plot_info = {
                "query_start": q_start,
                "query_end": q_end,
                "target_start": t_start,
                "target_end": t_end,
                "query_len": len(query_vecs),
                "target_len": len(target_vecs),
                "band_width": effective_band,
                "diag_offset": diag_offset,
                "diag_median": diag_median,
                "diag_span": diag_span,
                "diag_min": diag_min,
                "diag_max": diag_max,
            }

        dp_record = None
        if collect_trace and trace is not None:
            dp_record = {
                "query_id": query_id,
                "target_id": target_id,
                "cluster_id": int(cid),
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

        bundle = AlignmentBundle(
            score=float(score),
            record=enriched,
            plot_info=plot_info,
            dp_record=dp_record,
        )

        if not store_results:
            continue

        if use_heap:
            heapq.heappush(result_heap, (bundle.score, next(heap_counter), bundle))
            if len(result_heap) > top_n:
                heapq.heappop(result_heap)
        else:
            unlimited_results.append(bundle)

    if not store_results:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        if args.dp_output:
            Path(args.dp_output).write_text("")
        if args.alignment_text:
            Path(args.alignment_text).write_text("")
        if args.stats_json:
            Path(args.stats_json).write_text(json.dumps({"alignments": 0}, indent=2) + "\n")
        return

    if use_heap:
        selected_bundles = [entry[2] for entry in sorted(result_heap, key=lambda item: item[0], reverse=True)]
    else:
        selected_bundles = sorted(unlimited_results, key=lambda item: item.score, reverse=True)
        if top_n < 0:
            drop = abs(top_n)
            if drop >= len(selected_bundles):
                selected_bundles = []
            else:
                selected_bundles = selected_bundles[:-drop]

    if not selected_bundles:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        if args.dp_output:
            Path(args.dp_output).write_text("")
        if args.alignment_text:
            Path(args.alignment_text).write_text("")
        if args.stats_json:
            Path(args.stats_json).write_text(json.dumps({"alignments": 0}, indent=2) + "\n")
        return

    records = [bundle.record for bundle in selected_bundles]
    results_df = pd.DataFrame(records)
    results_df = results_df.sort_values("score", ascending=False)
    results_df.to_csv(args.output, sep="\t", index=False)

    plot_count = 0
    if plot_enabled and not results_df.empty:
        plot_dir = Path(args.plot_dir) if args.plot_dir else Path("alignment_plots")
        plot_dir.mkdir(parents=True, exist_ok=True)
        try:
            plt = prepare_matplotlib()
        except ImportError as exc:  # pragma: no cover - dependency should be present in env
            raise SystemExit(
                "matplotlib is required for --plot-scoring-matrices; ensure the environment includes matplotlib"
            ) from exc
        for rank, bundle in enumerate(selected_bundles, start=1):
            info = bundle.plot_info
            if not info:
                continue
            row = bundle.record
            query_id = str(row["query_id"])
            target_id = str(row["target_id"])
            q_slice = embeddings[query_id][info["query_start"] : info["query_end"]]
            t_slice = embeddings[target_id][info["target_start"] : info["target_end"]]
            _, path_plot, matrix = smith_waterman(
                q_slice,
                t_slice,
                mu0,
                sigma0,
                args.gamma,
                info["band_width"],
                info["diag_offset"],
                args.gap_open,
                args.gap_extend,
                args.xdrop,
                args.score_min,
                args.score_max,
                store_matrix=True,
            )
            if matrix is None or matrix.size == 0:
                continue
            metadata = {
                "rank": rank,
                "score": float(row["score"]),
                "query_id": query_id,
                "target_id": target_id,
                "cluster_id": int(row["cluster_id"]),
            }
            base_name = (
                f"rank_{rank:02d}_{safe_filename(query_id)}_vs_{safe_filename(target_id)}"
                f"_cluster_{int(row['cluster_id'])}_score_{metadata['score']:.2f}"
            )
            base_name = safe_filename(base_name)
            output_path = plot_dir / f"{base_name}.png"
            plot_scoring_matrix(matrix, path_plot, info, metadata, output_path, plt)
            plot_count += 1

    stats = {"alignments": int(len(results_df)), "mu0": mu0, "sigma0": sigma0}
    if lambda_param is not None and K_param is not None:
        stats["evd_lambda"] = float(lambda_param)
        stats["evd_K"] = float(K_param)
    if plot_enabled:
        stats["plots"] = int(plot_count)
    if args.stats_json:
        Path(args.stats_json).write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))

    if args.dp_output:
        wrote_trace = False
        with Path(args.dp_output).open("w", encoding="utf-8") as handle:
            for bundle in selected_bundles:
                if bundle.dp_record:
                    handle.write(json.dumps(bundle.dp_record) + "\n")
                    wrote_trace = True
        if not wrote_trace:
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


def plot_scoring_matrix(
    matrix: np.ndarray,
    path: List[Tuple[str, int, int]],
    info: Dict[str, int],
    metadata: Dict[str, object],
    output_path: Path,
    plt,
) -> None:
    if matrix.shape[0] <= 1 or matrix.shape[1] <= 1:
        return
    band_matrix = matrix[1:, 1:]
    if band_matrix.size == 0:
        return

    query_len = int(info.get("query_len", band_matrix.shape[0]))
    target_len = int(info.get("target_len", band_matrix.shape[1]))
    if query_len <= 0 or target_len <= 0:
        return

    full_matrix = np.full((query_len, target_len), np.nan, dtype=np.float32)
    query_start = max(0, min(int(info["query_start"]), query_len))
    query_end = max(query_start, min(int(info["query_end"]), query_len))
    target_start = max(0, min(int(info["target_start"]), target_len))
    target_end = max(target_start, min(int(info["target_end"]), target_len))

    slice_h = min(band_matrix.shape[0], query_end - query_start)
    slice_w = min(band_matrix.shape[1], target_end - target_start)
    if slice_h > 0 and slice_w > 0:
        full_matrix[query_start : query_start + slice_h, target_start : target_start + slice_w] = band_matrix[:slice_h, :slice_w]

    masked = np.ma.masked_invalid(full_matrix)
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(masked, origin="lower", aspect="auto", cmap="viridis")
    colorbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("SW score", fontsize=9)

    if path and slice_h > 0 and slice_w > 0:
        y_coords = [query_start + i - 1 for _, i, _ in path if i > 0]
        x_coords = [target_start + j - 1 for _, _, j in path if j > 0]
        if x_coords and y_coords and len(x_coords) == len(y_coords):
            ax.plot(x_coords, y_coords, color="red", linewidth=1.0, alpha=0.85)

    x_ticks = [0] + ([target_len - 1] if target_len > 1 else [])
    y_ticks = [0] + ([query_len - 1] if query_len > 1 else [])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks)

    if slice_h > 0 and slice_w > 0:
        from matplotlib.patches import Rectangle

        rect = Rectangle(
            (target_start - 0.5, query_start - 0.5),
            slice_w,
            slice_h,
            linewidth=1.0,
            edgecolor="white",
            facecolor="none",
            linestyle="--",
            alpha=0.7,
        )
        ax.add_patch(rect)

    ax.set_xlabel("Target window index")
    ax.set_ylabel("Query window index")
    ax.set_title(
        f"Rank {metadata['rank']}: {metadata['query_id']} vs {metadata['target_id']} (cluster {metadata['cluster_id']})"
    )

    diag_span = info.get("diag_span")
    diag_median = info.get("diag_median")
    subtitle = f"Score {metadata['score']:.3f}"
    band_width = info.get("band_width")
    subtitle += f" | Q[{query_start},{query_end}) vs T[{target_start},{target_end})"
    if diag_median is not None and diag_span is not None:
        subtitle += f" | Diagonal μ {diag_median}, span {diag_span}"
    if band_width is not None:
        subtitle += f" | Band {band_width}"
    ax.text(
        0.02,
        0.98,
        subtitle,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.6, "edgecolor": "none"},
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


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
