#!/usr/bin/env python3
"""Calculate EVD parameters (lambda and K) for E-value estimation.

This script calculates database-wide Extreme Value Distribution parameters
that can be reused across multiple query alignments. The parameters are
estimated by performing random alignments on shuffled sequences and fitting
a Gumbel distribution to the resulting score distribution.
"""
from __future__ import annotations

import argparse
import json
import sys
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate EVD parameters for E-value estimation"
    )
    parser.add_argument(
        "--embeddings", required=True, help="Node-level embeddings TSV"
    )
    parser.add_argument(
        "--id-column",
        default="transcript_id",
        help="Identifier column in embeddings table",
    )
    parser.add_argument(
        "--node-index-column",
        default="node_index",
        help="Column giving node order within each transcript",
    )
    parser.add_argument(
        "--embedding-column",
        default="embedding_vector",
        help="Column containing comma-separated embedding vectors",
    )
    parser.add_argument(
        "--background-samples",
        type=int,
        default=10000,
        help="Number of random pairs to sample for background statistics",
    )
    parser.add_argument(
        "--evd-samples",
        type=int,
        default=1000,
        help="Number of random alignment samples for EVD parameter estimation",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=1.5,
        help="Scaling factor for similarity scores",
    )
    parser.add_argument(
        "--band-width", type=int, default=96, help="Band width for banded alignment"
    )
    parser.add_argument(
        "--gap-open", type=float, default=12, help="Gap opening penalty"
    )
    parser.add_argument(
        "--gap-extend", type=float, default=2, help="Gap extension penalty"
    )
    parser.add_argument(
        "--xdrop", type=float, default=50, help="X-drop threshold for alignment"
    )
    parser.add_argument(
        "--score-min", type=float, default=-4, help="Minimum score clamp"
    )
    parser.add_argument(
        "--score-max", type=float, default=8, help="Maximum score clamp"
    )
    parser.add_argument(
        "--random-seed", type=int, default=42, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--output", required=True, help="Output JSON file with EVD parameters"
    )
    return parser.parse_args()


def load_embeddings(
    path: str,
    id_col: str,
    node_idx_col: str,
    emb_col: str,
    required_ids: Optional[Set[str]] = None,
) -> Dict[str, np.ndarray]:
    """Load embeddings from TSV into a dictionary keyed by transcript ID."""
    embeddings: Dict[str, List[np.ndarray]] = {}
    chunk_iter = pd.read_csv(
        path, sep="\t", chunksize=EMBEDDING_CHUNK_SIZE, dtype={id_col: str}
    )

    for chunk_df in chunk_iter:
        chunk_df = chunk_df.sort_values([id_col, node_idx_col])
        for tid, group in chunk_df.groupby(id_col, sort=False):
            tid_str = str(tid)
            if required_ids is not None and tid_str not in required_ids:
                continue
            vectors = [
                np.fromstring(row[emb_col], sep=",", dtype=np.float32)
                for _, row in group.iterrows()
            ]
            if tid_str not in embeddings:
                embeddings[tid_str] = []
            embeddings[tid_str].extend(vectors)

    # Convert lists to numpy arrays
    final_embeddings = {tid: np.array(vec_list) for tid, vec_list in embeddings.items()}
    return final_embeddings


def sample_background(
    embeddings: Dict[str, np.ndarray], sample_count: int, seed: int
) -> Tuple[float, float]:
    """Sample random node pairs to estimate background similarity distribution."""
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


def compute_score(
    dot: float, mu: float, sigma: float, gamma: float, smin: float, smax: float
) -> float:
    """Convert dot product to scaled score."""
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
    """Banded Smith-Waterman alignment with affine gaps.

    This is a simplified version for EVD parameter estimation that only
    needs to return the alignment score.

    Args:
        query: Query embedding array (m x d)
        target: Target embedding array (n x d)
        mu: Background mean for scoring
        sigma: Background std for scoring
        gamma: Scaling factor for scores
        band: Band width
        diag_offset: Diagonal offset for band center
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
        xdrop: X-drop threshold
        clamp_min: Minimum score clamp
        clamp_max: Maximum score clamp
        store_matrix: Whether to store full DP matrix (unused for EVD)

    Returns:
        Tuple of (score, path, matrix)
    """
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

            # Gap in target (insertion in query)
            up_m = M_prev[j] - gap_open
            up_y = Y_prev[j] - gap_extend
            if up_m >= up_y:
                Y_curr[j] = up_m
                trace_Y[i, local_idx] = TRACE_GAP_FROM_M
            else:
                Y_curr[j] = up_y
                trace_Y[i, local_idx] = TRACE_GAP_FROM_SELF

            # Track best score
            curr_best = max(M_curr[j], X_curr[j], Y_curr[j])
            if curr_best > best_score:
                best_score = curr_best
                if M_curr[j] >= X_curr[j] and M_curr[j] >= Y_curr[j]:
                    best_state = ("M", i, j)
                elif X_curr[j] >= Y_curr[j]:
                    best_state = ("X", i, j)
                else:
                    best_state = ("Y", i, j)

            if curr_best > row_best:
                row_best = curr_best

            if store_matrix:
                matrix_scores[i, j] = curr_best

        # X-drop pruning
        if best_score - row_best > xdrop and best_score > 0:
            # Continue for a few more rows to ensure we don't prematurely stop
            pass

        M_prev, M_curr = M_curr, M_prev
        X_prev, X_curr = X_curr, X_prev
        Y_prev, Y_curr = Y_curr, Y_prev

    # For EVD estimation, we only need the score, not the traceback
    return float(best_score), [], matrix_scores


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
        print(
            "Warning: Not enough sequences for EVD estimation. Using default parameters.",
            file=sys.stderr,
            flush=True,
        )
        return 0.1, 0.01

    null_scores: List[float] = []

    print(
        f"Estimating EVD parameters from {num_samples} random alignments...",
        file=sys.stderr,
        flush=True,
    )

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
            print(
                f"  Completed {i + 1}/{num_samples} random alignments...",
                file=sys.stderr,
                flush=True,
            )

    if len(null_scores) < 10:
        print(
            "Warning: Very few random alignments produced positive scores. Using default parameters.",
            file=sys.stderr,
            flush=True,
        )
        return 0.1, 0.01

    # Fit Gumbel distribution (right-skewed EVD)
    # The Gumbel distribution parameters are (loc=mu, scale=beta)
    fit_mu, fit_beta = gumbel_r.fit(null_scores)

    # Convert to Karlin-Altschul parameters
    lambda_param = 1.0 / fit_beta
    K = np.exp(-lambda_param * fit_mu)

    print(f"EVD fitting complete:", file=sys.stderr, flush=True)
    print(f"  Gumbel location (μ) = {fit_mu:.4f}", file=sys.stderr, flush=True)
    print(f"  Gumbel scale (β) = {fit_beta:.4f}", file=sys.stderr, flush=True)
    print(f"  Lambda (λ) = {lambda_param:.4f}", file=sys.stderr, flush=True)
    print(f"  K = {K:.6f}", file=sys.stderr, flush=True)
    print(f"  Mean null score = {np.mean(null_scores):.4f}", file=sys.stderr, flush=True)
    print(f"  Std null score = {np.std(null_scores):.4f}", file=sys.stderr, flush=True)

    return float(lambda_param), float(K)


def main() -> None:
    args = parse_args()

    print("Loading embeddings...", file=sys.stderr, flush=True)
    embeddings = load_embeddings(
        args.embeddings,
        args.id_column,
        args.node_index_column,
        args.embedding_column,
        required_ids=None,
    )

    print(
        f"Loaded {len(embeddings)} sequences with embeddings", file=sys.stderr, flush=True
    )

    print("Sampling background statistics...", file=sys.stderr, flush=True)
    mu0, sigma0 = sample_background(embeddings, args.background_samples, args.random_seed)
    print(f"  Background μ = {mu0:.4f}", file=sys.stderr, flush=True)
    print(f"  Background σ = {sigma0:.4f}", file=sys.stderr, flush=True)

    # Estimate EVD parameters
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

    # Write output JSON
    output_data = {
        "evd_lambda": float(lambda_param),
        "evd_K": float(K_param),
        "background_mu": float(mu0),
        "background_sigma": float(sigma0),
        "gamma": args.gamma,
        "band_width": args.band_width,
        "gap_open": args.gap_open,
        "gap_extend": args.gap_extend,
        "xdrop": args.xdrop,
        "score_min": args.score_min,
        "score_max": args.score_max,
        "evd_samples": args.evd_samples,
        "background_samples": args.background_samples,
        "random_seed": args.random_seed,
        "num_sequences": len(embeddings),
    }

    Path(args.output).write_text(json.dumps(output_data, indent=2) + "\n")

    print(f"\nEVD parameters written to {args.output}", file=sys.stderr, flush=True)
    print(json.dumps(output_data, indent=2))


if __name__ == "__main__":
    main()
