#!/usr/bin/env python3
"""Estimate Karlin–Altschul (λ, K) from a reverse-sequence SW null."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from ginfinity_sw import ScoringParameters, align_multiple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from record_pack import load_residue_embeddings  # noqa: E402


EULER_MASCHERONI = 0.5772156649015329
LN2 = math.log(2.0)


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def load_parameters(path: Path) -> tuple[ScoringParameters, str]:
    raw = path.read_bytes()
    payload = json.loads(raw.decode())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    params = ScoringParameters(**payload.get("scoring_parameters", payload))
    digest = hashlib.sha256(raw).hexdigest()
    return params, digest


def load_embeddings(path: Path) -> dict[str, np.ndarray]:
    arrays = {
        identifier: value
        for identifier, value in load_residue_embeddings(path).items()
        if value.ndim == 2 and value.shape[0] >= 2
    }
    if len(arrays) < 2:
        raise ValueError(f"{path} needs at least two residue embeddings of length >= 2")
    return arrays


def scoring_parameters_payload(params: ScoringParameters) -> dict:
    return {
        "mu": params.mu,
        "sigma": params.sigma,
        "gamma": params.gamma,
        "score_min": params.score_min,
        "score_max": params.score_max,
        "gap_open": params.gap_open,
        "gap_extend": params.gap_extend,
        "score_offset": params.score_offset,
    }


def random_crop(matrix: np.ndarray, max_length: int, rng: np.random.Generator) -> np.ndarray:
    length = int(matrix.shape[0])
    if max_length <= 0 or length <= max_length:
        return matrix
    start = int(rng.integers(0, length - max_length + 1))
    return matrix[start : start + max_length]


def sample_null_scores(
    embeddings: dict[str, np.ndarray],
    params: ScoringParameters,
    n_samples: int,
    max_length: int,
    max_cells: int,
    max_alignments: int,
    min_score: float,
    min_match_count: int,
    rng: np.random.Generator,
    workers: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    identifiers = list(embeddings)
    n_ids = len(identifiers)
    scores = []
    query_lengths = []
    target_lengths = []
    # Warm the Numba kernel on a tiny pair so the sample loop is not dominated
    # by compilation.
    probe_q = embeddings[identifiers[0]][: min(8, embeddings[identifiers[0]].shape[0])]
    probe_t = embeddings[identifiers[1]][: min(8, embeddings[identifiers[1]].shape[0])]
    align_multiple(
        probe_q,
        probe_t,
        params=params,
        max_alignments=max_alignments,
        min_score=min_score,
        min_match_count=min_match_count,
        max_cells=max_cells,
    )

    attempts = 0
    limit = max(n_samples * 8, n_samples + 16)
    pairs: list[tuple[np.ndarray, np.ndarray]] = []
    while len(pairs) < n_samples and attempts < limit:
        attempts += 1
        i, j = rng.choice(n_ids, size=2, replace=False)
        query = random_crop(embeddings[identifiers[i]][::-1], max_length, rng)
        target = random_crop(embeddings[identifiers[j]], max_length, rng)
        cells = int(query.shape[0]) * int(target.shape[0])
        if cells == 0 or cells > max_cells:
            continue
        pairs.append((query, target))

    def score_pair(pair: tuple[np.ndarray, np.ndarray]) -> tuple[float, int, int]:
        query, target = pair
        score = float(align_multiple(
            query,
            target,
            params=params,
            max_alignments=max_alignments,
            min_score=min_score,
            min_match_count=min_match_count,
            max_cells=max_cells,
        ).total_score)
        return score, int(query.shape[0]), int(target.shape[0])

    workers = max(1, int(workers))
    if workers == 1 or len(pairs) < 2:
        scored = [score_pair(pair) for pair in pairs]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            scored = list(pool.map(score_pair, pairs, chunksize=4))
    for score, query_length, target_length in scored:
        scores.append(score)
        query_lengths.append(query_length)
        target_lengths.append(target_length)
    if len(scores) < max(20, n_samples // 5):
        raise ValueError(f"only collected {len(scores)} null alignments; need more sequences or a larger --samples")
    return (
        np.asarray(scores, dtype=np.float64),
        np.asarray(query_lengths, dtype=np.int32),
        np.asarray(target_lengths, dtype=np.int32),
    )


def fit_gumbel_moments(values: np.ndarray) -> tuple[float, float]:
    mean = float(values.mean())
    var = float(values.var(ddof=1)) if values.size > 1 else 1.0
    scale = max(math.sqrt(6.0 * max(var, 1e-18)) / math.pi, 1e-8)
    loc = mean - EULER_MASCHERONI * scale
    return loc, scale


def fit_gumbel_mle(values: np.ndarray, max_iter: int = 80) -> tuple[float, float]:
    """Two-parameter Gumbel MLE. Scale is solved first, then location."""
    values = np.asarray(values, dtype=np.float64)
    loc, scale = fit_gumbel_moments(values)
    mean = float(values.mean())

    def scale_residual(candidate: float) -> float:
        weights = np.exp(-values / candidate)
        weighted = float(np.sum(values * weights) / np.sum(weights))
        return candidate - (mean - weighted)

    low, high = scale / 20.0, scale * 20.0
    low = max(low, 1e-8)
    f_low, f_high = scale_residual(low), scale_residual(high)
    if f_low * f_high > 0:
        return loc, scale
    for _ in range(max_iter):
        mid = 0.5 * (low + high)
        f_mid = scale_residual(mid)
        if abs(f_mid) < 1e-10 or high - low < 1e-10:
            scale = mid
            break
        if f_low * f_mid <= 0:
            high, f_high = mid, f_mid
        else:
            low, f_low = mid, f_mid
        scale = mid
    loc = -scale * math.log(float(np.mean(np.exp(-values / scale))))
    return loc, scale


def gumbel_nll(log_lambda: float, log_k: float, scores: np.ndarray, log_mn: np.ndarray) -> float:
    lam = math.exp(log_lambda)
    z = lam * scores - log_k - log_mn
    z = np.clip(z, -60.0, 60.0)
    return float(-np.sum(log_lambda - z - np.exp(-z)))


def nelder_mead(func, start: np.ndarray, steps: int = 120) -> np.ndarray:
    simplex = [start.copy()]
    for axis in range(start.size):
        vertex = start.copy()
        vertex[axis] += 0.15 if abs(start[axis]) < 1e-8 else 0.15 * abs(start[axis])
        simplex.append(vertex)
    simplex = [np.asarray(vertex, dtype=np.float64) for vertex in simplex]
    scores = [func(vertex) for vertex in simplex]
    for _ in range(steps):
        order = np.argsort(scores)
        simplex = [simplex[i] for i in order]
        scores = [scores[i] for i in order]
        centroid = np.mean(simplex[:-1], axis=0)
        reflected = centroid + (centroid - simplex[-1])
        reflected_score = func(reflected)
        if scores[0] <= reflected_score < scores[-2]:
            simplex[-1], scores[-1] = reflected, reflected_score
            continue
        if reflected_score < scores[0]:
            expanded = centroid + 2.0 * (reflected - centroid)
            expanded_score = func(expanded)
            simplex[-1], scores[-1] = (expanded, expanded_score) if expanded_score < reflected_score else (reflected, reflected_score)
            continue
        contracted = centroid + 0.5 * (simplex[-1] - centroid)
        contracted_score = func(contracted)
        if contracted_score < scores[-1]:
            simplex[-1], scores[-1] = contracted, contracted_score
            continue
        for index in range(1, len(simplex)):
            simplex[index] = simplex[0] + 0.5 * (simplex[index] - simplex[0])
            scores[index] = func(simplex[index])
    return simplex[int(np.argmin(scores))]


def fit_karlin_altschul(scores: np.ndarray, query_lengths: np.ndarray, target_lengths: np.ndarray) -> dict:
    positive = scores > 0
    if int(positive.sum()) < 20:
        raise ValueError(
            f"only {int(positive.sum())} positive null scores; "
            "the scoring system may be too harsh for EVD calibration"
        )
    fit_scores = scores[positive]
    log_mn = np.log(query_lengths[positive].astype(np.float64) * target_lengths[positive].astype(np.float64))
    loc, scale = fit_gumbel_mle(fit_scores)
    lam = 1.0 / scale
    # S - ln(mn)/λ ~ Gumbel(ln(K)/λ, 1/λ)
    for _ in range(6):
        residual = fit_scores - log_mn / lam
        loc_res, scale_res = fit_gumbel_mle(residual)
        lam = 1.0 / scale_res
    k_init = math.exp(lam * loc_res)

    def objective(theta: np.ndarray) -> float:
        return gumbel_nll(float(theta[0]), float(theta[1]), fit_scores, log_mn)

    theta = nelder_mead(objective, np.array([math.log(lam), math.log(max(k_init, 1e-12))]))
    lam = math.exp(float(theta[0]))
    k_value = math.exp(float(theta[1]))
    if not math.isfinite(lam) or not math.isfinite(k_value) or lam <= 0 or k_value <= 0:
        raise ValueError("Karlin–Altschul fit did not produce finite positive λ, K")

    z = lam * fit_scores - math.log(k_value) - log_mn
    # Kolmogorov–Smirnov against standard Gumbel: F(z) = exp(-exp(-z))
    ordered = np.sort(z)
    cdf = np.exp(-np.exp(-np.clip(ordered, -60.0, 60.0)))
    n = ordered.size
    empirical = np.arange(1, n + 1, dtype=np.float64) / n
    ks = float(np.max(np.abs(cdf - empirical)))
    return {
        "lambda": float(lam),
        "K": float(k_value),
        "gumbel_loc": float(loc),
        "gumbel_scale": float(scale),
        "n_positive": int(positive.sum()),
        "n_zero": int((~positive).sum()),
        "mean_score": float(fit_scores.mean()),
        "mean_log_mn": float(log_mn.mean()),
        "ks_statistic": ks,
    }


def database_residues(embeddings: dict[str, np.ndarray]) -> int:
    return int(sum(array.shape[0] for array in embeddings.values()))


def bit_score(raw_score: float, lam: float, k_value: float) -> float:
    return (lam * raw_score - math.log(k_value)) / LN2


def log_evalue(raw_score: float, query_length: int, search_residues: int, lam: float, k_value: float) -> float:
    return math.log(k_value) + math.log(max(query_length, 1)) + math.log(max(search_residues, 1)) - lam * raw_score


def evalue_from_log(log_e: float) -> float:
    if log_e > 700:
        return math.inf
    if log_e < -745:
        return 0.0
    return math.exp(log_e)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True, help="Packed faiss/ directory")
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--max-length", type=int, default=400)
    parser.add_argument("--max-cells", type=int, default=16_777_216)
    parser.add_argument("--max-alignments", type=int, default=16)
    parser.add_argument("--min-score", type=float, default=0.0)
    parser.add_argument("--min-match-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.samples < 20:
        print("error: --samples must be >= 20", file=sys.stderr)
        return 2
    if args.max_alignments < 1:
        print("error: --max-alignments must be >= 1", file=sys.stderr)
        return 2
    if not math.isfinite(args.min_score):
        print("error: --min-score must be finite", file=sys.stderr)
        return 2
    if args.min_match_count < 1:
        print("error: --min-match-count must be >= 1", file=sys.stderr)
        return 2
    try:
        if args.workers < 1:
            raise ValueError("--workers must be >= 1")
        embeddings = load_embeddings(args.database)
        params, scoring_sha256 = load_parameters(args.parameters)
        rng = np.random.default_rng(args.seed)
        scores, q_len, t_len = sample_null_scores(
            embeddings,
            params,
            args.samples,
            args.max_length,
            args.max_cells,
            args.max_alignments,
            args.min_score,
            args.min_match_count,
            rng,
            workers=args.workers,
        )
        fit = fit_karlin_altschul(scores, q_len, t_len)
    except (OSError, ValueError, TypeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    payload = {
        "lambda": fit["lambda"],
        "K": fit["K"],
        "null_model": "reverse_sequence",
        "fit_method": "length_aware_gumbel_mle",
        "score_definition": "sum_of_disjoint_local_hsp_scores",
        "max_alignments": args.max_alignments,
        "min_score": args.min_score,
        "min_match_count": args.min_match_count,
        "n_samples": int(scores.size),
        "n_positive": fit["n_positive"],
        "n_zero": fit["n_zero"],
        "database_sequences": len(embeddings),
        "database_residues": database_residues(embeddings),
        "max_length": args.max_length,
        "random_seed": args.seed,
        "scoring_sha256": scoring_sha256,
        "scoring_parameters": scoring_parameters_payload(params),
        "gumbel_loc": fit["gumbel_loc"],
        "gumbel_scale": fit["gumbel_scale"],
        "mean_score": fit["mean_score"],
        "mean_log_mn": fit["mean_log_mn"],
        "ks_statistic": fit["ks_statistic"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "lambda": payload["lambda"],
        "K": payload["K"],
        "n_positive": payload["n_positive"],
        "database_residues": payload["database_residues"],
        "ks_statistic": payload["ks_statistic"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
