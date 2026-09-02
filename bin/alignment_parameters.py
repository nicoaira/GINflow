"""Build GINFINITY-SW scoring parameters from command-line values."""
from __future__ import annotations

import argparse
import hashlib
import json

from ginfinity_sw import ScoringParameters


def add_scoring_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mu", type=float, required=True)
    parser.add_argument("--sigma", type=float, required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--score-min", type=float, required=True)
    parser.add_argument("--score-max", type=float, required=True)
    parser.add_argument("--gap-open", type=float, required=True)
    parser.add_argument("--gap-extend", type=float, required=True)
    parser.add_argument("--score-offset", type=float, required=True)


def scoring_parameters_from_args(args: argparse.Namespace) -> ScoringParameters:
    return ScoringParameters(
        mu=args.mu,
        sigma=args.sigma,
        gamma=args.gamma,
        score_min=args.score_min,
        score_max=args.score_max,
        gap_open=args.gap_open,
        gap_extend=args.gap_extend,
        score_offset=args.score_offset,
    )


def scoring_parameters_payload(params: ScoringParameters) -> dict[str, float]:
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


def scoring_parameters_digest(params: ScoringParameters) -> str:
    serialized = json.dumps(
        scoring_parameters_payload(params), sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(serialized).hexdigest()
