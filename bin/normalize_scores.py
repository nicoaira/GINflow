#!/usr/bin/env python3
import argparse
import pandas as pd
import math


def norm_sf(z: float) -> float:
    """Survival function (1 - CDF) for standard normal.

    Uses erfc for numerical stability in the extreme tail, avoiding
    catastrophic cancellation that makes small p-values become 0.0.
    """
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def main():
    parser = argparse.ArgumentParser(
        description="Normalize aggregated scores against a null distribution"
    )
    parser.add_argument(
        "--scores",
        required=True,
        help="TSV file with a 'score' column to normalize",
    )
    parser.add_argument(
        "--null-distribution",
        required=True,
        help="TSV file containing a 'score' column representing the null distribution",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV with added z_score and p_value columns",
    )
    args = parser.parse_args()

    scores_df = pd.read_csv(args.scores, sep="\t")
    null_df = pd.read_csv(args.null_distribution, sep="\t")

    if "score" not in scores_df.columns:
        raise ValueError("scores file must contain a 'score' column")
    if "score" not in null_df.columns:
        raise ValueError("null distribution file must contain a 'score' column")

    mu = null_df["score"].mean()
    sigma = null_df["score"].std(ddof=0)
    if sigma == 0:
        raise ValueError("null distribution standard deviation is zero")

    z = (scores_df["score"] - mu) / sigma
    scores_df["z_score"] = z
    scores_df["p_value"] = z.map(norm_sf)

    scores_df.to_csv(args.output, sep="\t", index=False)
    print(f"Written normalized scores to {args.output}")


if __name__ == "__main__":
    main()
