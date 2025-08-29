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
        help="Output TSV with added z_score, p_value, and q_value (FDR) columns",
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

    # Benjamini–Hochberg FDR (q-values)
    def fdr_bh(pvals: pd.Series) -> pd.Series:
        # Treat missing as 1.0 (conservative)
        s = pvals.fillna(1.0)
        m = len(s)
        if m == 0:
            return s
        # Sort ascending, compute adjusted, enforce monotonicity
        sorted_idx = s.sort_values(kind="mergesort").index
        sorted_p = s.loc[sorted_idx].tolist()
        adj = [min(p * m / (i + 1), 1.0) for i, p in enumerate(sorted_p)]
        for i in range(m - 2, -1, -1):
            adj[i] = min(adj[i], adj[i + 1])
        # Map back to original order
        q_map = {idx: val for idx, val in zip(sorted_idx, adj)}
        return s.index.to_series().map(q_map)

    scores_df["q_value"] = fdr_bh(scores_df["p_value"])  # FDR-adjusted p-values

    scores_df.to_csv(args.output, sep="\t", index=False)
    print(f"Written normalized scores with FDR to {args.output}")


if __name__ == "__main__":
    main()
