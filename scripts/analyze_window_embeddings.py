#!/usr/bin/env python3
"""Analyse database window embeddings for diversity and create diagnostic plots."""

from __future__ import annotations

import argparse
import json
import logging
import math
import warnings
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Inspect window embeddings exported by the pipeline, estimate their diversity, "
            "and produce plots that summarise the distribution."
        )
    )
    parser.add_argument(
        "--embeddings",
        type=Path,
        required=True,
        help="Path to database_windows.npy (rows correspond to window embeddings).",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Path to database_windows.tsv describing each window.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("window_embedding_report"),
        help="Directory where plots and summary.json are written (created if missing).",
    )
    parser.add_argument(
        "--label-column",
        type=str,
        default="transcript_id",
        help=(
            "Metadata column used to colour PCA scatter points and compute per-label dispersion. "
            "Falls back to unlabelled behaviour if the column is absent."
        ),
    )
    parser.add_argument(
        "--pairwise-sample-size",
        type=int,
        default=2000,
        help="Number of windows sampled to estimate pairwise cosine similarities (set <= 0 to use all).",
    )
    parser.add_argument(
        "--pca-sample-size",
        type=int,
        default=2000,
        help="Number of windows sampled for PCA scatter (set <= 0 to use all).",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=50000,
        help="Maximum number of windows sampled for overall analysis (set <= 0 to use all).",
    )
    parser.add_argument(
        "--dispersion-sample-size",
        type=int,
        default=20000,
        help="Number of windows sampled when estimating per-label dispersion (set <= 0 to use all sampled windows).",
    )
    parser.add_argument(
        "--dispersion-top",
        type=int,
        default=10,
        help="Maximum number of labels to include in the dispersion bar plot.",
    )
    parser.add_argument(
        "--hist-bins",
        type=int,
        default=60,
        help="Histogram bins for the pairwise cosine similarity plot.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=7,
        help="Random seed used for sampling windows.",
    )
    return parser.parse_args()


def _load_embeddings(path: Path) -> np.ndarray:
    """Load embeddings via memory mapping to limit RAM usage."""
    embeddings = np.load(path, mmap_mode="r")
    if embeddings.ndim != 2:
        raise ValueError(f"Expected 2D array in {path}, found shape {embeddings.shape!r}")
    return embeddings


def _load_metadata(path: Path, columns: Sequence[str] | None = None) -> pd.DataFrame:
    """Load the metadata TSV."""
    read_kwargs: dict[str, object] = {"sep": "\t"}
    if columns is not None:
        read_kwargs["usecols"] = list(columns)
    return pd.read_csv(path, **read_kwargs)


def _safe_tight_layout(fig: plt.Figure) -> None:
    """Attempt tight_layout while silencing benign warnings."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        fig.tight_layout()


def _sample_indices(
    total: int,
    requested: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample unique row indices."""
    if requested <= 0 or requested >= total:
        return np.arange(total)
    return rng.choice(total, size=requested, replace=False)


def _select_analysis_indices(
    total: int,
    max_windows: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Select the subset of windows to keep for the downstream analysis."""
    if max_windows <= 0 or max_windows >= total:
        return np.arange(total)
    return np.sort(rng.choice(total, size=max_windows, replace=False))


def _pairwise_cosine(sample_vectors: np.ndarray) -> np.ndarray:
    """Compute cosine similarities for all unique pairs in the sample."""
    if sample_vectors.shape[0] < 2:
        return np.empty(0, dtype=np.float64)
    gram = sample_vectors @ sample_vectors.T
    iu = np.triu_indices(sample_vectors.shape[0], k=1)
    return gram[iu]


def _run_pca(vectors: np.ndarray, components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Project vectors to the leading principal components."""
    if vectors.shape[0] < 2:
        raise ValueError("PCA needs at least two vectors.")

    mean = vectors.mean(axis=0, dtype=np.float64, keepdims=True)
    centered = vectors - mean
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    coords = centered @ vt[:components].T
    explained = (s**2) / max(vectors.shape[0] - 1, 1)
    total_var = explained.sum()
    if total_var <= 0:
        ratios = np.zeros(components, dtype=np.float64)
    else:
        ratios = explained[:components] / total_var
    return coords, ratios


def _plot_similarity_hist(
    similarities: np.ndarray,
    bins: int,
    path: Path,
) -> None:
    """Plot histogram of sampled pairwise cosine similarities."""
    if similarities.size == 0:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(similarities, bins=bins, color="#4C72B0", alpha=0.85)
    ax.set_xlabel("Cosine similarity")
    ax.set_ylabel("Frequency")
    ax.set_title("Pairwise cosine similarities (sampled windows)")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    _safe_tight_layout(fig)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_pca_scatter(
    coords: np.ndarray,
    labels: Sequence[str],
    ratios: Sequence[float],
    path: Path,
) -> None:
    """Plot PCA scatter coloured by categorical labels."""
    fig, ax = plt.subplots(figsize=(6.4, 5.2))

    labels = np.asarray(labels, dtype=object)
    ratios_arr = np.asarray(ratios, dtype=float)
    unique, inverse, counts = np.unique(labels, return_inverse=True, return_counts=True)
    order = np.argsort(-counts)
    unique = unique[order]
    inverse = np.array([np.where(order == idx)[0][0] for idx in inverse])

    cmap = plt.get_cmap("tab20", len(unique))
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=inverse,
        cmap=cmap,
        s=12,
        linewidths=0,
        alpha=0.75,
    )
    ax.set_xlabel(
        f"PC1 ({ratios_arr[0] * 100:.1f}% var)" if ratios_arr.size > 0 else "PC1"
    )
    ax.set_ylabel(
        f"PC2 ({ratios_arr[1] * 100:.1f}% var)" if ratios_arr.size > 1 else "PC2"
    )
    ax.set_title("PCA projection of sampled windows")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.3)

    handles, legend_labels = scatter.legend_elements(num=len(unique))
    legend = ax.legend(
        handles,
        [f"{label} ({counts[idx]})" for idx, label in enumerate(unique)],
        loc="best",
        fontsize="small",
        title="Label (count)",
        frameon=True,
    )
    legend.get_title().set_fontsize("small")

    _safe_tight_layout(fig)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _compute_label_dispersion(
    vectors: np.ndarray,
    labels: Iterable[str],
) -> list[dict[str, float | int | str]]:
    """Estimate how tightly grouped embeddings are for each label."""
    label_series = pd.Series(labels, dtype="object")
    groups = label_series.groupby(label_series).groups
    results: list[dict[str, float | int | str]] = []
    for label, idx in groups.items():
        indices = np.asarray(idx, dtype=int)
        if indices.size < 2:
            continue
        subset = vectors[indices]
        centroid = subset.mean(axis=0, dtype=np.float64)
        norm = np.linalg.norm(centroid)
        if not math.isfinite(norm) or norm == 0:
            continue
        centroid /= norm
        cosines = subset @ centroid
        dispersion = 1.0 - float(np.mean(cosines))
        results.append({"label": str(label), "count": int(indices.size), "dispersion": dispersion})
    results.sort(key=lambda item: (item["dispersion"], -item["count"]))
    return results


def _plot_dispersion_bars(
    dispersion_stats: Sequence[dict[str, float | int | str]],
    top_k: int,
    path: Path,
) -> None:
    """Plot bar chart for per-label dispersion."""
    if not dispersion_stats:
        return
    subset = dispersion_stats[:top_k]
    labels = [item["label"] for item in subset]
    values = [item["dispersion"] for item in subset]
    counts = [item["count"] for item in subset]

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    bars = ax.bar(range(len(labels)), values, color="#55A868", alpha=0.85)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize="small")
    ax.set_ylabel("1 - mean cosine to centroid")
    ax.set_title("Within-label dispersion (lower = tighter)")
    ax.grid(True, axis="y", linestyle="--", linewidth=0.4, alpha=0.3)
    for bar, count in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"n={count}",
            ha="center",
            va="bottom",
            fontsize="x-small",
            rotation=0,
        )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    rng = np.random.default_rng(args.random_seed)

    logging.info("Loading embeddings from %s", args.embeddings)
    embeddings = _load_embeddings(args.embeddings)
    total_windows, embedding_dim = embeddings.shape
    logging.info("Loaded embeddings with shape %d x %d", total_windows, embedding_dim)

    logging.info("Loading metadata from %s", args.metadata)
    try:
        metadata = _load_metadata(args.metadata, columns=[args.label_column])
    except ValueError:
        logging.warning(
            "Column %r not found when reading metadata selectively; loading full table.",
            args.label_column,
        )
        metadata = _load_metadata(args.metadata)
    logging.info("Loaded %d metadata rows", len(metadata))

    if len(metadata) != total_windows:
        raise ValueError(
            f"Metadata rows ({len(metadata)}) do not match embeddings ({total_windows})."
        )

    analysis_indices = _select_analysis_indices(total_windows, args.max_windows, rng)
    analysis_vectors = np.asarray(embeddings[analysis_indices], dtype=np.float32)
    analysis_metadata = metadata.iloc[analysis_indices].reset_index(drop=True)
    logging.info(
        "Using %d/%d windows for analysis",
        analysis_vectors.shape[0],
        total_windows,
    )

    logging.info("Preparing output directory at %s", args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    norms = np.linalg.norm(analysis_vectors, axis=1)
    logging.info(
        "Sampled embedding L2 norms: min=%.4f max=%.4f mean=%.4f std=%.4f",
        float(norms.min()),
        float(norms.max()),
        float(norms.mean()),
        float(norms.std()),
    )
    summary: dict[str, object] = {
        "num_windows": int(total_windows),
        "analysis_sample_size": int(analysis_vectors.shape[0]),
        "embedding_dimension": int(embedding_dim),
        "norm_min": float(norms.min()),
        "norm_max": float(norms.max()),
        "norm_mean": float(norms.mean()),
        "norm_std": float(norms.std()),
    }

    pairwise_indices = _sample_indices(
        analysis_vectors.shape[0],
        args.pairwise_sample_size,
        rng,
    )
    pairwise_vectors = analysis_vectors[pairwise_indices].astype(np.float64, copy=False)
    logging.info(
        "Computing pairwise cosine similarities on %d sampled windows",
        pairwise_vectors.shape[0],
    )
    pairwise_similarities = _pairwise_cosine(pairwise_vectors)
    summary["pairwise_sample_size"] = int(pairwise_vectors.shape[0])
    summary["pairwise_similarity_count"] = int(pairwise_similarities.size)
    if pairwise_similarities.size:
        summary["pairwise_similarity_mean"] = float(pairwise_similarities.mean())
        summary["pairwise_similarity_std"] = float(pairwise_similarities.std())
        summary["pairwise_similarity_min"] = float(pairwise_similarities.min())
        summary["pairwise_similarity_max"] = float(pairwise_similarities.max())
    logging.info(
        "Pairwise similarities computed (%d values)", pairwise_similarities.size
    )

    similarity_hist_path = args.outdir / "pairwise_cosine_hist.png"
    logging.info("Rendering histogram to %s", similarity_hist_path)
    _plot_similarity_hist(pairwise_similarities, args.hist_bins, similarity_hist_path)

    pca_indices = _sample_indices(
        analysis_vectors.shape[0],
        args.pca_sample_size,
        rng,
    )
    pca_vectors = analysis_vectors[pca_indices].astype(np.float64, copy=False)
    if pca_vectors.shape[0] >= 2:
        logging.info(
            "Running PCA on %d sampled windows of dimension %d",
            pca_vectors.shape[0],
            pca_vectors.shape[1],
        )
        pca_coords, ratios = _run_pca(pca_vectors, components=2)
        label_col = args.label_column
        if label_col in analysis_metadata.columns:
            labels = analysis_metadata.iloc[pca_indices][label_col].astype(str).to_list()
        else:
            labels = ["unlabelled"] * pca_vectors.shape[0]
        pca_path = args.outdir / "pca_scatter.png"
        logging.info("Rendering PCA scatter plot to %s", pca_path)
        _plot_pca_scatter(pca_coords, labels, ratios, pca_path)
        summary["pca_sample_size"] = int(pca_vectors.shape[0])
        summary["pca_explained_variance_ratio"] = [float(r) for r in ratios[:2]]
    else:
        summary["pca_sample_size"] = int(pca_vectors.shape[0])

    label_col = args.label_column
    if label_col in analysis_metadata.columns:
        dispersion_indices = _sample_indices(
            analysis_vectors.shape[0],
            args.dispersion_sample_size,
            rng,
        )
        dispersion_vectors = analysis_vectors[dispersion_indices].astype(
            np.float32, copy=False
        )
        dispersion_labels = (
            analysis_metadata.iloc[dispersion_indices][label_col].astype(str).to_list()
        )
        summary["dispersion_sample_size"] = int(dispersion_vectors.shape[0])
        dispersion_stats = _compute_label_dispersion(
            dispersion_vectors,
            dispersion_labels,
        )
        if dispersion_stats:
            logging.info(
                "Computed dispersion for %d labels, exporting top %d",
                len(dispersion_stats),
                args.dispersion_top,
            )
            summary["dispersion_top"] = dispersion_stats[: args.dispersion_top]
            dispersion_path = args.outdir / "label_dispersion.png"
            logging.info("Rendering dispersion bar chart to %s", dispersion_path)
            _plot_dispersion_bars(dispersion_stats, args.dispersion_top, dispersion_path)
        else:
            summary["dispersion_top"] = []
    else:
        summary["dispersion_sample_size"] = 0
        summary["dispersion_top"] = []

    summary_path = args.outdir / "summary.json"
    logging.info("Writing summary to %s", summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    logging.info("Analysis complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
