#!/usr/bin/env python3
"""ScaNN searcher construction, serialization, and FAISS-shaped search."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from faiss_index import IndexOptions

SCANN_DIRNAME = "scann"
BRUTE_FORCE_MAX = 20_000
AH_ONLY_MAX = 100_000
DEFAULT_REORDER = 100
DEFAULT_AH_DIM = 2
DEFAULT_ANISOTROPIC = 0.2
DEFAULT_NUM_NEIGHBORS = 100
SOAR_LAMBDA = 1.5


def is_scann_type(name: str | None) -> bool:
    if not name:
        return False
    key = name.strip().replace("-", "").replace("_", "").upper()
    return key in {"SCANN", "SCAN"}


def is_scann_database(directory: Path, meta: dict[str, Any] | None = None) -> bool:
    if is_scann_type(str((meta or {}).get("index_type") or "")):
        return True
    if str((meta or {}).get("backend") or "").strip().lower() == "scann":
        return True
    return (directory / SCANN_DIRNAME).is_dir()


def _import_scann():
    try:
        import scann
    except ImportError as exc:
        raise ValueError(
            "ScaNN was requested but the scann package is not installed. "
            "Re-run with --index scann so BUILD_FAISS_INDEX and SEARCH_FAISS "
            "use the ScaNN image."
        ) from exc
    return scann


class ScannIndex:
    """Minimal FAISS-compatible wrapper: ``ntotal`` and ``search(xq, k)``."""

    def __init__(
        self,
        searcher: Any,
        ntotal: int,
        *,
        leaves_to_search: int | None = None,
        reorder: int | None = None,
    ) -> None:
        self.searcher = searcher
        self.ntotal = int(ntotal)
        self.leaves_to_search = leaves_to_search
        self.reorder = reorder

    def search(self, queries: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
        xb = np.ascontiguousarray(queries, dtype=np.float32)
        if xb.ndim == 1:
            xb = xb.reshape(1, -1)
        if xb.ndim != 2:
            raise ValueError(f"ScaNN queries must be 1-D or 2-D, got {xb.shape}")
        search_k = max(1, min(int(k), self.ntotal))
        kwargs: dict[str, Any] = {"final_num_neighbors": search_k}
        if self.reorder is not None:
            kwargs["pre_reorder_num_neighbors"] = max(int(self.reorder), search_k)
        if self.leaves_to_search is not None:
            kwargs["leaves_to_search"] = int(self.leaves_to_search)
        idx, dist = self.searcher.search_batched(xb, **kwargs)
        return _faiss_shaped(idx, dist, xb.shape[0], search_k)


def _faiss_shaped(
    idx: Any,
    dist: Any,
    n_queries: int,
    k: int,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.full((n_queries, k), -1, dtype=np.int64)
    distances = np.full((n_queries, k), np.float32("-inf"), dtype=np.float32)
    labels_in = np.asarray(idx)
    if labels_in.dtype == object or (isinstance(idx, (list, tuple)) and idx and not isinstance(idx[0], (int, np.integer))):
        for i, (row_idx, row_dist) in enumerate(zip(idx, dist)):
            if i >= n_queries:
                break
            row_i = np.asarray(row_idx, dtype=np.int64).reshape(-1)
            row_d = np.asarray(row_dist, dtype=np.float32).reshape(-1)
            take = min(k, row_i.size, row_d.size)
            labels[i, :take] = row_i[:take]
            distances[i, :take] = row_d[:take]
        return distances, labels
    scores_in = np.asarray(dist, dtype=np.float32)
    if labels_in.size == 0:
        return distances, labels
    if labels_in.ndim == 1:
        labels_in = labels_in.reshape(1, -1)
        scores_in = scores_in.reshape(1, -1)
    take = min(k, labels_in.shape[1], scores_in.shape[1])
    rows = min(n_queries, labels_in.shape[0])
    labels[:rows, :take] = labels_in[:rows, :take]
    distances[:rows, :take] = scores_in[:rows, :take]
    return distances, labels


def resolve_num_leaves(n_vectors: int, requested: int | None) -> int:
    if n_vectors < 1:
        raise ValueError("cannot partition a ScaNN index with no vectors")
    if requested is not None:
        if requested < 1:
            raise ValueError("--scann_leaves must be >= 1")
        if requested > n_vectors:
            raise ValueError(
                f"--scann_leaves {requested} exceeds n_windows={n_vectors}. "
                "Use a smaller value or leave it unset for sqrt(n) leaves."
            )
        return requested
    return max(1, min(n_vectors, int(round(math.sqrt(n_vectors)))))


def resolve_leaves_to_search(num_leaves: int, requested: int | None) -> int:
    if requested is None:
        return max(1, min(8, num_leaves))
    if requested < 1:
        raise ValueError("--scann_leaves_to_search must be >= 1")
    return min(requested, num_leaves)


def resolve_scann_plan(n_vectors: int, options: IndexOptions) -> dict[str, Any]:
    force_tree = options.scann_leaves is not None
    if n_vectors < BRUTE_FORCE_MAX and not force_tree:
        return {
            "scoring": "brute_force",
            "partitioned": False,
            "num_leaves": None,
            "leaves_to_search": None,
        }
    use_tree = force_tree or n_vectors >= AH_ONLY_MAX
    num_leaves = resolve_num_leaves(n_vectors, options.scann_leaves) if use_tree else None
    leaves_to_search = (
        resolve_leaves_to_search(num_leaves, options.scann_leaves_to_search) if num_leaves else None
    )
    return {
        "scoring": "ah",
        "partitioned": use_tree,
        "num_leaves": num_leaves,
        "leaves_to_search": leaves_to_search,
    }


def build_populated_searcher(
    vectors: np.ndarray,
    options: IndexOptions,
) -> tuple[ScannIndex, dict[str, Any]]:
    if options.gpu:
        raise ValueError(
            "--faiss_gpu is not supported for ScaNN. ScaNN is CPU-only (AVX/FMA)."
        )
    xb = np.ascontiguousarray(vectors, dtype=np.float32)
    if xb.ndim != 2 or xb.shape[0] == 0:
        raise ValueError("need a non-empty 2-D vector matrix to build a ScaNN index")
    if options.scann_ah_dim < 1:
        raise ValueError("--scann_ah_dim must be >= 1")
    if options.scann_reorder < 1:
        raise ValueError("--scann_reorder must be >= 1")
    if options.scann_num_neighbors < 1:
        raise ValueError("--seed_k / ScaNN num_neighbors must be >= 1")

    scann = _import_scann()
    n_vectors, _dim = xb.shape
    num_neighbors = int(options.scann_num_neighbors)
    reorder = max(int(options.scann_reorder), num_neighbors)
    plan = resolve_scann_plan(n_vectors, options)

    builder = scann.scann_ops_pybind.builder(xb, num_neighbors, "dot_product")
    if plan["partitioned"]:
        tree_kwargs: dict[str, Any] = {
            "num_leaves": plan["num_leaves"],
            "num_leaves_to_search": plan["leaves_to_search"],
            "spherical": True,
        }
        min_partition = max(1, n_vectors // int(plan["num_leaves"]))
        if min_partition < 50:
            tree_kwargs["min_partition_size"] = min_partition
        if options.scann_soar:
            tree_kwargs["soar_lambda"] = SOAR_LAMBDA
        builder = builder.tree(**tree_kwargs)

    if plan["scoring"] == "brute_force":
        builder = builder.score_brute_force()
    else:
        builder = builder.score_ah(
            int(options.scann_ah_dim),
            anisotropic_quantization_threshold=float(options.scann_anisotropic),
        ).reorder(reorder)

    searcher = builder.build()
    index = ScannIndex(
        searcher,
        n_vectors,
        leaves_to_search=plan["leaves_to_search"],
        reorder=reorder if plan["scoring"] == "ah" else None,
    )
    details = {
        "index_type": "ScaNN",
        "index_class": "ScannSearcher",
        "metric": "inner_product",
        "gpu": False,
        "gpu_capable": False,
        "backend": "scann",
        "nlist": plan["num_leaves"],
        "nprobe": plan["leaves_to_search"],
        "pq_m": None,
        "pq_nbits": None,
        "pq_m_refine": None,
        "hnsw_m": None,
        "hnsw_ef_construction": None,
        "hnsw_ef_search": None,
        "lsh_nbits": None,
        "sq_type": None,
        "scann_scoring": plan["scoring"],
        "scann_partitioned": plan["partitioned"],
        "scann_reorder": reorder if plan["scoring"] == "ah" else None,
        "scann_ah_dim": int(options.scann_ah_dim) if plan["scoring"] == "ah" else None,
        "scann_anisotropic": (
            float(options.scann_anisotropic) if plan["scoring"] == "ah" else None
        ),
        "scann_soar": bool(options.scann_soar) if plan["partitioned"] else False,
        "scann_num_neighbors": num_neighbors,
    }
    return index, details


def serialize_index(index: ScannIndex, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    # relative_path=True stores asset names (dataset.npy) instead of the build
    # workdir path, so SEARCH_FAISS can load the staged copy.
    index.searcher.serialize(str(directory.resolve()), True)


def load_index(
    directory: Path,
    ntotal: int | None = None,
    *,
    leaves_to_search: int | None = None,
    reorder: int | None = None,
) -> ScannIndex:
    if not directory.is_dir():
        raise ValueError(f"ScaNN artifacts directory not found: {directory}")
    scann = _import_scann()
    searcher = scann.scann_ops_pybind.load_searcher(str(directory))
    size = ntotal
    if size is None:
        getter = getattr(searcher, "size", None)
        size = int(getter()) if getter is not None else 0
    return ScannIndex(
        searcher,
        int(size),
        leaves_to_search=leaves_to_search,
        reorder=reorder,
    )


def apply_search_params(
    index: ScannIndex,
    *,
    nprobe: int | None = None,
    reorder: int | None = None,
) -> ScannIndex:
    if nprobe is not None:
        if nprobe < 1:
            raise ValueError("--scann_leaves_to_search must be >= 1")
        index.leaves_to_search = int(nprobe)
    if reorder is not None:
        if reorder < 1:
            raise ValueError("--scann_reorder must be >= 1")
        index.reorder = int(reorder)
    return index
