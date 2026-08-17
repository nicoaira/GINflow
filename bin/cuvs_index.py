#!/usr/bin/env python3
"""GPU-only cuVS index support for CAGRA, IVF-Flat, and IVF-PQ."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np


CUVS_INDEX_TYPES = ("CAGRA", "IVF", "IVF-PQ")
_ALIASES = {
    "CAGRA": "CAGRA",
    "IVF": "IVF",
    "IVFFLAT": "IVF",
    "IVFPQ": "IVF-PQ",
}


def normalize_index_type(name: str) -> str:
    key = name.strip().replace("_", "").replace("-", "").upper()
    if key not in _ALIASES:
        allowed = ", ".join(kind.lower() for kind in CUVS_INDEX_TYPES)
        raise ValueError(f"unknown cuVS index type {name!r}. Choose one of: {allowed}")
    return _ALIASES[key]


def is_cuvs_type(name: str | None) -> bool:
    if not name:
        return False
    try:
        normalize_index_type(name)
    except ValueError:
        return False
    return True


def is_cuvs_database(directory: Path, meta: dict[str, Any] | None = None) -> bool:
    if str((meta or {}).get("backend") or "").strip().lower() == "cuvs":
        return True
    if is_cuvs_type(str((meta or {}).get("index_type") or "")):
        return True
    return (directory / "cuvs").is_dir()


def gpu_available() -> bool:
    try:
        import cupy as cp

        return int(cp.cuda.runtime.getDeviceCount()) > 0
    except Exception:
        return False


def require_gpu() -> None:
    try:
        import cupy  # noqa: F401
    except ImportError as exc:
        raise ValueError(
            "cuVS was requested but CuPy is not installed. Re-run with --index cuvs "
            "so BUILD_FAISS_INDEX and SEARCH_FAISS use the cuVS GPU image."
        ) from exc
    if not gpu_available():
        raise ValueError(
            "cuVS requires a visible NVIDIA GPU. Use -profile gpu on a host with CUDA."
        )


def resolve_n_lists(n_vectors: int, requested: int | None) -> int:
    if n_vectors < 1:
        raise ValueError("cannot build a cuVS IVF index with no vectors")
    if requested is not None:
        if requested < 1:
            raise ValueError("--cuvs_n_lists must be >= 1")
        if requested > n_vectors:
            raise ValueError(
                f"--cuvs_n_lists {requested} exceeds n_windows={n_vectors}. "
                "Use a smaller value."
            )
        return requested
    return max(1, min(n_vectors, int(4 * math.sqrt(n_vectors))))


def resolve_n_probes(n_lists: int, requested: int | None) -> int:
    if requested is not None and requested < 1:
        raise ValueError("--cuvs_n_probes must be >= 1")
    return max(1, min(n_lists, requested or min(20, n_lists)))


def resolve_degree(n_vectors: int, requested: int, name: str) -> int:
    if n_vectors < 2:
        raise ValueError("cuVS graph indexes require at least 2 vectors")
    if requested < 1:
        raise ValueError(f"--{name} must be >= 1")
    return max(1, min(n_vectors - 1, requested))


def _module(kind: str) -> Any:
    from cuvs.neighbors import cagra, ivf_flat, ivf_pq

    return {"CAGRA": cagra, "IVF": ivf_flat, "IVF-PQ": ivf_pq}[kind]


class CuvsIndex:
    """FAISS-shaped wrapper around a serialized cuVS index."""

    def __init__(
        self,
        index: Any,
        kind: str,
        ntotal: int,
        nprobe: int | None = None,
        itopk_size: int = 64,
    ) -> None:
        self.index = index
        self.index_type = normalize_index_type(kind)
        self.ntotal = int(ntotal)
        self.metric = "cuvs_cosine"
        self.nprobe = nprobe
        self.itopk_size = int(itopk_size)
        module = _module(self.index_type)
        if self.index_type == "CAGRA":
            self.search_params = module.SearchParams(itopk_size=self.itopk_size)
        else:
            self.search_params = module.SearchParams(n_probes=int(nprobe or 1))

    def search(self, queries: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
        import cupy as cp

        xb = np.ascontiguousarray(queries, dtype=np.float32)
        if xb.ndim == 1:
            xb = xb.reshape(1, -1)
        if xb.ndim != 2:
            raise ValueError(f"cuVS queries must be 1-D or 2-D, got {xb.shape}")
        search_k = max(1, min(int(k), self.ntotal))
        module = _module(self.index_type)
        search_params = self.search_params
        if self.index_type == "CAGRA" and search_k > self.itopk_size:
            search_params = module.SearchParams(itopk_size=search_k)
        distances, labels = module.search(
            search_params,
            self.index,
            cp.asarray(xb),
            search_k,
        )
        return (
            np.asarray(cp.asnumpy(cp.asarray(distances)), dtype=np.float32),
            np.asarray(cp.asnumpy(cp.asarray(labels)), dtype=np.int64),
        )

    def close(self) -> None:
        self.search_params = None
        self.index = None


def build_populated_index(
    vectors: np.ndarray,
    index_type: str = "cagra",
    *,
    n_lists: int | None = None,
    n_probes: int | None = None,
    pq_bits: int = 8,
    pq_dim: int = 0,
    intermediate_graph_degree: int = 128,
    graph_degree: int = 64,
    build_algo: str = "nn_descent",
    itopk_size: int = 64,
) -> tuple[CuvsIndex, dict[str, Any]]:
    require_gpu()
    xb = np.ascontiguousarray(vectors, dtype=np.float32)
    if xb.ndim != 2 or xb.shape[0] == 0:
        raise ValueError("need a non-empty 2-D vector matrix to build a cuVS index")
    kind = normalize_index_type(index_type)
    import cupy as cp

    dataset = cp.asarray(xb)
    resolved_lists = None
    resolved_probes = None
    if kind == "CAGRA":
        module = _module(kind)
        resolved_intermediate = resolve_degree(
            xb.shape[0], intermediate_graph_degree, "cuvs_intermediate_graph_degree"
        )
        resolved_graph = resolve_degree(xb.shape[0], graph_degree, "cuvs_graph_degree")
        if resolved_graph > resolved_intermediate:
            raise ValueError(
                "--cuvs_graph_degree cannot exceed --cuvs_intermediate_graph_degree"
            )
        params = module.IndexParams(
            metric="cosine",
            intermediate_graph_degree=resolved_intermediate,
            graph_degree=resolved_graph,
            build_algo=build_algo,
        )
        index = module.build(params, dataset)
        resolved_itopk = max(int(itopk_size), 1)
    elif kind == "IVF":
        module = _module(kind)
        resolved_lists = resolve_n_lists(xb.shape[0], n_lists)
        resolved_probes = resolve_n_probes(resolved_lists, n_probes)
        params = module.IndexParams(n_lists=resolved_lists, metric="cosine")
        index = module.build(params, dataset)
        resolved_itopk = 64
    else:
        module = _module(kind)
        resolved_lists = resolve_n_lists(xb.shape[0], n_lists)
        resolved_probes = resolve_n_probes(resolved_lists, n_probes)
        if pq_bits not in {4, 5, 6, 7, 8}:
            raise ValueError("--cuvs_pq_bits must be one of 4, 5, 6, or 7, or 8")
        params = module.IndexParams(
            n_lists=resolved_lists,
            metric="cosine",
            pq_bits=pq_bits,
            pq_dim=pq_dim,
        )
        index = module.build(params, dataset)
        resolved_itopk = 64

    wrapper = CuvsIndex(
        index,
        kind,
        xb.shape[0],
        nprobe=resolved_probes,
        itopk_size=resolved_itopk,
    )
    details = {
        "backend": "cuvs",
        "index_type": kind,
        "index_class": kind,
        "metric": "cuvs_cosine",
        "gpu": True,
        "gpu_capable": True,
        "n_lists": resolved_lists,
        "n_probes": resolved_probes,
        "pq_bits": pq_bits if kind == "IVF-PQ" else None,
        "pq_dim": pq_dim if kind == "IVF-PQ" else None,
        "intermediate_graph_degree": (
            params.intermediate_graph_degree if kind == "CAGRA" else None
        ),
        "graph_degree": params.graph_degree if kind == "CAGRA" else None,
        "build_algo": build_algo if kind == "CAGRA" else None,
        "itopk_size": resolved_itopk if kind == "CAGRA" else None,
    }
    return wrapper, details


def load_index(path: Path, meta: dict[str, Any], n_probes: int | None = None) -> CuvsIndex:
    kind = normalize_index_type(str(meta.get("index_type") or "cagra"))
    ntotal = meta.get("n_windows")
    if ntotal is None:
        raise ValueError("cuVS database metadata is missing n_windows")
    module = _module(kind)
    index = module.load(str(path / "index.bin"))
    resolved_probes = None
    if kind != "CAGRA":
        n_lists = meta.get("n_lists")
        if n_lists is None:
            raise ValueError("cuVS IVF metadata is missing n_lists")
        requested_probes = n_probes if n_probes is not None else meta.get("n_probes")
        resolved_probes = resolve_n_probes(int(n_lists), requested_probes)
    return CuvsIndex(
        index,
        kind,
        int(ntotal),
        nprobe=resolved_probes,
        itopk_size=int(meta.get("itopk_size") or 64),
    )


def serialize_index(index: CuvsIndex, destination: Path) -> None:
    import cuvs  # noqa: F401

    destination.mkdir(parents=True, exist_ok=True)
    module = _module(index.index_type)
    module.save(str(destination / "index.bin"), index.index, include_dataset=True)
    index.close()


def meta_from_details(details: dict[str, Any]) -> dict[str, Any]:
    return {
        "backend": "cuvs",
        "index_type": details["index_type"],
        "index_class": details["index_class"],
        "metric": details["metric"],
        "gpu": True,
        "gpu_capable": True,
        "n_lists": details.get("n_lists"),
        "n_probes": details.get("n_probes"),
        "pq_bits": details.get("pq_bits"),
        "pq_dim": details.get("pq_dim"),
        "intermediate_graph_degree": details.get("intermediate_graph_degree"),
        "graph_degree": details.get("graph_degree"),
        "build_algo": details.get("build_algo"),
        "itopk_size": details.get("itopk_size"),
    }
