#!/usr/bin/env python3
"""FAISS index constructors, GPU helpers, and score conversion for GINflow."""
from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import faiss
except ImportError:  # pragma: no cover - ScaNN image has no FAISS
    faiss = None  # type: ignore[assignment]

import numpy as np


INDEX_TYPES = (
    "FlatIP",
    "FlatL2",
    "HNSW",
    "IVFFlat",
)

# Classic GPU FAISS implements Flat and IVFFlat. HNSW stays on CPU
# (https://github.com/facebookresearch/faiss/wiki/Faiss-indexes).
GPU_INDEX_TYPES = frozenset({"FlatIP", "FlatL2", "IVFFlat"})

_ALIASES = {
    "FLAT": "FlatIP",
    "FLATIP": "FlatIP",
    "INDEXFLATIP": "FlatIP",
    "FLATL2": "FlatL2",
    "INDEXFLATL2": "FlatL2",
    "HNSW": "HNSW",
    "HNSWFLAT": "HNSW",
    "INDEXHNSWFLAT": "HNSW",
    "IVFFLAT": "IVFFlat",
    "INDEXIVFFLAT": "IVFFlat",
}

_SQ_TYPES = {
    "8": "8bit",
    "8bit": "8bit",
    "qt_8bit": "8bit",
    "6": "6bit",
    "6bit": "6bit",
    "qt_6bit": "6bit",
    "4": "4bit",
    "4bit": "4bit",
    "qt_4bit": "4bit",
    "fp16": "fp16",
    "qt_fp16": "fp16",
}

DEFAULT_NPROBE = 8
DEFAULT_PQ_M = 16
DEFAULT_PQ_NBITS = 8
DEFAULT_PQ_M_REFINE = 4
DEFAULT_HNSW_M = 32
DEFAULT_HNSW_EF_CONSTRUCTION = 40
DEFAULT_HNSW_EF_SEARCH = 16


@dataclass(frozen=True)
class IndexOptions:
    index_type: str = "flatip"
    nlist: int | None = None
    nprobe: int | None = None
    pq_m: int = DEFAULT_PQ_M
    pq_nbits: int = DEFAULT_PQ_NBITS
    pq_m_refine: int = DEFAULT_PQ_M_REFINE
    hnsw_m: int = DEFAULT_HNSW_M
    hnsw_ef_construction: int = DEFAULT_HNSW_EF_CONSTRUCTION
    hnsw_ef_search: int | None = DEFAULT_HNSW_EF_SEARCH
    lsh_nbits: int | None = None
    sq_type: str = "8bit"
    gpu: bool = False
    gpu_device: int = 0
    scann_reorder: int = 100
    scann_ah_dim: int = 2
    scann_anisotropic: float = 0.2
    scann_soar: bool = False
    scann_num_neighbors: int = 100
    scann_leaves: int | None = None
    scann_leaves_to_search: int | None = None


def normalize_index_type(name: str) -> str:
    key = name.strip().replace("-", "").replace("_", "").replace(",", "").upper()
    if key not in _ALIASES:
        allowed = ", ".join(kind.lower() for kind in INDEX_TYPES)
        raise ValueError(
            f"unknown FAISS index type {name!r}. Choose one of: {allowed}"
        )
    return _ALIASES[key]


def supports_gpu(index_type: str) -> bool:
    return normalize_index_type(index_type) in GPU_INDEX_TYPES


def gpu_runtime_available() -> bool:
    return faiss is not None and hasattr(faiss, "StandardGpuResources")


def cuda_device_visible() -> bool:
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible in {"", "-1"}:
        return False
    return Path("/dev/nvidiactl").exists() or Path("/dev/nvidia0").exists()


def gpu_device_count() -> int:
    if not cuda_device_visible():
        return 0
    getter = getattr(faiss, "get_num_gpus", None) if faiss is not None else None
    if getter is None:
        return 0
    try:
        return int(getter())
    except Exception:
        return 0


def require_gpu(index_type: str) -> None:
    kind = normalize_index_type(index_type)
    if kind not in GPU_INDEX_TYPES:
        supported = ", ".join(sorted(kind.lower() for kind in GPU_INDEX_TYPES))
        raise ValueError(
            f"--faiss_gpu is not supported for index type {kind.lower()}. "
            f"GPU indexes: {supported}. Use CPU for {kind.lower()}."
        )
    if not gpu_runtime_available():
        raise ValueError(
            "FAISS GPU was requested but this environment is faiss-cpu. "
            "Re-run with -profile gpu so BUILD_FAISS_INDEX and SEARCH_FAISS "
            "use the faiss-gpu image."
        )
    if gpu_device_count() < 1:
        raise ValueError(
            "FAISS GPU was requested but no CUDA device is visible. "
            "Use -profile gpu on a machine with an NVIDIA GPU."
        )


def resolve_nlist(n_vectors: int, requested: int | None) -> int:
    if n_vectors < 1:
        raise ValueError("cannot build an IVF index with no vectors")
    if requested is not None:
        if requested < 1:
            raise ValueError("--faiss_nlist must be >= 1")
        if requested > n_vectors:
            raise ValueError(
                f"--faiss_nlist {requested} exceeds n_windows={n_vectors}. "
                "Use a smaller nlist or leave it unset for an automatic value."
            )
        return requested
    auto = max(1, int(4 * math.sqrt(n_vectors)))
    # FAISS IVF training wants several points per centroid when possible.
    trained = max(1, n_vectors // 39) if n_vectors >= 39 else n_vectors
    return max(1, min(auto, trained, n_vectors))


def resolve_nprobe(nlist: int, requested: int | None) -> int:
    if requested is None:
        return max(1, min(DEFAULT_NPROBE, nlist))
    if requested < 1:
        raise ValueError("--faiss_nprobe must be >= 1")
    return min(requested, nlist)


def resolve_pq_m(dim: int, requested: int) -> int:
    if requested < 1:
        raise ValueError("--faiss_pq_m must be >= 1")
    if dim % requested == 0:
        return requested
    divisors = [i for i in range(1, dim + 1) if dim % i == 0]
    preview = ", ".join(str(v) for v in divisors[:16])
    extra = ", ..." if len(divisors) > 16 else ""
    raise ValueError(
        f"--faiss_pq_m {requested} must divide window_dim={dim}. "
        f"Valid values: {preview}{extra}"
    )


def resolve_sq_type(name: str) -> str:
    key = name.strip().lower().replace("-", "")
    if key not in _SQ_TYPES:
        raise ValueError(
            f"unknown --faiss_sq_type {name!r}. Choose 8bit, 6bit, 4bit, or fp16"
        )
    return _SQ_TYPES[key]


def _sq_enum(name: str) -> int:
    mapping = {
        "8bit": faiss.ScalarQuantizer.QT_8bit,
        "6bit": faiss.ScalarQuantizer.QT_6bit,
        "4bit": faiss.ScalarQuantizer.QT_4bit,
        "fp16": faiss.ScalarQuantizer.QT_fp16,
    }
    return mapping[resolve_sq_type(name)]


def default_lsh_nbits(dim: int, requested: int | None) -> int:
    if requested is None:
        return 2 * dim
    if requested < 1:
        raise ValueError("--faiss_lsh_nbits must be >= 1")
    return requested


def make_cpu_index(dim: int, n_vectors: int, options: IndexOptions) -> tuple[Any, dict[str, Any]]:
    if faiss is None:
        raise ValueError(
            "FAISS is not installed in this environment. Use the FAISS container."
        )
    kind = normalize_index_type(options.index_type)
    if dim < 1:
        raise ValueError(f"embedding dimension must be >= 1, got {dim}")
    if options.pq_nbits not in {4, 8, 12, 16}:
        raise ValueError("--faiss_pq_nbits must be 4, 8, 12, or 16")
    if options.hnsw_m < 2:
        raise ValueError("--faiss_hnsw_m must be >= 2")
    ef_search = options.hnsw_ef_search if options.hnsw_ef_search is not None else DEFAULT_HNSW_EF_SEARCH
    if options.hnsw_ef_construction < 1 or ef_search < 1:
        raise ValueError("--faiss_hnsw_ef_construction and --faiss_hnsw_ef_search must be >= 1")
    if options.pq_m_refine < 1:
        raise ValueError("--faiss_pq_m_refine must be >= 1")

    metric = "inner_product"
    nlist = None
    nprobe = None
    pq_m = None
    lsh_nbits = None
    sq_type = None
    ip = faiss.METRIC_INNER_PRODUCT

    if kind == "FlatIP":
        index = faiss.IndexFlatIP(dim)
    elif kind == "FlatL2":
        index = faiss.IndexFlatL2(dim)
        metric = "l2"
    elif kind == "HNSW":
        index = faiss.IndexHNSWFlat(dim, options.hnsw_m, ip)
        index.hnsw.efConstruction = options.hnsw_ef_construction
        index.hnsw.efSearch = ef_search
    elif kind == "IVFFlat":
        nlist = resolve_nlist(n_vectors, options.nlist)
        nprobe = resolve_nprobe(nlist, options.nprobe)
        quantizer = faiss.IndexFlatIP(dim)
        index = faiss.IndexIVFFlat(quantizer, dim, nlist, ip)
        index.nprobe = nprobe
    else:
        raise ValueError(f"unhandled FAISS index type {kind}")

    details = {
        "index_type": kind,
        "index_class": type(index).__name__,
        "metric": metric,
        "gpu_capable": kind in GPU_INDEX_TYPES,
        "nlist": nlist,
        "nprobe": nprobe,
        "pq_m": pq_m,
        "pq_nbits": options.pq_nbits if pq_m is not None else None,
        "pq_m_refine": options.pq_m_refine if kind == "IVFPQR" else None,
        "hnsw_m": options.hnsw_m if kind == "HNSW" else None,
        "hnsw_ef_construction": options.hnsw_ef_construction if kind == "HNSW" else None,
        "hnsw_ef_search": ef_search if kind == "HNSW" else None,
        "lsh_nbits": lsh_nbits,
        "sq_type": sq_type,
    }
    return index, details


def index_to_gpu(index: faiss.Index, device: int = 0) -> faiss.Index:
    resources = faiss.StandardGpuResources()
    return faiss.index_cpu_to_gpu(resources, device, index)


def index_to_cpu(index: faiss.Index) -> faiss.Index:
    to_cpu = getattr(faiss, "index_gpu_to_cpu", None)
    if to_cpu is None:
        return index
    try:
        return to_cpu(index)
    except RuntimeError:
        return index


def min_train_vectors(kind: str, nlist: int | None, pq_nbits: int | None) -> int:
    needed = 1
    if kind == "IVFFlat" and nlist:
        needed = max(needed, int(nlist))
    return needed


def populate_index(index: faiss.Index, vectors: np.ndarray, details: dict[str, Any] | None = None) -> faiss.Index:
    xb = np.ascontiguousarray(vectors, dtype=np.float32)
    if xb.ndim != 2:
        raise ValueError(f"index vectors must be 2-D, got {xb.shape}")
    if details is not None:
        needed = min_train_vectors(details["index_type"], details.get("nlist"), details.get("pq_nbits"))
        if xb.shape[0] < needed:
            raise ValueError(
                f"{details['index_type']} training needs at least {needed} windows "
                f"(nlist and/or 2^pq_nbits); this database has {xb.shape[0]}. "
                "Use FlatIP/HNSW, lower --faiss_nlist / --faiss_pq_nbits, or index more sequences."
            )
    if not index.is_trained:
        index.train(xb)
    index.add(xb)
    return index


def build_populated_index(vectors: np.ndarray, options: IndexOptions) -> tuple[faiss.Index, dict[str, Any]]:
    xb = np.ascontiguousarray(vectors, dtype=np.float32)
    if xb.ndim != 2 or xb.shape[0] == 0:
        raise ValueError("need a non-empty 2-D vector matrix to build a FAISS index")
    kind = normalize_index_type(options.index_type)
    if options.gpu:
        require_gpu(kind)
    index, details = make_cpu_index(xb.shape[1], xb.shape[0], options)
    details["gpu"] = bool(options.gpu)
    if options.gpu:
        index = index_to_gpu(index, options.gpu_device)
    try:
        populate_index(index, xb, details)
    except ValueError:
        raise
    except Exception:
        if options.gpu:
            raise ValueError(
                f"failed to train/add on GPU for index type {kind}. "
                "This combination may not be implemented in GPU FAISS; "
                "retry without --faiss_gpu."
            ) from None
        raise
    if options.gpu:
        index = index_to_cpu(index)
        if details.get("nprobe"):
            apply_search_params(index, nprobe=details["nprobe"], ef_search=details.get("hnsw_ef_search"))
    return index, details


def _ivf_index(index: faiss.Index) -> Any | None:
    extractor = getattr(faiss, "extract_index_ivf", None)
    if extractor is None:
        return index if hasattr(index, "nprobe") else None
    try:
        return extractor(index)
    except RuntimeError:
        return index if hasattr(index, "nprobe") else None


def apply_search_params(
    index: faiss.Index,
    nprobe: int | None = None,
    ef_search: int | None = None,
) -> None:
    if nprobe is not None:
        ivf = _ivf_index(index)
        if ivf is not None and hasattr(ivf, "nprobe"):
            ivf.nprobe = int(nprobe)
    if ef_search is not None and hasattr(index, "hnsw"):
        index.hnsw.efSearch = int(ef_search)


def distances_to_similarity(
    distances: np.ndarray,
    metric: str,
    lsh_nbits: int | None = None,
) -> np.ndarray:
    scores = np.asarray(distances, dtype=np.float32)
    if metric == "inner_product":
        return scores
    if metric == "ngt_cosine":
        return 1.0 - scores
    if metric == "cuvs_cosine":
        return 1.0 - scores
    if metric == "l2":
        return 1.0 - (scores / np.float32(2.0))
    if metric == "hamming":
        bits = float(lsh_nbits or 1)
        return 1.0 - (scores / np.float32(bits))
    raise ValueError(f"unknown FAISS metric {metric!r}")


def meta_from_details(details: dict[str, Any]) -> dict[str, Any]:
    return {
        "index_type": details["index_type"],
        "index_class": details["index_class"],
        "metric": details["metric"],
        "gpu": bool(details.get("gpu", False)),
        "gpu_capable": bool(details.get("gpu_capable", False)),
        "nlist": details.get("nlist"),
        "nprobe": details.get("nprobe"),
        "pq_m": details.get("pq_m"),
        "pq_nbits": details.get("pq_nbits"),
        "pq_m_refine": details.get("pq_m_refine"),
        "hnsw_m": details.get("hnsw_m"),
        "hnsw_ef_construction": details.get("hnsw_ef_construction"),
        "hnsw_ef_search": details.get("hnsw_ef_search"),
        "lsh_nbits": details.get("lsh_nbits"),
        "sq_type": details.get("sq_type"),
        "backend": details.get("backend", "faiss"),
        "scann_scoring": details.get("scann_scoring"),
        "scann_partitioned": details.get("scann_partitioned"),
        "scann_reorder": details.get("scann_reorder"),
        "scann_ah_dim": details.get("scann_ah_dim"),
        "scann_anisotropic": details.get("scann_anisotropic"),
        "scann_soar": details.get("scann_soar"),
        "scann_num_neighbors": details.get("scann_num_neighbors"),
    }


def prepare_search_index(index: faiss.Index, meta: dict[str, Any], options: IndexOptions) -> tuple[faiss.Index, str, int | None]:
    kind = str(meta.get("index_type") or "FlatIP")
    metric = str(meta.get("metric") or "inner_product")
    lsh_nbits = meta.get("lsh_nbits")
    if options.gpu:
        require_gpu(kind)
        index = index_to_gpu(index, options.gpu_device)
    nprobe = options.nprobe if options.nprobe is not None else meta.get("nprobe")
    ef_search = options.hnsw_ef_search if options.hnsw_ef_search is not None else meta.get("hnsw_ef_search")
    apply_search_params(index, nprobe=nprobe, ef_search=ef_search)
    return index, metric, lsh_nbits
