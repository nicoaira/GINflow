#!/usr/bin/env python3
"""NGT, quantized graph, and quantized blob graph index support."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


NGT_INDEX_TYPES = ("NGT", "QG", "QBG")
_ALIASES = {
    "NGT": "NGT",
    "ANNG": "NGT",
    "QG": "QG",
    "QUANTIZEDGRAPH": "QG",
    "QBG": "QBG",
    "QUANTIZEDBLOBGRAPH": "QBG",
}


def normalize_index_type(name: str) -> str:
    key = name.strip().replace("-", "").replace("_", "").upper()
    if key not in _ALIASES:
        allowed = ", ".join(kind.lower() for kind in NGT_INDEX_TYPES)
        raise ValueError(f"unknown NGT index type {name!r}. Choose one of: {allowed}")
    return _ALIASES[key]


def is_ngt_type(name: str | None) -> bool:
    if not name:
        return False
    try:
        normalize_index_type(name)
    except ValueError:
        return False
    return True


def is_ngt_database(directory: Path, meta: dict[str, Any] | None = None) -> bool:
    if str((meta or {}).get("backend") or "").strip().lower() == "ngt":
        return True
    if is_ngt_type(str((meta or {}).get("index_type") or "")):
        return True
    return (directory / "ngt").is_dir()


class NgtIndex:
    """FAISS-shaped wrapper around one of the NGT Python index classes."""

    def __init__(self, path: Path, index_type: str, ntotal: int) -> None:
        try:
            import ngtpy
        except ImportError as exc:  # pragma: no cover - exercised in task images
            raise ValueError(
                "NGT was requested but ngtpy is not installed. Re-run with "
                "--index ngt so BUILD_NGT_INDEX and SEARCH_NGT use the NGT image."
            ) from exc

        self.path = Path(path)
        self.index_type = normalize_index_type(index_type)
        self.ntotal = int(ntotal)
        if self.index_type == "NGT":
            self._index = ngtpy.Index(str(self.path))
            self.metric = "ngt_cosine"
        elif self.index_type == "QG":
            self._index = ngtpy.QuantizedIndex(str(self.path), log_disabled=True)
            self.metric = "ngt_cosine"
        else:
            self._index = ngtpy.QuantizedBlobIndex(
                str(self.path),
                max_no_of_edges=128,
                log_disabled=True,
            )
            # Python's QBG two-step search otherwise starts with zero probes on
            # current NGT releases and returns no blobs for small collections.
            self._index.set(
                num_of_search_objects=20,
                result_expansion=0.0,
                exploration_size=256,
                num_of_probes=1,
            )
            self.metric = "l2"

    def search(self, queries: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
        xb = np.ascontiguousarray(queries, dtype=np.float32)
        if xb.ndim == 1:
            xb = xb.reshape(1, -1)
        if xb.ndim != 2:
            raise ValueError(f"NGT queries must be 1-D or 2-D, got {xb.shape}")
        search_k = max(1, min(int(k), self.ntotal))
        labels = np.full((xb.shape[0], search_k), -1, dtype=np.int64)
        distances = np.full((xb.shape[0], search_k), np.float32("inf"), dtype=np.float32)
        for row, query in enumerate(xb):
            results = self._index.search(query, search_k)
            for rank, (identifier, distance) in enumerate(results[:search_k]):
                labels[row, rank] = int(identifier)
                distances[row, rank] = float(distance)
        return distances, labels

    def close(self) -> None:
        index = getattr(self, "_index", None)
        if index is None:
            return
        close = getattr(index, "close", None)
        if close is not None:
            close()
        self._index = None


def _run_qbg(*args: str) -> None:
    subprocess.run(["qbg", *args], check=True)


def _build_regular(path: Path, vectors: np.ndarray) -> None:
    import ngtpy

    ngtpy.create(
        str(path),
        int(vectors.shape[1]),
        edge_size_for_creation=10,
        edge_size_for_search=40,
        distance_type="Cosine",
        object_type="Float",
    )
    index = ngtpy.Index(str(path))
    try:
        index.batch_insert(vectors)
        index.save()
    finally:
        index.close()


def build_populated_index(
    vectors: np.ndarray,
    index_type: str = "ngt",
) -> tuple[NgtIndex, dict[str, Any]]:
    xb = np.ascontiguousarray(vectors, dtype=np.float32)
    if xb.ndim != 2 or xb.shape[0] == 0:
        raise ValueError("need a non-empty 2-D vector matrix to build an NGT index")
    kind = normalize_index_type(index_type)
    root = Path(tempfile.mkdtemp(prefix="ginflow-ngt-"))
    path = root / "index"
    try:
        if kind in {"NGT", "QG"}:
            _build_regular(path, xb)
            if kind == "QG":
                _run_qbg("create-qg", str(path))
                _run_qbg("build-qg", "-o", str(xb.shape[0]), str(path))
        else:
            data = root / "vectors.tsv"
            np.savetxt(data, xb, fmt="%.9g")
            _run_qbg("create", "-d", str(xb.shape[1]), "-D", "2", str(path))
            _run_qbg("append", str(path), str(data))
            _run_qbg("build", "-o", str(xb.shape[0]), str(path))
        index = NgtIndex(path, kind, xb.shape[0])
    except Exception:
        shutil.rmtree(root, ignore_errors=True)
        raise

    details = {
        "backend": "ngt",
        "index_type": kind,
        "index_class": kind,
        "metric": index.metric,
        "gpu": False,
        "gpu_capable": False,
        "path": path,
    }
    return index, details


def load_index(path: Path, meta: dict[str, Any]) -> NgtIndex:
    kind = normalize_index_type(str(meta.get("index_type") or "ngt"))
    ntotal = meta.get("n_windows")
    if ntotal is None:
        raise ValueError("NGT database metadata is missing n_windows")
    return NgtIndex(path, kind, int(ntotal))


def serialize_index(index: NgtIndex, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(index.path, destination)
    root = index.path.parent
    index.close()
    shutil.rmtree(root, ignore_errors=True)


def meta_from_details(details: dict[str, Any]) -> dict[str, Any]:
    return {
        "backend": "ngt",
        "index_type": details["index_type"],
        "index_class": details["index_class"],
        "metric": details["metric"],
        "gpu": False,
        "gpu_capable": False,
        "ngt_index_type": details["index_type"],
    }
