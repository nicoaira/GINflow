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
NGT_OPTION_KEYS = (
    "edge_size_for_creation",
    "edge_size_for_search",
    "num_threads",
    "max_no_of_edges",
    "num_of_search_objects",
    "search_range_coefficient",
    "blob_search_range_coefficient",
    "search_radius",
    "result_expansion",
    "exploration_size",
    "exact_result_expansion",
    "num_of_probes",
    "qg_subvector_dimensions",
    "qbg_subvectors",
    "qbg_cluster_data_type",
)
_NGT_DEFAULTS = {
    "edge_size_for_creation": 10,
    "edge_size_for_search": 40,
    "num_threads": 8,
    "max_no_of_edges": None,
    "num_of_search_objects": 20,
    "search_range_coefficient": None,
    "blob_search_range_coefficient": None,
    "search_radius": None,
    "result_expansion": None,
    "exploration_size": 256,
    "exact_result_expansion": 0.0,
    "num_of_probes": 1,
    "qg_subvector_dimensions": None,
    "qbg_subvectors": None,
    "qbg_cluster_data_type": "PQ4",
}
_NGT_QBG_DEFAULT_MAX_NO_OF_EDGES = 128
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


def _normalise_options(options: dict[str, Any] | None) -> dict[str, Any]:
    unknown = set(options or ()) - set(NGT_OPTION_KEYS)
    if unknown:
        raise ValueError(f"unknown NGT option(s): {', '.join(sorted(unknown))}")
    resolved = dict(_NGT_DEFAULTS)
    resolved.update({key: value for key, value in (options or {}).items() if value is not None})
    resolved["qbg_cluster_data_type"] = str(resolved["qbg_cluster_data_type"]).upper()
    if resolved["qbg_cluster_data_type"] not in {"PQ4", "SQ8"}:
        raise ValueError(
            "NGT qbg_cluster_data_type must be pq4 or sq8, "
            f"got {resolved['qbg_cluster_data_type']}"
        )
    positive_ints = (
        "edge_size_for_creation",
        "edge_size_for_search",
        "num_threads",
        "max_no_of_edges",
        "num_of_search_objects",
        "exploration_size",
        "num_of_probes",
        "qg_subvector_dimensions",
        "qbg_subvectors",
    )
    for key in positive_ints:
        value = resolved[key]
        if value is not None and int(value) < 1:
            raise ValueError(f"NGT {key} must be >= 1, got {value}")
        if value is not None:
            resolved[key] = int(value)
    nonnegative_floats = (
        "search_range_coefficient",
        "blob_search_range_coefficient",
        "search_radius",
        "result_expansion",
        "exact_result_expansion",
    )
    for key in nonnegative_floats:
        value = resolved[key]
        if value is not None and float(value) < 0:
            raise ValueError(f"NGT {key} must be >= 0, got {value}")
        if value is not None:
            resolved[key] = float(value)
    return resolved


def _options_from_meta(meta: dict[str, Any]) -> dict[str, Any]:
    return {
        key: meta[f"ngt_{key}"]
        for key in NGT_OPTION_KEYS
        if f"ngt_{key}" in meta and meta[f"ngt_{key}"] is not None
    }


def _search_options(options: dict[str, Any]) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "num_of_search_objects": options["num_of_search_objects"],
    }
    if options["search_range_coefficient"] is not None:
        kwargs["epsilon"] = options["search_range_coefficient"]
    if options["search_radius"] is not None:
        kwargs["search_radius"] = options["search_radius"]
    return kwargs


class NgtIndex:
    """FAISS-shaped wrapper around one of the NGT Python index classes."""

    def __init__(
        self,
        path: Path,
        index_type: str,
        ntotal: int,
        options: dict[str, Any] | None = None,
    ) -> None:
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
        self.options = _normalise_options(options)
        if self.index_type == "NGT":
            self._index = ngtpy.Index(str(self.path))
            self.metric = "ngt_cosine"
        elif self.index_type == "QG":
            self._index = ngtpy.QuantizedIndex(str(self.path), log_disabled=True)
            self.metric = "ngt_cosine"
        else:
            self._index = ngtpy.QuantizedBlobIndex(
                str(self.path),
                max_no_of_edges=(
                    self.options["max_no_of_edges"] or _NGT_QBG_DEFAULT_MAX_NO_OF_EDGES
                ),
                log_disabled=True,
            )
            self.metric = "l2"
        self._configure_search()

    def _configure_search(self) -> None:
        edge_size = (
            self.options["max_no_of_edges"]
            if self.options["max_no_of_edges"] is not None
            else self.options["edge_size_for_search"]
        )
        if self.index_type == "NGT":
            kwargs = _search_options(self.options)
            kwargs["edge_size"] = edge_size
        elif self.index_type == "QG":
            kwargs = _search_options(self.options)
            if self.options["result_expansion"] is not None:
                kwargs["result_expansion"] = self.options["result_expansion"]
            kwargs["edge_size"] = edge_size
        else:
            kwargs = {
                "num_of_search_objects": self.options["num_of_search_objects"],
                "result_expansion": (
                    self.options["result_expansion"]
                    if self.options["result_expansion"] is not None
                    else 0.0
                ),
                "exploration_size": self.options["exploration_size"],
                "exact_result_expansion": self.options["exact_result_expansion"],
                "num_of_probes": self.options["num_of_probes"],
                "edge_size": edge_size,
            }
            if self.options["search_range_coefficient"] is not None:
                kwargs["epsilon"] = self.options["search_range_coefficient"]
            if self.options["blob_search_range_coefficient"] is not None:
                kwargs["blob_epsilon"] = self.options["blob_search_range_coefficient"]
            if self.options["search_radius"] is not None:
                kwargs["radius"] = self.options["search_radius"]
        self._index.set(**kwargs)

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


def _build_regular(path: Path, vectors: np.ndarray, options: dict[str, Any]) -> None:
    import ngtpy

    ngtpy.create(
        str(path),
        int(vectors.shape[1]),
        edge_size_for_creation=options["edge_size_for_creation"],
        edge_size_for_search=options["edge_size_for_search"],
        distance_type="Cosine",
        object_type="Float",
    )
    index = ngtpy.Index(str(path))
    try:
        index.batch_insert(vectors, num_threads=options["num_threads"])
        index.save()
    finally:
        index.close()


def build_populated_index(
    vectors: np.ndarray,
    index_type: str = "ngt",
    options: dict[str, Any] | None = None,
) -> tuple[NgtIndex, dict[str, Any]]:
    xb = np.ascontiguousarray(vectors, dtype=np.float32)
    if xb.ndim != 2 or xb.shape[0] == 0:
        raise ValueError("need a non-empty 2-D vector matrix to build an NGT index")
    kind = normalize_index_type(index_type)
    resolved_options = _normalise_options(options)
    root = Path(tempfile.mkdtemp(prefix="ginflow-ngt-"))
    path = root / "index"
    try:
        if kind in {"NGT", "QG"}:
            _build_regular(path, xb, resolved_options)
            if kind == "QG":
                create_args = ["create-qg"]
                if resolved_options["qg_subvector_dimensions"] is not None:
                    create_args.extend(["-Q", str(resolved_options["qg_subvector_dimensions"])])
                create_args.append(str(path))
                _run_qbg(*create_args)
                build_args = ["build-qg"]
                if resolved_options["max_no_of_edges"] is not None:
                    build_args.extend(["-E", str(resolved_options["max_no_of_edges"])])
                build_args.extend(["-o", str(xb.shape[0]), str(path)])
                _run_qbg(*build_args)
        else:
            data = root / "vectors.tsv"
            np.savetxt(data, xb, fmt="%.9g")
            create_args = ["create", "-d", str(xb.shape[1]), "-D", "2"]
            create_args.extend(["-C", resolved_options["qbg_cluster_data_type"]])
            if resolved_options["qbg_subvectors"] is not None:
                create_args.extend(["-N", str(resolved_options["qbg_subvectors"])])
            create_args.append(str(path))
            _run_qbg(*create_args)
            _run_qbg("append", str(path), str(data))
            build_args = ["build"]
            if resolved_options["max_no_of_edges"] is not None:
                build_args.extend(["-E", str(resolved_options["max_no_of_edges"])])
            build_args.extend(["-o", str(xb.shape[0]), str(path)])
            _run_qbg(*build_args)
        index = NgtIndex(path, kind, xb.shape[0], resolved_options)
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
        "ngt_options": dict(resolved_options),
    }
    return index, details


def load_index(
    path: Path,
    meta: dict[str, Any],
    options: dict[str, Any] | None = None,
) -> NgtIndex:
    kind = normalize_index_type(str(meta.get("index_type") or "ngt"))
    ntotal = meta.get("n_windows")
    if ntotal is None:
        raise ValueError("NGT database metadata is missing n_windows")
    stored_options = _options_from_meta(meta)
    stored_options.update({key: value for key, value in (options or {}).items() if value is not None})
    return NgtIndex(path, kind, int(ntotal), stored_options)


def serialize_index(index: NgtIndex, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(index.path, destination)
    root = index.path.parent
    index.close()
    shutil.rmtree(root, ignore_errors=True)


def meta_from_details(details: dict[str, Any]) -> dict[str, Any]:
    metadata = {
        "backend": "ngt",
        "index_type": details["index_type"],
        "index_class": details["index_class"],
        "metric": details["metric"],
        "gpu": False,
        "gpu_capable": False,
        "ngt_index_type": details["index_type"],
    }
    for key, value in (details.get("ngt_options") or {}).items():
        metadata[f"ngt_{key}"] = value
    return metadata
