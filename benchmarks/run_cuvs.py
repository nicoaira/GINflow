#!/usr/bin/env python3
"""Benchmark GINflow's GPU-only cuVS indexes against a prepared cache.

This runner intentionally uses :mod:`bin.cuvs_index`, rather than calling
cuVS directly.  Builds therefore exercise the same CAGRA, IVF-Flat, IVF-PQ,
cosine, serialization, and search-wrapper semantics as ``BUILD_CUVS_INDEX``
and ``SEARCH_CUVS``.  Run it in the pinned cuVS 24.10 image with a visible
NVIDIA GPU:

```
docker run --rm --gpus all -v "$PWD":/work -w /work \\
  -v /mnt/ssd_samsung/ginflow-benchmarks:/bench-cache \\
  community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4 \\
  python3 benchmarks/run_cuvs.py --cache-dir /bench-cache --dataset rouskin_6k
```

cuVS construction materializes the complete matrix on the GPU in the current
production adapter.  The runner consequently has deliberately conservative
host-RAM and free-VRAM preflight gates.  A skipped row is evidence of an
unsafe configuration, not a zero-performance measurement.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np

from benchmark_utils import (
    RUNNER_VERSION,
    PeakResourceSampler,
    hardware_snapshot,
    make_result_record,
    measure_search_repeats,
    sha256_array,
    software_snapshot,
    stable_id,
    validate_result_record,
    write_result,
)


CUVS_IMAGE = "community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4"
CUVS_VERSION_PREFIX = "24.10"
BACKEND = "cuvs"
K = 100
DEFAULT_REPEATS = 3
DEFAULT_WARMUP_QUERIES = 64
DEFAULT_QUERY_BATCH_SIZE = 32
DEFAULT_MAX_RAM_FRACTION = 0.80
DEFAULT_MAX_VRAM_FRACTION = 0.90
GIB = 1024**3


@dataclass(frozen=True)
class CacheInputs:
    """Validated backend-neutral cache inputs used by one cuVS run."""

    cache_dir: Path
    dataset_id: str
    root: Path
    vectors_path: Path
    query_path: Path
    ground_truth_path: Path
    vectors: np.ndarray
    queries: np.ndarray
    ground_truth_ids: np.ndarray
    dataset_window_count: int
    dimension: int
    embedding_cache_id: str
    query_selection_id: str
    ground_truth_cache_id: str
    query_ids_sha256: str
    ground_truth_ids_sha256: str


@dataclass(frozen=True)
class IndexSpec:
    """One GINflow-exposed cuVS build/query configuration.

    ``build_n_probes`` is deliberately separate from ``n_probes``.  IVF
    stores the former in its database metadata, while ``SEARCH_CUVS`` may
    override the latter when it reloads that database.  CAGRA does not expose
    a corresponding query-time override: its ``itopk_size`` is persisted in
    metadata and requires a rebuilt GINflow database to change.
    """

    index_type: str
    n_lists: int | None = None
    n_probes: int | None = None
    build_n_probes: int | None = None
    pq_bits: int | None = None
    pq_dim: int | None = None
    intermediate_graph_degree: int | None = None
    graph_degree: int | None = None
    build_algo: str | None = None
    itopk_size: int | None = None

    def build_key(self) -> tuple[Any, ...]:
        return (
            self.index_type,
            self.n_lists,
            self.build_n_probes,
            self.pq_bits,
            self.pq_dim,
            self.intermediate_graph_degree,
            self.graph_degree,
            self.build_algo,
            self.itopk_size,
        )

    def build_options(self) -> dict[str, Any]:
        return {
            "n_lists": self.n_lists,
            "n_probes": self.build_n_probes,
            "pq_bits": int(self.pq_bits or 8),
            "pq_dim": int(self.pq_dim or 0),
            "intermediate_graph_degree": int(self.intermediate_graph_degree or 128),
            "graph_degree": int(self.graph_degree or 64),
            "build_algo": self.build_algo or "nn_descent",
            "itopk_size": int(self.itopk_size or 64),
        }

    def parameters(self, *, query_batch_size: int) -> dict[str, Any]:
        common = {
            "cuvs_index": self.index_type.lower(),
            "cuvs_n_lists": self.n_lists,
            "cuvs_pq_bits": self.pq_bits,
            "cuvs_pq_dim": self.pq_dim,
            "cuvs_intermediate_graph_degree": self.intermediate_graph_degree,
            "cuvs_graph_degree": self.graph_degree,
            "cuvs_build_algo": self.build_algo,
            "cuvs_itopk_size": self.itopk_size,
            "query_batch_size": int(query_batch_size),
        }
        if self.index_type == "CAGRA":
            common.update(
                {
                    "cuvs_n_probes": None,
                    "cuvs_build_n_probes": None,
                    "build_persisted_controls": [
                        "cuvs_intermediate_graph_degree",
                        "cuvs_graph_degree",
                        "cuvs_build_algo",
                        "cuvs_itopk_size",
                    ],
                    "search_override_controls": [],
                }
            )
        else:
            common.update(
                {
                    "cuvs_n_probes": self.n_probes,
                    "cuvs_build_n_probes": self.build_n_probes,
                    "build_persisted_controls": [
                        "cuvs_n_lists",
                        "cuvs_pq_bits",
                        "cuvs_pq_dim",
                        "cuvs_build_n_probes",
                    ],
                    "search_override_controls": ["cuvs_n_probes"],
                }
            )
        return common

    def label(self) -> str:
        pieces = [self.index_type.lower()]
        if self.n_lists is not None:
            pieces.append(f"nlist{self.n_lists}")
        if self.pq_bits is not None:
            pieces.append(f"bits{self.pq_bits}")
        if self.pq_dim is not None:
            pieces.append(f"pqdim{self.pq_dim}")
        if self.intermediate_graph_degree is not None:
            pieces.append(f"intermediate{self.intermediate_graph_degree}")
        if self.graph_degree is not None:
            pieces.append(f"degree{self.graph_degree}")
        if self.itopk_size is not None:
            pieces.append(f"itopk{self.itopk_size}")
        if self.n_probes is not None:
            pieces.append(f"nprobe{self.n_probes}")
        return "-".join(pieces)


@dataclass(frozen=True)
class CapacityCheck:
    """Conservative adapter-aware capacity evidence for one configuration."""

    feasible: bool
    estimated_vram_bytes: int
    available_vram_bytes: int | None
    vram_safe_limit_bytes: int | None
    estimated_ram_bytes: int
    available_ram_bytes: int | None
    ram_safe_limit_bytes: int | None
    raw_vector_bytes: int
    reason: str | None = None


@dataclass
class BuildArtifact:
    """Serialized production index and build-time measurements."""

    directory: Path
    meta: dict[str, Any]
    details: dict[str, Any]
    build_seconds: float
    serialize_seconds: float
    index_bytes: int
    peak_rss_bytes: int | None
    peak_vram_bytes: int | None


def cache_root_from_args(value: str | None) -> Path:
    """Resolve an external cache without making an implicit local corpus."""

    configured = value or os.environ.get("GINFLOW_BENCHMARK_CACHE")
    return Path(configured).expanduser().resolve() if configured else (Path.cwd() / ".benchmark-cache").resolve()


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _required_text(mapping: Mapping[str, Any], key: str, source: Path) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{source} is missing a non-empty {key!r}")
    return value


def _artifact_path(directory: Path, configured: Any, fallback: str) -> Path:
    value = str(configured or fallback)
    candidate = (directory / value).resolve()
    root = directory.resolve()
    if candidate != root and root not in candidate.parents:
        raise ValueError(f"cache manifest path escapes {directory}: {value}")
    return candidate


def load_cache_inputs(cache_dir: Path, dataset_id: str) -> CacheInputs:
    """Load cache artifacts and reject incompatible query/ground-truth data."""

    root = cache_dir / dataset_id
    flat_dir = root / "flat"
    query_dir = root / "queries"
    truth_dir = root / "ground-truth"
    flatten_path = flat_dir / "flatten-manifest.json"
    selection_path = query_dir / "query-selection.json"
    truth_path = truth_dir / "ground-truth.json"
    flatten = _read_json(flatten_path)
    selection = _read_json(selection_path)
    truth = _read_json(truth_path)
    vectors_path = _artifact_path(flat_dir, flatten.get("vectors", {}).get("path"), "vectors.npy")
    query_path = _artifact_path(query_dir, selection.get("query_vectors", {}).get("path"), "queries.npy")
    ground_truth_path = _artifact_path(
        truth_dir, truth.get("ground_truth", {}).get("path"), "ground-truth.npz"
    )
    for path in (vectors_path, query_path, ground_truth_path):
        if not path.is_file():
            raise ValueError(f"prepared cache artifact is missing: {path}")

    vectors = np.load(vectors_path, mmap_mode="r", allow_pickle=False)
    queries = np.load(query_path, allow_pickle=False)
    with np.load(ground_truth_path, allow_pickle=False) as archive:
        if "ids" not in archive:
            raise ValueError(f"{ground_truth_path} lacks exact neighbour ids")
        ground_truth_ids = np.asarray(archive["ids"], dtype=np.int64)
    if vectors.ndim != 2 or vectors.shape[0] < K or vectors.dtype != np.float32:
        raise ValueError(f"{vectors_path} must be a float32 2-D database with at least {K} rows")
    if queries.ndim != 2 or queries.shape[0] < 1 or queries.shape[1] != vectors.shape[1]:
        raise ValueError("prepared query vectors do not match the flattened vector dimension")
    if ground_truth_ids.shape != (queries.shape[0], K):
        raise ValueError(
            f"exact ground truth shape {ground_truth_ids.shape} must be ({queries.shape[0]}, {K})"
        )
    if np.any(ground_truth_ids < 0) or np.any(ground_truth_ids >= vectors.shape[0]):
        raise ValueError("exact ground truth contains invalid vector identifiers")
    sample_rows = min(vectors.shape[0], 1_024)
    if not np.isfinite(queries).all() or not np.isfinite(vectors[:sample_rows]).all():
        raise ValueError("prepared cache contains non-finite vectors")

    embedding_cache_id = _required_text(flatten, "embedding_cache_id", flatten_path)
    query_selection_id = _required_text(selection, "query_selection_id", selection_path)
    ground_truth_cache_id = _required_text(truth, "ground_truth_cache_id", truth_path)
    if selection.get("embedding_cache_id") != embedding_cache_id:
        raise ValueError("query selection and flattened vectors have different embedding_cache_id values")
    if truth.get("embedding_cache_id") != embedding_cache_id:
        raise ValueError("exact ground truth and flattened vectors have different embedding_cache_id values")
    if truth.get("query_selection_id") != query_selection_id:
        raise ValueError("exact ground truth and query selection have different query_selection_id values")
    if truth.get("metric") != "cosine" or int(truth.get("k") or 0) != K:
        raise ValueError("exact ground truth must be cosine top-100")
    if int(truth.get("database_window_count") or 0) != vectors.shape[0]:
        raise ValueError("exact ground truth database count does not match flattened vectors")
    if int(truth.get("dimension") or 0) != vectors.shape[1]:
        raise ValueError("exact ground truth dimension does not match flattened vectors")
    query_ids_sha256 = _required_text(selection, "query_ids_sha256", selection_path)
    if truth.get("query_ids_sha256") != query_ids_sha256:
        raise ValueError("exact ground truth and query selection have different query id digests")
    ground_truth_ids_sha256 = _required_text(truth.get("ground_truth", {}), "ids_sha256", truth_path)
    if sha256_array(ground_truth_ids) != ground_truth_ids_sha256:
        raise ValueError("exact ground truth id digest does not match its manifest")
    if any(str(item.get("dataset_id")) != dataset_id for item in (flatten, selection, truth)):
        raise ValueError("prepared cache manifests do not match the requested dataset id")
    return CacheInputs(
        cache_dir=cache_dir,
        dataset_id=dataset_id,
        root=root,
        vectors_path=vectors_path,
        query_path=query_path,
        ground_truth_path=ground_truth_path,
        vectors=vectors,
        queries=np.ascontiguousarray(queries, dtype=np.float32),
        ground_truth_ids=ground_truth_ids,
        dataset_window_count=int(vectors.shape[0]),
        dimension=int(vectors.shape[1]),
        embedding_cache_id=embedding_cache_id,
        query_selection_id=query_selection_id,
        ground_truth_cache_id=ground_truth_cache_id,
        query_ids_sha256=query_ids_sha256,
        ground_truth_ids_sha256=ground_truth_ids_sha256,
    )


def _unique_positive(values: Iterable[int], upper: int) -> list[int]:
    return sorted({max(1, min(int(value), upper)) for value in values})


def _default_nlists(n_vectors: int, *, plan: str) -> list[int]:
    if plan == "smoke":
        return _unique_positive((min(32, n_vectors),), n_vectors)
    values = (2_048, 4_096) if n_vectors < 1_500_000 else (4_096, 8_192)
    return _unique_positive(values, n_vectors)


def _probe_ladder(n_lists: int, *, plan: str) -> list[int]:
    values = (1, min(4, n_lists)) if plan == "smoke" else (8, 16, 32, 64, 128)
    return _unique_positive(values, n_lists)


def _dedupe_specs(specs: Iterable[IndexSpec]) -> list[IndexSpec]:
    unique: list[IndexSpec] = []
    seen: set[tuple[Any, ...]] = set()
    for spec in specs:
        key = tuple(asdict(spec).items())
        if key not in seen:
            seen.add(key)
            unique.append(spec)
    return unique


def frontier_specs(
    n_vectors: int,
    dimension: int,
    *,
    plan: str,
    cagra_build_algo: str = "nn_descent",
    ivf_nlists: Sequence[int] | None = None,
    pq_bits: Sequence[int] = (4, 6, 8),
    pq_dims: Sequence[int] = (0, 256),
) -> list[IndexSpec]:
    """Return the exposed cuVS recall/QPS frontier without hidden knobs.

    The frontier varies each expensive CAGRA control around one baseline rather
    than making a 3x3x3 cartesian product.  IVF-PQ tests the full bit ladder
    at the first list count, a list-count sensitivity at 8 bits, and one
    explicit ``pq_dim`` sensitivity.  This keeps the overnight plan finite
    while covering every pipeline-exposed control.
    """

    if n_vectors < K:
        raise ValueError(f"cuVS benchmark requires at least k={K} vectors")
    if plan not in {"frontier", "smoke"}:
        raise ValueError("plan must be frontier or smoke")
    if cagra_build_algo not in {"ivf_pq", "nn_descent", "iterative_cagra_search", "ace"}:
        raise ValueError("unknown CAGRA build algorithm")
    list_counts = _unique_positive(ivf_nlists or _default_nlists(n_vectors, plan=plan), n_vectors)
    if not list_counts:
        raise ValueError("at least one IVF list count is required")
    bits = sorted({int(value) for value in pq_bits})
    if any(value not in {4, 5, 6, 7, 8} for value in bits):
        raise ValueError("--pq-bits values must be in 4,5,6,7,8")
    dimensions = sorted({int(value) for value in pq_dims})
    if any(value < 0 or value > dimension for value in dimensions):
        raise ValueError(f"--pq-dims values must be between 0 and {dimension}")
    if 0 not in dimensions:
        dimensions.insert(0, 0)

    specs: list[IndexSpec] = []
    if plan == "smoke":
        specs.append(
            IndexSpec(
                "CAGRA",
                intermediate_graph_degree=64,
                graph_degree=32,
                build_algo=cagra_build_algo,
                itopk_size=100,
            )
        )
    else:
        # ``itopk`` is a persisted CAGRA search setting in GINflow.  It needs
        # a separate build even though cuVS itself constructs it at load time.
        specs.extend(
            (
                IndexSpec("CAGRA", intermediate_graph_degree=128, graph_degree=64, build_algo=cagra_build_algo, itopk_size=100),
                IndexSpec("CAGRA", intermediate_graph_degree=128, graph_degree=64, build_algo=cagra_build_algo, itopk_size=128),
                IndexSpec("CAGRA", intermediate_graph_degree=128, graph_degree=64, build_algo=cagra_build_algo, itopk_size=256),
                IndexSpec("CAGRA", intermediate_graph_degree=256, graph_degree=64, build_algo=cagra_build_algo, itopk_size=128),
                IndexSpec("CAGRA", intermediate_graph_degree=256, graph_degree=128, build_algo=cagra_build_algo, itopk_size=128),
            )
        )

    for n_lists in list_counts:
        probes = _probe_ladder(n_lists, plan=plan)
        build_probe = probes[0]
        specs.extend(
            IndexSpec("IVF", n_lists=n_lists, build_n_probes=build_probe, n_probes=nprobe)
            for nprobe in probes
        )

    primary_nlist = list_counts[0]
    primary_probes = _probe_ladder(primary_nlist, plan=plan)
    primary_bits = (8,) if plan == "smoke" else tuple(bits)
    for bits_value in primary_bits:
        specs.extend(
            IndexSpec(
                "IVF-PQ",
                n_lists=primary_nlist,
                build_n_probes=primary_probes[0],
                n_probes=nprobe,
                pq_bits=bits_value,
                pq_dim=0,
            )
            for nprobe in primary_probes
        )
    # A second list count is a build sensitivity, not a full duplicated grid.
    if plan == "frontier" and len(list_counts) > 1:
        n_lists = list_counts[1]
        probes = _probe_ladder(n_lists, plan=plan)
        specs.extend(
            IndexSpec("IVF-PQ", n_lists=n_lists, build_n_probes=probes[0], n_probes=nprobe, pq_bits=8, pq_dim=0)
            for nprobe in probes
        )
    # Explicit pq_dim is adapter-visible but version-sensitive.  Include one
    # 8-bit sensitivity; a cuVS 24.10 rejection becomes an auditable error.
    explicit_dims = [value for value in dimensions if value != 0]
    if plan == "frontier" and explicit_dims:
        pq_dim = explicit_dims[0]
        specs.extend(
            IndexSpec(
                "IVF-PQ",
                n_lists=primary_nlist,
                build_n_probes=primary_probes[0],
                n_probes=nprobe,
                pq_bits=8,
                pq_dim=pq_dim,
            )
            for nprobe in primary_probes
        )
    return _dedupe_specs(specs)


def group_build_specs(specs: Iterable[IndexSpec]) -> Iterator[tuple[IndexSpec, list[IndexSpec]]]:
    groups: dict[tuple[Any, ...], list[IndexSpec]] = {}
    for spec in specs:
        groups.setdefault(spec.build_key(), []).append(spec)
    for values in groups.values():
        yield values[0], values


def _available_ram_bytes() -> int | None:
    meminfo = Path("/proc/meminfo")
    if not meminfo.is_file():
        return None
    match = re.search(r"^MemAvailable:\s+(\d+)\s+kB", meminfo.read_text(errors="replace"), re.MULTILINE)
    return int(match.group(1)) * 1024 if match else None


def _available_vram_bytes(device: int) -> int | None:
    try:
        output = subprocess.check_output(
            [
                "nvidia-smi",
                "--query-gpu=memory.total,memory.used",
                "--format=csv,noheader,nounits",
            ],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=10,
        )
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    rows = [line.strip() for line in output.splitlines() if line.strip()]
    if not 0 <= device < len(rows):
        return None
    try:
        total_mib, used_mib = (float(value.strip()) for value in rows[device].split(",", 1))
    except (IndexError, ValueError):
        return None
    return max(0, int((total_mib - used_mib) * 1024 * 1024))


def capacity_check(
    spec: IndexSpec,
    *,
    n_vectors: int,
    dimension: int,
    available_vram_bytes: int | None,
    available_ram_bytes: int | None,
    max_vram_fraction: float,
    max_ram_fraction: float,
) -> CapacityCheck:
    """Refuse configurations the current full-matrix adapter cannot safely hold.

    ``build_populated_index`` makes a contiguous NumPy copy from the memmap,
    uploads the complete copy with ``cp.asarray``, and cuVS then builds its
    native structure.  The estimates intentionally retain both adapter copies
    and a structure/workspace allowance.  They are a safety boundary, not a
    claimed cuVS memory formula.
    """

    raw = int(n_vectors) * int(dimension) * np.dtype(np.float32).itemsize
    host_estimate = (2 * raw) + (512 * 1024**2)
    if spec.index_type == "CAGRA":
        intermediate = int(spec.intermediate_graph_degree or 128)
        graph = int(spec.graph_degree or 64)
        edge_bytes = n_vectors * (intermediate + graph) * np.dtype(np.int32).itemsize
        # NN-descent and graph construction need a non-trivial candidate/work
        # buffer in addition to the final adjacency.  This is intentionally
        # modest enough to let a baseline 6k pilot proceed on an idle 8 GiB GPU.
        workspace = max(int(raw * 0.35), 512 * 1024**2)
        vram_estimate = raw + edge_bytes + workspace + (384 * 1024**2)
    elif spec.index_type == "IVF":
        centroids = int(spec.n_lists or 1) * dimension * np.dtype(np.float32).itemsize
        vram_estimate = raw + centroids + max(int(raw * 0.18), 384 * 1024**2) + (256 * 1024**2)
    elif spec.index_type == "IVF-PQ":
        centroids = int(spec.n_lists or 1) * dimension * np.dtype(np.float32).itemsize
        codebook = (2 ** int(spec.pq_bits or 8)) * dimension * np.dtype(np.float32).itemsize
        vram_estimate = raw + centroids + codebook + max(int(raw * 0.28), 512 * 1024**2) + (256 * 1024**2)
    else:  # pragma: no cover - IndexSpec construction is local
        raise ValueError(f"unknown cuVS index type {spec.index_type!r}")

    ram_limit = int(available_ram_bytes * max_ram_fraction) if available_ram_bytes is not None else None
    vram_limit = int(available_vram_bytes * max_vram_fraction) if available_vram_bytes is not None else None
    reason: str | None = None
    if available_vram_bytes is None:
        reason = "no visible NVIDIA free-VRAM information; refusing an unsafe GPU-only cuVS build"
    elif raw > vram_limit:
        reason = (
            "the production cuVS adapter must upload the complete raw float32 matrix "
            f"({raw} bytes), which exceeds the configured free-VRAM limit ({vram_limit} bytes)"
        )
    elif vram_estimate > vram_limit:
        reason = (
            "adapter-aware cuVS build estimate exceeds the configured free-VRAM limit "
            f"({vram_estimate} > {vram_limit} bytes)"
        )
    elif ram_limit is not None and host_estimate > ram_limit:
        reason = (
            "adapter-aware cuVS host materialization estimate exceeds the configured MemAvailable limit "
            f"({host_estimate} > {ram_limit} bytes)"
        )
    return CapacityCheck(
        feasible=reason is None,
        estimated_vram_bytes=int(vram_estimate),
        available_vram_bytes=available_vram_bytes,
        vram_safe_limit_bytes=vram_limit,
        estimated_ram_bytes=int(host_estimate),
        available_ram_bytes=available_ram_bytes,
        ram_safe_limit_bytes=ram_limit,
        raw_vector_bytes=raw,
        reason=reason,
    )


def _pipeline_cuvs_adapter() -> tuple[Any, Any]:
    """Load the pinned cuVS runtime and GINflow's production adapter."""

    try:
        installed = package_version("cuvs")
    except PackageNotFoundError as exc:  # pragma: no cover - pinned-image check
        raise RuntimeError("cuVS is not installed; run this runner in the pinned cuVS 24.10 container") from exc
    if not installed.startswith(f"{CUVS_VERSION_PREFIX}."):
        raise RuntimeError(f"this benchmark requires cuVS {CUVS_VERSION_PREFIX}.x, found {installed}")
    try:
        import cupy as cp
    except ImportError as exc:  # pragma: no cover - pinned-image check
        raise RuntimeError("CuPy is not installed; run this runner in the pinned cuVS container") from exc
    bin_dir = Path(__file__).resolve().parents[1] / "bin"
    if str(bin_dir) not in sys.path:
        sys.path.insert(0, str(bin_dir))
    adapter = importlib.import_module("cuvs_index")
    adapter.require_gpu()
    return adapter, cp


def _directory_bytes(directory: Path) -> int:
    return sum(path.stat().st_size for path in directory.rglob("*") if path.is_file())


def _safe_remove_scratch(path: Path, scratch_root: Path) -> None:
    """Remove only a fresh temporary directory created by this runner."""

    if not path.exists():
        return
    if path.parent.resolve() != scratch_root.resolve() or not path.name.startswith("cuvs-"):
        raise RuntimeError(f"refusing to remove non-runner cuVS scratch path: {path}")
    shutil.rmtree(path)


def _empty_cupy_pools(cp: Any) -> None:
    """Release cached allocations between separate large GPU construction runs."""

    try:
        cp.cuda.Stream.null.synchronize()
        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    except (AttributeError, RuntimeError):
        pass


def _cupy_free_vram_bytes(cp: Any, gpu_device: int) -> int | None:
    """Use the CUDA runtime when an image does not expose ``nvidia-smi``."""

    try:
        cp.cuda.Device(gpu_device).use()
        free_bytes, _ = cp.cuda.runtime.memGetInfo()
        return max(0, int(free_bytes))
    except (AttributeError, RuntimeError):
        return None


def build_index(
    adapter: Any,
    cp: Any,
    cache: CacheInputs,
    spec: IndexSpec,
    *,
    gpu_device: int,
    scratch_root: Path,
) -> BuildArtifact:
    """Build, serialize, and close a production cuVS index once per build key."""

    scratch_root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="cuvs-", dir=scratch_root))
    index: Any | None = None
    cp.cuda.Device(gpu_device).use()
    try:
        with PeakResourceSampler(gpu_device=gpu_device) as sampler:
            started = time.perf_counter()
            index, details = adapter.build_populated_index(
                cache.vectors,
                spec.index_type,
                **spec.build_options(),
            )
            cp.cuda.Device(gpu_device).synchronize()
            build_seconds = time.perf_counter() - started
            if int(index.ntotal) != cache.dataset_window_count:
                raise RuntimeError(
                    f"cuVS built {index.ntotal} vectors, expected {cache.dataset_window_count}"
                )
            serialize_started = time.perf_counter()
            adapter.serialize_index(index, directory)
            index = None
            cp.cuda.Device(gpu_device).synchronize()
            serialize_seconds = time.perf_counter() - serialize_started
        meta = adapter.meta_from_details(dict(details))
        meta["n_windows"] = cache.dataset_window_count
        return BuildArtifact(
            directory=directory,
            meta=meta,
            details=dict(details),
            build_seconds=float(build_seconds),
            serialize_seconds=float(serialize_seconds),
            index_bytes=_directory_bytes(directory),
            peak_rss_bytes=sampler.peak_rss_bytes,
            peak_vram_bytes=sampler.peak_vram_bytes,
        )
    except Exception:
        if index is not None:
            index.close()
        _empty_cupy_pools(cp)
        _safe_remove_scratch(directory, scratch_root)
        raise


def _load_for_search(adapter: Any, cp: Any, artifact: BuildArtifact, spec: IndexSpec, gpu_device: int) -> tuple[Any, float]:
    cp.cuda.Device(gpu_device).use()
    started = time.perf_counter()
    index = adapter.load_index(artifact.directory, artifact.meta, n_probes=spec.n_probes)
    cp.cuda.Device(gpu_device).synchronize()
    return index, float(time.perf_counter() - started)


def _run_id(cache: CacheInputs, spec: IndexSpec, parameters: Mapping[str, Any], hardware_id: str, image: str) -> str:
    return stable_id(
        "cuvs-run",
        {
            "dataset_id": cache.dataset_id,
            "embedding_cache_id": cache.embedding_cache_id,
            "query_selection_id": cache.query_selection_id,
            "ground_truth_cache_id": cache.ground_truth_cache_id,
            "hardware_id": hardware_id,
            "image": image,
            "spec": asdict(spec),
            "parameters": dict(parameters),
        },
    )


def _result_path(output_dir: Path, label: str, run_id: str, repeat: int, *, status: str = "ok") -> Path:
    slug = re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")
    suffix = f"-repeat-{repeat}" if status == "ok" else f"-repeat-{repeat}-{status}"
    return output_dir / f"{slug}-{run_id}{suffix}.json"


def _is_complete_result(path: Path, run_id: str) -> bool:
    if not path.is_file():
        return False
    try:
        record = _read_json(path)
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    return record.get("status") == "ok" and record.get("run_id") == run_id and not validate_result_record(record)


def _cuvs_runtime_metadata(cp: Any) -> dict[str, Any]:
    runtime = None
    try:
        runtime = int(cp.cuda.runtime.runtimeGetVersion())
    except RuntimeError:
        pass
    try:
        cuvs_version = package_version("cuvs")
    except PackageNotFoundError:  # pragma: no cover - checked before this call
        cuvs_version = "unknown"
    return {"cuvs_version": cuvs_version, "cupy_version": str(cp.__version__), "cuda_runtime_version": runtime}


def _base_provenance(
    cache: CacheInputs,
    *,
    repository_root: Path,
    hardware: Mapping[str, Any],
    image: str,
    container: str,
    gpu_device: int,
    runtime: Mapping[str, Any] | None,
) -> dict[str, Any]:
    software = software_snapshot(repository_root)
    provenance: dict[str, Any] = {
        "git_commit": software.get("git_commit") or "unknown",
        "runner_version": RUNNER_VERSION,
        "hardware_id": str(hardware["id"]),
        "embedding_cache_id": cache.embedding_cache_id,
        "query_selection_id": cache.query_selection_id,
        "ground_truth_cache_id": cache.ground_truth_cache_id,
        "container": container,
        "image": image,
        "hardware": dict(hardware),
        "software": software,
        "gpu_device": int(gpu_device),
        "index_constructor": "bin/cuvs_index.py::build_populated_index",
        "index_serializer": "bin/cuvs_index.py::serialize_index",
        "index_loader": "bin/cuvs_index.py::load_index",
        "search_adapter": "bin/cuvs_index.py::CuvsIndex.search",
        "benchmark_build_method": "full-memmap-through-production-cuvs-adapter-v1",
        "production_build_difference": (
            "The benchmark passes flat/vectors.npy as a float32 memmap to the production adapter; "
            "BUILD_CUVS_INDEX first concatenates window shard arrays. Both paths then invoke the "
            "adapter's full contiguous host materialization and GPU upload. Search, serialization, "
            "and reload use production code; benchmark build RSS is not a direct pipeline peak-RSS requirement."
        ),
        "vram_measurement": "GPU-wide device usage sampled with nvidia-smi; not per-process attribution",
        "source_references": [
            "https://github.com/NVIDIA/cuvs/tree/branch-24.10",
            "https://arxiv.org/abs/2308.15136",
        ],
    }
    if runtime is not None:
        provenance.update(dict(runtime))
    return provenance


def _unavailable_record(
    cache: CacheInputs,
    spec: IndexSpec,
    *,
    parameters: Mapping[str, Any],
    run_id: str,
    provenance: Mapping[str, Any],
    status: str,
    reason: str,
    capacity: CapacityCheck | None,
) -> dict[str, Any]:
    measurement: dict[str, Any] = {"timing_scope": "not measured"}
    if capacity is not None:
        measurement["capacity_preflight"] = {
            **asdict(capacity),
            "method": (
                "production adapter full host contiguous copy plus full GPU upload, "
                "native structure estimate, and construction workspace reserve"
            ),
        }
    return make_result_record(
        backend=BACKEND,
        dataset_id=cache.dataset_id,
        dataset_window_count=cache.dataset_window_count,
        dimension=cache.dimension,
        parameter_label=spec.label(),
        parameters=parameters,
        run_id=run_id,
        repeat=0,
        warmup_queries=0,
        timed_queries=0,
        query_ids_sha256=cache.query_ids_sha256,
        ground_truth_ids_sha256=cache.ground_truth_ids_sha256,
        provenance=provenance,
        status=status,
        error=reason,
        measurement=measurement,
    )


def _measurement_payload(measured: Any, artifact: BuildArtifact, load_seconds: float, spec: IndexSpec) -> dict[str, Any]:
    payload = dict(getattr(measured, "measurement", {}) or {})
    payload.update(
        {
            "timing_scope": "warm serialized-and-reloaded cuVS search; build, serialization, and load excluded from QPS",
            "index_constructor": "bin/cuvs_index.py::build_populated_index",
            "index_search": "bin/cuvs_index.py::CuvsIndex.search",
            "build_peak_rss_bytes": artifact.peak_rss_bytes,
            "build_peak_vram_bytes": artifact.peak_vram_bytes,
            "build_seconds_scope": "build_populated_index including full host materialization and GPU upload",
            "serialize_seconds": artifact.serialize_seconds,
            "load_seconds": load_seconds,
            "serialized_index_bytes_method": "sum of files from bin/cuvs_index.py::serialize_index",
            "build_persisted_controls": spec.parameters(query_batch_size=payload["query_batch_size"])["build_persisted_controls"],
            "search_override_controls": spec.parameters(query_batch_size=payload["query_batch_size"])["search_override_controls"],
            "production_build_rss_comparable": False,
        }
    )
    return payload


def _only_filter(specs: Sequence[IndexSpec], requested: str | None) -> list[IndexSpec]:
    if not requested:
        return list(specs)
    aliases = {"cagra": "CAGRA", "ivf": "IVF", "ivf-pq": "IVF-PQ", "ivfpq": "IVF-PQ"}
    selected = {item.strip().lower() for item in requested.split(",") if item.strip()}
    unknown = sorted(selected.difference(aliases))
    if unknown:
        raise ValueError(f"unknown --only values: {', '.join(unknown)}")
    kinds = {aliases[item] for item in selected}
    return [spec for spec in specs if spec.index_type in kinds]


def _parse_int_list(value: str, *, option: str) -> list[int]:
    try:
        values = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise ValueError(f"{option} must be comma-separated integers") from exc
    if not values:
        raise ValueError(f"{option} must contain at least one integer")
    return values


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    parser.add_argument("--dataset", required=True, help="prepared cache dataset id")
    parser.add_argument("--output-dir", type=Path, help="defaults to <cache>/<dataset>/results/cuvs")
    parser.add_argument("--gpu-device", type=int, default=0)
    parser.add_argument("--plan", choices=("frontier", "smoke"), default="frontier")
    parser.add_argument("--only", help="comma-separated types: cagra,ivf,ivf-pq")
    parser.add_argument("--cagra-build-algo", default="nn_descent", choices=("ivf_pq", "nn_descent", "iterative_cagra_search", "ace"))
    parser.add_argument("--ivf-nlists", help="optional comma-separated IVF build list counts")
    parser.add_argument("--pq-bits", default="4,6,8", help="comma-separated IVF-PQ bit widths")
    parser.add_argument(
        "--pq-dims",
        help=(
            "comma-separated IVF-PQ dimensions; 0 uses the cuVS heuristic. "
            "Defaults to 0,256 when the embedding dimension permits 256, otherwise 0."
        ),
    )
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--warmup-queries", type=int, default=DEFAULT_WARMUP_QUERIES)
    parser.add_argument("--query-batch-size", type=int, default=DEFAULT_QUERY_BATCH_SIZE)
    parser.add_argument("--max-ram-fraction", type=float, default=DEFAULT_MAX_RAM_FRACTION)
    parser.add_argument("--max-vram-fraction", type=float, default=DEFAULT_MAX_VRAM_FRACTION)
    parser.add_argument("--image", default=CUVS_IMAGE, help="pinned runtime image recorded in provenance")
    parser.add_argument("--container", default="docker", help="runtime type recorded in provenance")
    parser.add_argument("--keep-index-artifacts", action="store_true", help="retain runner-owned serialized indexes under cache/scratch")
    parser.add_argument("--no-resume", action="store_true", help="re-measure settings with complete result records")
    parser.add_argument("--dry-run", action="store_true", help="validate cache and print capacity plan without importing cuVS")
    return parser.parse_args(argv)


def _validate_args(args: argparse.Namespace, cache: CacheInputs) -> None:
    for name in ("repeats", "warmup_queries", "query_batch_size"):
        if int(getattr(args, name)) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 1")
    for name in ("max_ram_fraction", "max_vram_fraction"):
        value = float(getattr(args, name))
        if not 0.0 < value <= 1.0:
            raise ValueError(f"--{name.replace('_', '-')} must be in (0, 1]")
    if args.gpu_device < 0:
        raise ValueError("--gpu-device must be >= 0")
    if cache.queries.shape[0] % args.query_batch_size:
        raise ValueError("prepared query count must be an exact multiple of --query-batch-size")
    if args.warmup_queries > cache.queries.shape[0]:
        raise ValueError("--warmup-queries must not exceed the prepared query count")
    if args.warmup_queries % args.query_batch_size:
        raise ValueError("--warmup-queries must be an exact multiple of --query-batch-size")
    if args.plan == "frontier" and args.repeats < DEFAULT_REPEATS:
        raise ValueError("frontier runs require at least three timed repeats; use --plan smoke for an integration check")
    if args.dataset in {"rouskin_6k", "rouskin_30k"} and args.repeats < DEFAULT_REPEATS:
        raise ValueError("rouskin benchmark datasets require --repeats >= 3")


def _plan_payload(
    cache: CacheInputs,
    specs: Sequence[IndexSpec],
    capacities: Mapping[IndexSpec, CapacityCheck],
    *,
    hardware: Mapping[str, Any],
    image: str,
    args: argparse.Namespace,
) -> dict[str, Any]:
    return {
        "schema_version": "ginflow-cuvs-benchmark-plan-v1",
        "dataset_id": cache.dataset_id,
        "dataset_window_count": cache.dataset_window_count,
        "dimension": cache.dimension,
        "image": image,
        "hardware_id": hardware["id"],
        "gpu_device": args.gpu_device,
        "plan": args.plan,
        "settings": [
            {
                "label": spec.label(),
                "spec": asdict(spec),
                "parameters": spec.parameters(query_batch_size=args.query_batch_size),
                "capacity": asdict(capacities[spec]),
            }
            for spec in specs
        ],
    }


def run(args: argparse.Namespace) -> int:
    repository_root = Path(__file__).resolve().parents[1]
    cache = load_cache_inputs(cache_root_from_args(args.cache_dir), args.dataset)
    _validate_args(args, cache)
    ivf_nlists = _parse_int_list(args.ivf_nlists, option="--ivf-nlists") if args.ivf_nlists else None
    pq_dims = (
        _parse_int_list(args.pq_dims, option="--pq-dims")
        if args.pq_dims
        else ([0, 256] if cache.dimension >= 256 else [0])
    )
    specs = _only_filter(
        frontier_specs(
            cache.dataset_window_count,
            cache.dimension,
            plan=args.plan,
            cagra_build_algo=args.cagra_build_algo,
            ivf_nlists=ivf_nlists,
            pq_bits=_parse_int_list(args.pq_bits, option="--pq-bits"),
            pq_dims=pq_dims,
        ),
        args.only,
    )
    if not specs:
        raise ValueError("the selected --only filter produced no benchmark settings")
    hardware = hardware_snapshot()
    available_ram = _available_ram_bytes()
    available_vram = _available_vram_bytes(args.gpu_device)
    capacities = {
        spec: capacity_check(
            spec,
            n_vectors=cache.dataset_window_count,
            dimension=cache.dimension,
            available_vram_bytes=available_vram,
            available_ram_bytes=available_ram,
            max_vram_fraction=args.max_vram_fraction,
            max_ram_fraction=args.max_ram_fraction,
        )
        for spec in specs
    }
    if args.dry_run:
        print(json.dumps(_plan_payload(cache, specs, capacities, hardware=hardware, image=args.image, args=args), indent=2, sort_keys=True))
        return 0

    adapter, cp = _pipeline_cuvs_adapter()
    if int(cp.cuda.runtime.getDeviceCount()) <= args.gpu_device:
        raise RuntimeError(f"requested CUDA device {args.gpu_device} is not visible to cuVS")
    cp.cuda.Device(args.gpu_device).use()
    if available_vram is None:
        # Keep dry-run dependency-free, but do not turn a container without
        # nvidia-smi into a false blanket skip when CUDA can report free VRAM.
        available_vram = _cupy_free_vram_bytes(cp, args.gpu_device)
        capacities = {
            spec: capacity_check(
                spec,
                n_vectors=cache.dataset_window_count,
                dimension=cache.dimension,
                available_vram_bytes=available_vram,
                available_ram_bytes=available_ram,
                max_vram_fraction=args.max_vram_fraction,
                max_ram_fraction=args.max_ram_fraction,
            )
            for spec in specs
        }
    output_dir = (args.output_dir or (cache.root / "results" / BACKEND)).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = _base_provenance(
        cache,
        repository_root=repository_root,
        hardware=hardware,
        image=args.image,
        container=args.container,
        gpu_device=args.gpu_device,
        runtime=_cuvs_runtime_metadata(cp),
    )
    scratch_root = cache.root / "scratch" / BACKEND

    try:
        for build_spec, search_specs in group_build_specs(specs):
            active: list[tuple[IndexSpec, dict[str, Any], str, list[Path]]] = []
            for spec in search_specs:
                parameters = spec.parameters(query_batch_size=args.query_batch_size)
                run_id = _run_id(cache, spec, parameters, str(hardware["id"]), args.image)
                paths = [_result_path(output_dir, spec.label(), run_id, repeat) for repeat in range(args.repeats)]
                if not args.no_resume and all(_is_complete_result(path, run_id) for path in paths):
                    continue
                capacity = capacities[spec]
                if not capacity.feasible:
                    record = _unavailable_record(
                        cache,
                        spec,
                        parameters=parameters,
                        run_id=run_id,
                        provenance=provenance,
                        status="skipped",
                        reason=str(capacity.reason),
                        capacity=capacity,
                    )
                    write_result(_result_path(output_dir, spec.label(), run_id, 0, status="skipped"), record)
                    continue
                active.append((spec, parameters, run_id, paths))
            if not active:
                continue

            artifact: BuildArtifact | None = None
            try:
                artifact = build_index(
                    adapter,
                    cp,
                    cache,
                    build_spec,
                    gpu_device=args.gpu_device,
                    scratch_root=scratch_root,
                )
            except Exception as exc:
                reason = f"{type(exc).__name__}: {exc}"
                for spec, parameters, run_id, _ in active:
                    record = _unavailable_record(
                        cache,
                        spec,
                        parameters=parameters,
                        run_id=run_id,
                        provenance=provenance,
                        status="error",
                        reason=reason,
                        capacity=capacities[spec],
                    )
                    write_result(_result_path(output_dir, spec.label(), run_id, 0, status="error"), record)
                continue

            try:
                for spec, parameters, run_id, paths in active:
                    index: Any | None = None
                    try:
                        index, load_seconds = _load_for_search(adapter, cp, artifact, spec, args.gpu_device)
                        measurements = measure_search_repeats(
                            index.search,
                            cache.queries,
                            cache.ground_truth_ids,
                            k=K,
                            warmup_queries=args.warmup_queries,
                            repeats=args.repeats,
                            query_batch_size=args.query_batch_size,
                            gpu_device=args.gpu_device,
                            synchronize=lambda: cp.cuda.Device(args.gpu_device).synchronize(),
                        )
                        for measured in measurements:
                            repeat = int(measured.repeat)
                            if not args.no_resume and _is_complete_result(paths[repeat], run_id):
                                continue
                            peak_rss_values = [value for value in (artifact.peak_rss_bytes, measured.peak_rss_bytes) if value is not None]
                            peak_vram_values = [value for value in (artifact.peak_vram_bytes, measured.peak_vram_bytes) if value is not None]
                            record = make_result_record(
                                backend=BACKEND,
                                dataset_id=cache.dataset_id,
                                dataset_window_count=cache.dataset_window_count,
                                dimension=cache.dimension,
                                parameter_label=spec.label(),
                                parameters=parameters,
                                run_id=run_id,
                                repeat=repeat,
                                warmup_queries=args.warmup_queries,
                                timed_queries=cache.queries.shape[0],
                                query_ids_sha256=cache.query_ids_sha256,
                                ground_truth_ids_sha256=cache.ground_truth_ids_sha256,
                                provenance=provenance,
                                qps=float(measured.qps),
                                latency_ms=measured.latency_ms,
                                recall=float(measured.recall_at_100),
                                index_bytes=artifact.index_bytes,
                                build_seconds=artifact.build_seconds,
                                peak_rss_bytes=max(peak_rss_values) if peak_rss_values else None,
                                peak_vram_bytes=max(peak_vram_values) if peak_vram_values else None,
                                measurement=_measurement_payload(measured, artifact, load_seconds, spec),
                            )
                            write_result(paths[repeat], record)
                    except Exception as exc:
                        reason = f"{type(exc).__name__}: {exc}"
                        record = _unavailable_record(
                            cache,
                            spec,
                            parameters=parameters,
                            run_id=run_id,
                            provenance=provenance,
                            status="error",
                            reason=reason,
                            capacity=capacities[spec],
                        )
                        write_result(_result_path(output_dir, spec.label(), run_id, 0, status="error"), record)
                    finally:
                        if index is not None:
                            index.close()
                        _empty_cupy_pools(cp)
            finally:
                if artifact is not None and not args.keep_index_artifacts:
                    _safe_remove_scratch(artifact.directory, scratch_root)
                _empty_cupy_pools(cp)
    finally:
        _empty_cupy_pools(cp)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
