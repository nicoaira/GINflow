#!/usr/bin/env python3
"""Benchmark the pinned ScaNN frontier against a prepared GINflow cache.

The runner uses GINflow's production ``bin/scann_index.py`` constructor,
serializer, loader, and ``ScannIndex.search`` wrapper.  It therefore measures
the same dot-product / cosine semantics used by ``BUILD_SCANN_INDEX`` and
``SEARCH_SCANN`` rather than a parallel direct-ScaNN implementation.

The default ``auto`` mode follows the production adapter's documented plan:
brute force below 20k vectors, asymmetric hashing below 100k, and spherical
tree + asymmetric hashing + exact reorder above that.  The supplied Rouskin
datasets are in the latter category.  The query-time tree coverage and reorder
ladders follow the official ScaNN algorithms/configuration guidance, which is
also cited in ``docs/indexes.md``.

Run this script in the pinned ScaNN image.  For example:

```
docker run --rm -v "$PWD":/work -w /work \\
  -v /mnt/ssd_samsung/ginflow-benchmarks:/bench-cache \\
  community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91 \\
  python3 benchmarks/run_scann.py --cache-dir /bench-cache --dataset rouskin_6k
```
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import os
import re
import shutil
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


SCANN_IMAGE = "community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91"
SCANN_VERSION = "1.4.2"
BACKEND = "scann"
K = 100
BRUTE_FORCE_MAX = 20_000
AH_ONLY_MAX = 100_000
DEFAULT_REPEATS = 3
DEFAULT_WARMUP_QUERIES = 64
DEFAULT_QUERY_BATCH_SIZE = 32
DEFAULT_MAX_RAM_FRACTION = 0.80
DEFAULT_THREADS = min(6, os.cpu_count() or 1)
DEFAULT_AH_DIM = 2
DEFAULT_ANISOTROPIC = 0.2
DEFAULT_REORDER = 100
SOAR_LAMBDA = 1.5


@dataclass(frozen=True)
class CacheInputs:
    """Validated backend-neutral cache inputs needed for a ScaNN run."""

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
    """One production ScaNN build/search setting.

    ``build_*`` values are serialized in the ScaNN artifact.  The pipeline's
    search adapter can then apply ``leaves_to_search`` and ``reorder`` at query
    time; retaining that split makes a recall/QPS ladder attributable to a
    single built index.
    """

    mode: str
    num_leaves: int | None = None
    leaves_to_search: int | None = None
    reorder: int | None = None
    ah_dim: int | None = None
    anisotropic: float | None = None
    soar: bool = False
    build_leaves_to_search: int | None = None
    build_reorder: int | None = None

    def build_key(self) -> tuple[Any, ...]:
        return (
            self.mode,
            self.num_leaves,
            self.ah_dim,
            self.anisotropic,
            self.soar,
            self.build_leaves_to_search,
            self.build_reorder,
        )

    def parameters(self, *, query_batch_size: int, threads: int) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "scann_leaves": self.num_leaves,
            "scann_leaves_to_search": self.leaves_to_search,
            "scann_reorder": self.reorder,
            "scann_ah_dim": self.ah_dim,
            "scann_anisotropic": self.anisotropic,
            "scann_soar": self.soar,
            "scann_soar_lambda": SOAR_LAMBDA if self.soar else None,
            "build_scann_leaves_to_search": self.build_leaves_to_search,
            "build_scann_reorder": self.build_reorder,
            "scann_num_neighbors": K,
            "query_batch_size": int(query_batch_size),
            "threads_requested": int(threads),
        }

    def label(self) -> str:
        pieces = [self.mode.replace("_", "-")]
        if self.num_leaves is not None:
            pieces.append(f"leaves{self.num_leaves}")
        if self.leaves_to_search is not None:
            pieces.append(f"lts{self.leaves_to_search}")
        if self.reorder is not None:
            pieces.append(f"reorder{self.reorder}")
        if self.ah_dim is not None:
            pieces.append(f"ah{self.ah_dim}")
        if self.soar:
            pieces.append("soar")
        return "-".join(pieces)


@dataclass(frozen=True)
class CapacityCheck:
    feasible: bool
    estimated_bytes: int
    available_bytes: int | None
    safe_limit_bytes: int | None
    reason: str | None = None


@dataclass(frozen=True)
class ThreadControl:
    requested: int
    active_cpus: tuple[int, ...] | None
    initial_cpus: tuple[int, ...] | None
    affinity_applied: bool
    error: str | None = None


@dataclass
class BuildArtifact:
    """A loaded production ScaNN index and its one-time build evidence."""

    index: Any
    details: dict[str, Any]
    build_seconds: float
    serialize_seconds: float
    load_seconds: float
    index_bytes: int
    peak_rss_bytes: int | None
    serialized_directory: Path


def cache_root_from_args(value: str | None) -> Path:
    """Resolve the external cache without creating it implicitly."""

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
    """Resolve a manifest-relative artifact path without escaping its cache phase."""

    value = str(configured or fallback)
    candidate = (directory / value).resolve()
    resolved_directory = directory.resolve()
    if candidate != resolved_directory and resolved_directory not in candidate.parents:
        raise ValueError(f"cache manifest path escapes {directory}: {value}")
    return candidate


def load_cache_inputs(cache_dir: Path, dataset_id: str) -> CacheInputs:
    """Load and cross-check a prepared vector/query/exact-neighbour cache."""

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
    if vectors.ndim != 2 or vectors.shape[0] < K:
        raise ValueError(f"{vectors_path} must be a 2-D database with at least {K} rows")
    if vectors.dtype != np.float32:
        raise ValueError(f"{vectors_path} must contain float32 vectors, got {vectors.dtype}")
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
    if str(flatten.get("dataset_id")) != dataset_id or str(selection.get("dataset_id")) != dataset_id:
        raise ValueError("prepared cache manifests do not match the requested dataset id")
    if str(truth.get("dataset_id")) != dataset_id:
        raise ValueError("exact ground truth manifest does not match the requested dataset id")
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


def _unique_positive(values: Iterable[int], upper: int, *, lower: int = 1) -> list[int]:
    return sorted({max(lower, min(int(value), upper)) for value in values})


def _resolved_mode(mode: str, n_vectors: int) -> str:
    if mode == "auto":
        if n_vectors < BRUTE_FORCE_MAX:
            return "brute_force"
        if n_vectors < AH_ONLY_MAX:
            return "ah"
        return "tree_ah"
    return mode.replace("-", "_")


def _default_num_leaves(n_vectors: int) -> int:
    return max(1, min(n_vectors, int(round(math.sqrt(n_vectors)))))


def _tree_leaf_variants(n_vectors: int, *, plan: str, requested: Sequence[int] | None) -> list[int]:
    if requested:
        invalid = [value for value in requested if value < 1 or value > n_vectors]
        if invalid:
            raise ValueError(f"--tree-leaves values must be between 1 and {n_vectors}: {invalid}")
        return _unique_positive(requested, n_vectors)
    default = _default_num_leaves(n_vectors)
    if plan == "smoke":
        return [default]
    lower = 2 ** max(0, int(math.floor(math.log2(default))))
    upper_sensitivity = min(n_vectors, 2 ** (int(math.ceil(math.log2(default))) + 1))
    # These values evaluate the production sqrt(n) default and the same useful
    # 512/2048 (6k) or 1024/4096 (30k) scale sensitivities in docs/indexes.md.
    return _unique_positive((default, lower, upper_sensitivity), n_vectors)


def _coverage_ladder(num_leaves: int, *, plan: str) -> list[int]:
    values = (4, 8) if plan == "smoke" else (8, 16, 32, 64, 128)
    return _unique_positive(values, num_leaves)


def _reorder_ladder(n_vectors: int, *, plan: str) -> list[int]:
    values = (DEFAULT_REORDER,) if plan == "smoke" else (100, 200, 400, 800)
    return _unique_positive(values, n_vectors, lower=K)


def _dedupe_specs(specs: Iterable[IndexSpec]) -> list[IndexSpec]:
    result: list[IndexSpec] = []
    seen: set[tuple[Any, ...]] = set()
    for spec in specs:
        key = tuple(asdict(spec).items())
        if key not in seen:
            seen.add(key)
            result.append(spec)
    return result


def frontier_specs(
    n_vectors: int,
    *,
    plan: str,
    mode: str = "auto",
    include_brute_force: bool = False,
    include_soar: bool = False,
    tree_leaves: Sequence[int] | None = None,
) -> list[IndexSpec]:
    """Return a compact ScaNN recall/QPS frontier using pipeline-exposed knobs."""

    if n_vectors < K:
        raise ValueError(f"ScaNN benchmark requires at least k={K} vectors")
    if plan not in {"frontier", "smoke"}:
        raise ValueError("plan must be frontier or smoke")
    if mode not in {"auto", "brute-force", "ah", "tree-ah"}:
        raise ValueError("mode must be auto, brute-force, ah, or tree-ah")
    resolved = _resolved_mode(mode, n_vectors)
    # The production adapter has no separate command-line switch for a scoring
    # mode.  It chooses brute/AH/tree from N, except that supplying
    # --scann_leaves forces the tree path.  Do not silently benchmark a direct
    # ScaNN configuration that the pipeline itself cannot construct.
    if mode == "auto" and tree_leaves:
        resolved = "tree_ah"
    if mode == "brute-force" and n_vectors >= BRUTE_FORCE_MAX:
        raise ValueError(
            "--mode brute-force is not exposed by GINflow's production ScaNN adapter "
            f"for n_windows >= {BRUTE_FORCE_MAX}; use the exact shared ground truth instead"
        )
    if mode == "ah" and not (BRUTE_FORCE_MAX <= n_vectors < AH_ONLY_MAX):
        raise ValueError(
            "--mode ah is only production-representable for "
            f"{BRUTE_FORCE_MAX} <= n_windows < {AH_ONLY_MAX}; use --mode tree-ah above that range"
        )
    if include_brute_force and n_vectors >= BRUTE_FORCE_MAX:
        raise ValueError(
            "--include-brute-force is not exposed by GINflow's production ScaNN adapter "
            f"for n_windows >= {BRUTE_FORCE_MAX}; the shared exact ground truth is the reference"
        )
    if include_soar and resolved != "tree_ah":
        raise ValueError("--include-soar requires a tree index; set --mode tree-ah or --tree-leaves")
    specs: list[IndexSpec] = []
    if resolved == "brute_force":
        specs.append(IndexSpec(mode="brute_force"))
    elif resolved == "ah":
        build_reorder = DEFAULT_REORDER
        for reorder in _reorder_ladder(n_vectors, plan=plan):
            specs.append(
                IndexSpec(
                    mode="ah",
                    reorder=reorder,
                    ah_dim=DEFAULT_AH_DIM,
                    anisotropic=DEFAULT_ANISOTROPIC,
                    build_reorder=build_reorder,
                )
            )
    elif resolved == "tree_ah":
        leaves = _tree_leaf_variants(n_vectors, plan=plan, requested=tree_leaves)
        default_leaves = _default_num_leaves(n_vectors)
        primary_leaves = default_leaves if default_leaves in leaves else leaves[0]
        for leaf_count in leaves:
            coverage = _coverage_ladder(leaf_count, plan=plan)
            for leaves_to_search in coverage:
                specs.append(
                    IndexSpec(
                        mode="tree_ah",
                        num_leaves=leaf_count,
                        leaves_to_search=leaves_to_search,
                        reorder=DEFAULT_REORDER,
                        ah_dim=DEFAULT_AH_DIM,
                        anisotropic=DEFAULT_ANISOTROPIC,
                        build_leaves_to_search=coverage[0],
                        build_reorder=DEFAULT_REORDER,
                    )
                )
        primary_coverage = _coverage_ladder(primary_leaves, plan=plan)
        reorder_coverage = min(max(primary_coverage), 64)
        for reorder in _reorder_ladder(n_vectors, plan=plan)[1:]:
            specs.append(
                IndexSpec(
                    mode="tree_ah",
                    num_leaves=primary_leaves,
                    leaves_to_search=reorder_coverage,
                    reorder=reorder,
                    ah_dim=DEFAULT_AH_DIM,
                    anisotropic=DEFAULT_ANISOTROPIC,
                    build_leaves_to_search=primary_coverage[0],
                    build_reorder=DEFAULT_REORDER,
                )
            )
        if include_soar:
            specs.append(
                IndexSpec(
                    mode="tree_ah",
                    num_leaves=primary_leaves,
                    leaves_to_search=reorder_coverage,
                    reorder=DEFAULT_REORDER,
                    ah_dim=DEFAULT_AH_DIM,
                    anisotropic=DEFAULT_ANISOTROPIC,
                    soar=True,
                    build_leaves_to_search=primary_coverage[0],
                    build_reorder=DEFAULT_REORDER,
                )
            )
    else:  # pragma: no cover - kept for future mode additions
        raise AssertionError(f"unexpected ScaNN mode {resolved!r}")

    if include_brute_force and resolved != "brute_force":
        specs.insert(0, IndexSpec(mode="brute_force"))
    return _dedupe_specs(specs)


def group_build_specs(specs: Iterable[IndexSpec]) -> Iterator[tuple[IndexSpec, list[IndexSpec]]]:
    groups: dict[tuple[Any, ...], list[IndexSpec]] = {}
    for spec in specs:
        groups.setdefault(spec.build_key(), []).append(spec)
    for values in groups.values():
        yield values[0], values


def _raw_vector_bytes(n_vectors: int, dimension: int) -> int:
    return int(n_vectors) * int(dimension) * np.dtype(np.float32).itemsize


def _available_ram_bytes() -> int | None:
    meminfo = Path("/proc/meminfo")
    if not meminfo.is_file():
        return None
    match = re.search(r"^MemAvailable:\s+(\d+)\s+kB", meminfo.read_text(errors="replace"), re.MULTILINE)
    return int(match.group(1)) * 1024 if match else None


def capacity_check(
    spec: IndexSpec,
    *,
    n_vectors: int,
    dimension: int,
    available_bytes: int | None,
    max_fraction: float,
) -> CapacityCheck:
    """Conservatively refuse unsafe full-matrix ScaNN construction.

    Production ScaNN accepts the complete float32 matrix and exact reorder
    retains raw vectors.  The estimate reserves three raw matrix equivalents
    plus 1 GiB for hashing/tree/workspace.  It is deliberately a safety gate,
    not a claim about the exact final on-disk ScaNN representation.
    """

    _ = spec
    raw = _raw_vector_bytes(n_vectors, dimension)
    estimated = (3 * raw) + (1024**3)
    if available_bytes is None:
        return CapacityCheck(True, estimated, None, None)
    safe_limit = int(available_bytes * max_fraction)
    if estimated > safe_limit:
        return CapacityCheck(
            False,
            estimated,
            available_bytes,
            safe_limit,
            "estimated ScaNN full-matrix build peak exceeds the configured safe MemAvailable fraction",
        )
    return CapacityCheck(True, estimated, available_bytes, safe_limit)


def _pipeline_scann_adapter() -> Any:
    """Import the production ScaNN adapter without requiring package install."""

    try:
        installed = package_version("scann")
    except PackageNotFoundError as exc:
        raise RuntimeError("ScaNN is not installed; run this runner in the pinned ScaNN container") from exc
    if installed != SCANN_VERSION:
        raise RuntimeError(f"this benchmark requires pinned ScaNN {SCANN_VERSION}, found {installed}")
    bin_dir = Path(__file__).resolve().parents[1] / "bin"
    if str(bin_dir) not in sys.path:
        sys.path.insert(0, str(bin_dir))
    return importlib.import_module("scann_index")


def _configure_threads(requested: int) -> ThreadControl:
    """Apply best-effort process affinity before the native ScaNN runtime loads."""

    if requested < 1:
        raise ValueError("--threads must be >= 1")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = str(requested)
    try:
        initial = tuple(sorted(os.sched_getaffinity(0)))
        selected = tuple(initial[:requested])
        if not selected:
            return ThreadControl(requested, None, initial, False, "no CPU affinity entries are available")
        os.sched_setaffinity(0, set(selected))
        active = tuple(sorted(os.sched_getaffinity(0)))
        return ThreadControl(requested, active, initial, active == selected)
    except (AttributeError, OSError) as exc:
        return ThreadControl(requested, None, None, False, f"CPU affinity unavailable: {exc}")


def _restore_threads(control: ThreadControl) -> None:
    if control.affinity_applied and control.initial_cpus is not None:
        try:
            os.sched_setaffinity(0, set(control.initial_cpus))
        except (AttributeError, OSError):
            pass


def _production_options(adapter: Any, spec: IndexSpec) -> Any:
    """Construct exactly the IndexOptions shape consumed by the pipeline adapter."""

    return adapter.IndexOptions(
        index_type="scann",
        gpu=False,
        scann_reorder=int(spec.build_reorder or DEFAULT_REORDER),
        scann_ah_dim=int(spec.ah_dim or DEFAULT_AH_DIM),
        scann_anisotropic=float(spec.anisotropic or DEFAULT_ANISOTROPIC),
        scann_soar=bool(spec.soar),
        scann_num_neighbors=K,
        scann_leaves=spec.num_leaves,
        scann_leaves_to_search=spec.build_leaves_to_search,
    )


def _directory_bytes(directory: Path) -> int:
    return sum(path.stat().st_size for path in directory.rglob("*") if path.is_file())


def build_index(
    adapter: Any,
    cache: CacheInputs,
    spec: IndexSpec,
    *,
    scratch_root: Path,
) -> BuildArtifact:
    """Build, serialize, and reload a production ScaNN index for warm search."""

    options = _production_options(adapter, spec)
    scratch_root.mkdir(parents=True, exist_ok=True)
    serialized_directory = Path(tempfile.mkdtemp(prefix="scann-", dir=scratch_root))
    built_index: Any | None = None
    try:
        with PeakResourceSampler() as sampler:
            started = time.perf_counter()
            built_index, details = adapter.build_populated_searcher(cache.vectors, options)
            build_seconds = time.perf_counter() - started
        if int(built_index.ntotal) != cache.dataset_window_count:
            raise RuntimeError(
                f"ScaNN built {built_index.ntotal} vectors, expected {cache.dataset_window_count}"
            )
        serialize_started = time.perf_counter()
        adapter.serialize_index(built_index, serialized_directory)
        serialize_seconds = time.perf_counter() - serialize_started
        index_bytes = _directory_bytes(serialized_directory)
        load_started = time.perf_counter()
        loaded_index = adapter.load_index(
            serialized_directory,
            ntotal=cache.dataset_window_count,
            leaves_to_search=spec.build_leaves_to_search,
            reorder=spec.build_reorder,
        )
        load_seconds = time.perf_counter() - load_started
        del built_index
        return BuildArtifact(
            index=loaded_index,
            details=dict(details),
            build_seconds=float(build_seconds),
            serialize_seconds=float(serialize_seconds),
            load_seconds=float(load_seconds),
            index_bytes=int(index_bytes),
            peak_rss_bytes=sampler.peak_rss_bytes,
            serialized_directory=serialized_directory,
        )
    except Exception:
        if built_index is not None:
            del built_index
        _remove_scratch_directory(serialized_directory, scratch_root)
        raise


def _remove_scratch_directory(path: Path, scratch_root: Path) -> None:
    """Remove only a fresh runner-owned temporary serialized index directory."""

    if not path.exists():
        return
    if path.parent.resolve() != scratch_root.resolve() or not path.name.startswith("scann-"):
        raise RuntimeError(f"refusing to remove a non-runner ScaNN scratch path: {path}")
    shutil.rmtree(path)


def _set_search_parameters(adapter: Any, index: Any, spec: IndexSpec) -> None:
    adapter.apply_search_params(index, nprobe=spec.leaves_to_search, reorder=spec.reorder)


def _run_id(
    cache: CacheInputs,
    spec: IndexSpec,
    parameters: Mapping[str, Any],
    hardware_id: str,
    image: str,
) -> str:
    return stable_id(
        "scann-run",
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


def _base_provenance(
    cache: CacheInputs,
    *,
    repository_root: Path,
    hardware: Mapping[str, Any],
    image: str,
    container: str,
    thread_control: ThreadControl | None,
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
        "scann_version": SCANN_VERSION,
        "index_constructor": "bin/scann_index.py::build_populated_searcher",
        "index_serializer": "bin/scann_index.py::serialize_index",
        "index_loader": "bin/scann_index.py::load_index",
        "search_adapter": "bin/scann_index.py::ScannIndex.search",
        "benchmark_build_method": "memmap-through-production-scann-adapter-v1",
        "production_build_difference": (
            "The benchmark passes flat/vectors.npy as a float32 memmap to the production "
            "adapter; BUILD_SCANN_INDEX currently receives a pre-concatenated in-memory matrix. "
            "Search semantics and serialized/load behavior are production code, but benchmark "
            "build RSS is not a direct pipeline peak-RSS requirement."
        ),
        "vram_measurement": "not applicable: ScaNN is CPU-only",
        "source_references": [
            "https://github.com/google-research/google-research/blob/master/scann/docs/algorithms.md",
            "https://arxiv.org/abs/1908.10396",
        ],
    }
    if thread_control is not None:
        provenance["thread_control"] = {
            "requested": thread_control.requested,
            "active_cpu_ids": list(thread_control.active_cpus) if thread_control.active_cpus is not None else None,
            "affinity_applied": thread_control.affinity_applied,
            "note": (
                "best-effort process affinity and thread environment; the current production "
                "ScaNN adapter does not expose ScaNN's set_n_training_threads()"
            ),
            "error": thread_control.error,
        }
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
            "estimated_bytes": capacity.estimated_bytes,
            "available_bytes": capacity.available_bytes,
            "safe_limit_bytes": capacity.safe_limit_bytes,
            "method": "3x raw float32 matrix plus 1 GiB safety reserve",
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


def _measurement_payload(measured: Any, build: BuildArtifact) -> dict[str, Any]:
    payload = dict(getattr(measured, "measurement", {}) or {})
    payload.update(
        {
            "timing_scope": "warm serialized-and-reloaded ScaNN search; build, serialization, and load excluded from QPS",
            "index_constructor": "bin/scann_index.py::build_populated_searcher",
            "index_search": "bin/scann_index.py::ScannIndex.search -> search_batched",
            "build_peak_rss_bytes": build.peak_rss_bytes,
            "build_seconds_scope": "build_populated_searcher only",
            "serialize_seconds": build.serialize_seconds,
            "load_seconds": build.load_seconds,
            "serialized_index_bytes_method": "sum of files produced by bin/scann_index.py::serialize_index",
            "production_build_rss_comparable": False,
            "peak_vram_scope": "not applicable: ScaNN is CPU-only",
        }
    )
    return payload


def _parse_tree_leaves(value: str | None) -> list[int] | None:
    if value is None:
        return None
    try:
        values = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise ValueError("--tree-leaves must be comma-separated integers") from exc
    if not values:
        raise ValueError("--tree-leaves must contain at least one integer")
    return values


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    parser.add_argument("--dataset", required=True, help="prepared cache dataset id")
    parser.add_argument("--output-dir", type=Path, help="defaults to <cache>/<dataset>/results/scann")
    parser.add_argument("--plan", choices=("frontier", "smoke"), default="frontier")
    parser.add_argument("--mode", choices=("auto", "brute-force", "ah", "tree-ah"), default="auto")
    parser.add_argument("--include-brute-force", action="store_true", help="add an explicit exact ScaNN baseline")
    parser.add_argument("--include-soar", action="store_true", help="add one separate production SOAR tree variant")
    parser.add_argument("--tree-leaves", help="override tree build leaf counts with comma-separated values")
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--warmup-queries", type=int, default=DEFAULT_WARMUP_QUERIES)
    parser.add_argument("--query-batch-size", type=int, default=DEFAULT_QUERY_BATCH_SIZE)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="best-effort CPU affinity/thread environment")
    parser.add_argument("--max-ram-fraction", type=float, default=DEFAULT_MAX_RAM_FRACTION)
    parser.add_argument("--image", default=SCANN_IMAGE, help="pinned runtime image recorded in provenance")
    parser.add_argument("--container", default="docker", help="runtime type recorded in provenance")
    parser.add_argument("--keep-index-artifacts", action="store_true", help="retain runner-owned serialized indexes under cache/scratch")
    parser.add_argument("--no-resume", action="store_true", help="re-measure settings with complete result records")
    parser.add_argument("--dry-run", action="store_true", help="validate the cache and print the plan without importing ScaNN")
    return parser.parse_args(argv)


def _validate_args(args: argparse.Namespace, cache: CacheInputs) -> None:
    for name in ("repeats", "warmup_queries", "query_batch_size", "threads"):
        if int(getattr(args, name)) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 1")
    if not 0.0 < float(args.max_ram_fraction) <= 1.0:
        raise ValueError("--max-ram-fraction must be in (0, 1]")
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
        "schema_version": "ginflow-scann-benchmark-plan-v1",
        "dataset_id": cache.dataset_id,
        "dataset_window_count": cache.dataset_window_count,
        "dimension": cache.dimension,
        "image": image,
        "hardware_id": hardware["id"],
        "plan": args.plan,
        "mode": args.mode,
        "settings": [
            {
                "label": spec.label(),
                "spec": asdict(spec),
                "parameters": spec.parameters(
                    query_batch_size=args.query_batch_size,
                    threads=args.threads,
                ),
                "capacity": asdict(capacities[spec]),
            }
            for spec in specs
        ],
    }


def run(args: argparse.Namespace) -> int:
    repository_root = Path(__file__).resolve().parents[1]
    cache_dir = cache_root_from_args(args.cache_dir)
    cache = load_cache_inputs(cache_dir, args.dataset)
    _validate_args(args, cache)
    tree_leaves = _parse_tree_leaves(args.tree_leaves)
    specs = frontier_specs(
        cache.dataset_window_count,
        plan=args.plan,
        mode=args.mode,
        include_brute_force=args.include_brute_force,
        include_soar=args.include_soar,
        tree_leaves=tree_leaves,
    )
    if not specs:
        raise ValueError("the selected ScaNN plan produced no benchmark settings")
    hardware = hardware_snapshot()
    available_ram = _available_ram_bytes()
    capacities = {
        spec: capacity_check(
            spec,
            n_vectors=cache.dataset_window_count,
            dimension=cache.dimension,
            available_bytes=available_ram,
            max_fraction=args.max_ram_fraction,
        )
        for spec in specs
    }
    if args.dry_run:
        print(json.dumps(_plan_payload(cache, specs, capacities, hardware=hardware, image=args.image, args=args), indent=2, sort_keys=True))
        return 0

    output_dir = (args.output_dir or (cache.root / "results" / BACKEND)).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    initial_provenance = _base_provenance(
        cache,
        repository_root=repository_root,
        hardware=hardware,
        image=args.image,
        container=args.container,
        thread_control=None,
    )
    active_by_build: dict[tuple[Any, ...], list[tuple[IndexSpec, dict[str, Any], str, list[Path]]]] = {}
    for spec in specs:
        parameters = spec.parameters(query_batch_size=args.query_batch_size, threads=args.threads)
        parameters["max_ram_fraction"] = float(args.max_ram_fraction)
        run_id = _run_id(cache, spec, parameters, str(hardware["id"]), args.image)
        paths = [_result_path(output_dir, spec.label(), run_id, repeat) for repeat in range(args.repeats)]
        if not args.no_resume and all(_is_complete_result(path, run_id) for path in paths):
            continue
        capacity = capacities[spec]
        if not capacity.feasible:
            skipped_path = _result_path(output_dir, spec.label(), run_id, 0, status="skipped")
            if args.no_resume or not skipped_path.is_file():
                record = _unavailable_record(
                    cache,
                    spec,
                    parameters=parameters,
                    run_id=run_id,
                    provenance=initial_provenance,
                    status="skipped",
                    reason=str(capacity.reason),
                    capacity=capacity,
                )
                write_result(skipped_path, record)
            continue
        active_by_build.setdefault(spec.build_key(), []).append((spec, parameters, run_id, paths))
    if not active_by_build:
        return 0

    thread_control = _configure_threads(args.threads)
    provenance = _base_provenance(
        cache,
        repository_root=repository_root,
        hardware=hardware,
        image=args.image,
        container=args.container,
        thread_control=thread_control,
    )
    try:
        try:
            adapter = _pipeline_scann_adapter()
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            for settings in active_by_build.values():
                for spec, parameters, run_id, _paths in settings:
                    error_path = _result_path(output_dir, spec.label(), run_id, 0, status="error")
                    write_result(
                        error_path,
                        _unavailable_record(
                            cache,
                            spec,
                            parameters=parameters,
                            run_id=run_id,
                            provenance=provenance,
                            status="error",
                            reason=reason,
                            capacity=None,
                        ),
                    )
            return 0

        scratch_root = cache.root / "scratch" / BACKEND
        for settings in active_by_build.values():
            build_spec = settings[0][0]
            build: BuildArtifact | None = None
            try:
                build = build_index(adapter, cache, build_spec, scratch_root=scratch_root)
            except Exception as exc:
                reason = f"{type(exc).__name__}: {exc}"
                for spec, parameters, run_id, _paths in settings:
                    error_path = _result_path(output_dir, spec.label(), run_id, 0, status="error")
                    write_result(
                        error_path,
                        _unavailable_record(
                            cache,
                            spec,
                            parameters=parameters,
                            run_id=run_id,
                            provenance=provenance,
                            status="error",
                            reason=reason,
                            capacity=None,
                        ),
                    )
                continue
            try:
                for spec, parameters, run_id, paths in settings:
                    try:
                        _set_search_parameters(adapter, build.index, spec)
                        samples = measure_search_repeats(
                            build.index.search,
                            cache.queries,
                            cache.ground_truth_ids,
                            k=K,
                            warmup_queries=args.warmup_queries,
                            repeats=args.repeats,
                            query_batch_size=args.query_batch_size,
                        )
                        for measured in samples:
                            repeat = int(measured.repeat)
                            if not args.no_resume and _is_complete_result(paths[repeat], run_id):
                                continue
                            peak_values = [value for value in (build.peak_rss_bytes, measured.peak_rss_bytes) if value is not None]
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
                                index_bytes=build.index_bytes,
                                build_seconds=build.build_seconds,
                                peak_rss_bytes=max(peak_values) if peak_values else None,
                                peak_vram_bytes=None,
                                measurement=_measurement_payload(measured, build),
                            )
                            write_result(paths[repeat], record)
                    except Exception as exc:
                        reason = f"{type(exc).__name__}: {exc}"
                        error_path = _result_path(output_dir, spec.label(), run_id, 0, status="error")
                        write_result(
                            error_path,
                            _unavailable_record(
                                cache,
                                spec,
                                parameters=parameters,
                                run_id=run_id,
                                provenance=provenance,
                                status="error",
                                reason=reason,
                                capacity=None,
                            ),
                        )
            finally:
                if build is not None:
                    del build.index
                    if not args.keep_index_artifacts:
                        _remove_scratch_directory(build.serialized_directory, scratch_root)
    finally:
        _restore_threads(thread_control)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
