#!/usr/bin/env python3
"""Benchmark GINflow's pinned NGT, QG, and QBG index implementations.

The default frontier deliberately invokes ``bin/ngt_index.py`` unchanged.  It
therefore measures the three structures exposed by GINflow, rather than a
nearby NGT configuration that users cannot reproduce through the pipeline.
NGT's native epsilon, radius, edge, and quantization controls are available
only with ``--include-native-frontier``.  Those records are visibly marked as
experimental and are never combined with pipeline-faithful points.

Run inside the pinned NGT 2.3.12 image, with the prepared cache mounted on a
large external volume.  For example:

    docker run --rm -v "$PWD":/work -w /work \\
      -v /mnt/ssd_samsung/ginflow-benchmarks:/bench-cache \\
      community.wave.seqera.io/library/python_numpy_ngt:9a0ca7a46e9c18b2 \\
      python3 benchmarks/run_ngt.py --cache-dir /bench-cache --dataset rouskin_6k

NGT's quantized builders can reserve substantial native workspace.  The
conservative preflight records an explicit skipped result instead of allowing
an unsafe build to be killed by the operating system.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np

from benchmark_utils import (
    RUNNER_VERSION,
    PeakResourceSampler,
    hardware_snapshot,
    make_result_record,
    measure_search_repeats,
    software_snapshot,
    stable_id,
    validate_result_record,
    write_result,
)
from run_faiss import CacheInputs, cache_root_from_args, load_cache_inputs


NGT_IMAGE = "community.wave.seqera.io/library/python_numpy_ngt:9a0ca7a46e9c18b2"
BACKEND = "ngt"
K = 100
DEFAULT_REPEATS = 3
DEFAULT_WARMUP_QUERIES = 64
DEFAULT_QUERY_BATCH_SIZE = 32
DEFAULT_MAX_RAM_FRACTION = 0.70
DEFAULT_MAX_DISK_FRACTION = 0.80

# The qbg executable reports a large process virtual-memory footprint even on
# tiny input.  These terms are intentionally capacity-planning reservations,
# not claims about a serialized index's size.  They protect an unattended run
# on the stated 62 GiB workstation when NGT's native allocation is opaque.
NATIVE_WORKSPACE_BYTES = {
    "NGT": 1 * 1024**3,
    "QG": 10 * 1024**3,
    "QBG": 24 * 1024**3,
}

OFFICIAL_NGT_SOURCES = (
    "https://github.com/NGT-labs/NGT/blob/v2.3.12/bin/ngt/README.md",
    "https://github.com/NGT-labs/NGT/blob/v2.3.12/bin/qbg/README.md",
    "https://github.com/NGT-labs/NGT/blob/v2.3.12/python/README.md",
)


@dataclass(frozen=True)
class NgtSpec:
    """One NGT build/search point.

    ``mode=pipeline`` means every non-type field is ``None`` and construction
    is delegated to ``bin/ngt_index.py::build_populated_index``.  ``native``
    specs deliberately opt into documented NGT controls outside GINflow's
    current command-line surface.
    """

    index_type: str
    mode: str = "pipeline"
    edge_size_for_creation: int | None = None
    edge_size_for_search: int | None = None
    qg_subvector_dimension: int | None = None
    qg_max_edges: int | None = None
    qg_trials: int | None = None
    qbg_max_edges: int | None = None
    qbg_trials: int | None = None
    qbg_local_divisions: int | None = None
    epsilon: float | None = None
    search_radius: float | None = None
    search_edge_size: int | None = None
    expected_accuracy: float | None = None
    result_expansion: float | None = None
    blob_epsilon: float | None = None
    exploration_size: int | None = None
    num_of_probes: int | None = None
    exact_result_expansion: float | None = None

    def __post_init__(self) -> None:
        kind = self.index_type.upper()
        if kind not in {"NGT", "QG", "QBG"}:
            raise ValueError("index_type must be NGT, QG, or QBG")
        if self.mode not in {"pipeline", "native"}:
            raise ValueError("mode must be pipeline or native")
        if self.mode == "pipeline" and any(
            value is not None for key, value in asdict(self).items() if key not in {"index_type", "mode"}
        ):
            raise ValueError("pipeline NGT specs must not override hidden NGT controls")
        if self.qg_max_edges is not None and self.qg_max_edges % 16:
            raise ValueError("qg_max_edges must be a multiple of 16")

    @property
    def kind(self) -> str:
        return self.index_type.upper()

    def build_key(self) -> tuple[Any, ...]:
        return (
            self.kind,
            self.mode,
            self.edge_size_for_creation,
            self.edge_size_for_search,
            self.qg_subvector_dimension,
            self.qg_max_edges,
            self.qg_trials,
            self.qbg_max_edges,
            self.qbg_trials,
            self.qbg_local_divisions,
        )

    def parameters(self, *, query_batch_size: int) -> dict[str, Any]:
        values = asdict(self)
        values["index_type"] = self.kind
        values["query_batch_size"] = int(query_batch_size)
        values["pipeline_faithful"] = self.mode == "pipeline"
        return values

    def label(self) -> str:
        parts = [self.mode, self.kind.lower()]
        if self.edge_size_for_creation is not None:
            parts.append(f"build-e{self.edge_size_for_creation}")
        if self.edge_size_for_search is not None:
            parts.append(f"build-s{self.edge_size_for_search}")
        if self.qg_subvector_dimension is not None:
            parts.append(f"q{self.qg_subvector_dimension}")
        if self.qg_max_edges is not None:
            parts.append(f"qg-e{self.qg_max_edges}")
        if self.qbg_max_edges is not None:
            parts.append(f"qbg-e{self.qbg_max_edges}")
        if self.epsilon is not None:
            parts.append(f"eps{self.epsilon:g}")
        if self.search_edge_size is not None:
            parts.append(f"search-e{self.search_edge_size}")
        if self.result_expansion is not None:
            parts.append(f"expand{self.result_expansion:g}")
        if self.exploration_size is not None:
            parts.append(f"explore{self.exploration_size}")
        if self.num_of_probes is not None:
            parts.append(f"probes{self.num_of_probes}")
        return "-".join(parts)


@dataclass(frozen=True)
class CapacityCheck:
    feasible: bool
    estimated_ram_bytes: int
    available_ram_bytes: int | None
    safe_ram_bytes: int | None
    estimated_disk_bytes: int
    free_disk_bytes: int | None
    safe_disk_bytes: int | None
    reason: str | None = None


@dataclass
class BuildArtifact:
    index: Any
    details: dict[str, Any]
    build_seconds: float
    serialize_seconds: float
    load_seconds: float
    index_bytes: int
    peak_rss_bytes: int | None
    peak_vram_bytes: int | None
    child_max_rss_bytes: int | None
    temporary_root: Path


def _unique(values: Iterable[int]) -> list[int]:
    return sorted({int(value) for value in values if int(value) > 0})


def _native_specs(plan: str, dimension: int) -> list[NgtSpec]:
    """Return documented, explicitly non-pipeline NGT sensitivity points."""

    if dimension % 8:
        qg_subvector = 1
    else:
        qg_subvector = 8
    if plan == "smoke":
        return [
            NgtSpec("NGT", "native", edge_size_for_creation=10, edge_size_for_search=40, epsilon=0.0),
            NgtSpec("NGT", "native", edge_size_for_creation=10, edge_size_for_search=40, epsilon=0.1),
            NgtSpec(
                "QG",
                "native",
                edge_size_for_creation=10,
                edge_size_for_search=40,
                qg_subvector_dimension=qg_subvector,
                epsilon=0.02,
                result_expansion=3.0,
            ),
            NgtSpec(
                "QBG",
                "native",
                epsilon=0.02,
                exploration_size=256,
                num_of_probes=1,
            ),
        ]

    specs: list[NgtSpec] = []
    for build_edges, search_edges in ((10, 40), (20, 80)):
        for epsilon in (0.0, 0.05, 0.1, 0.2):
            specs.append(
                NgtSpec(
                    "NGT",
                    "native",
                    edge_size_for_creation=build_edges,
                    edge_size_for_search=search_edges,
                    epsilon=epsilon,
                    search_edge_size=search_edges,
                )
            )
    for epsilon, expansion in ((0.0, 1.0), (0.02, 3.0), (0.05, 3.0), (0.1, 5.0)):
        specs.append(
            NgtSpec(
                "QG",
                "native",
                edge_size_for_creation=10,
                edge_size_for_search=40,
                qg_subvector_dimension=qg_subvector,
                qg_max_edges=16,
                epsilon=epsilon,
                result_expansion=expansion,
            )
        )
    for epsilon, exploration, probes in ((0.0, 128, 1), (0.02, 256, 1), (0.05, 256, 2), (0.1, 512, 2)):
        specs.append(
            NgtSpec(
                "QBG",
                "native",
                epsilon=epsilon,
                exploration_size=exploration,
                num_of_probes=probes,
            )
        )
    return specs


def frontier_specs(
    n_vectors: int,
    dimension: int,
    *,
    plan: str,
    include_native_frontier: bool = False,
) -> list[NgtSpec]:
    """Return faithful GINflow variants and optional labelled native variants."""

    if n_vectors < K:
        raise ValueError("NGT benchmark needs at least 100 vectors")
    if dimension < 1:
        raise ValueError("dimension must be positive")
    if plan not in {"frontier", "smoke"}:
        raise ValueError("plan must be frontier or smoke")
    specs = [NgtSpec("NGT"), NgtSpec("QG"), NgtSpec("QBG")]
    if include_native_frontier:
        specs.extend(_native_specs(plan, dimension))
    return specs


def group_build_specs(specs: Iterable[NgtSpec]) -> Iterator[tuple[NgtSpec, list[NgtSpec]]]:
    grouped: dict[tuple[Any, ...], list[NgtSpec]] = defaultdict(list)
    for spec in specs:
        grouped[spec.build_key()].append(spec)
    for group in grouped.values():
        yield group[0], group


def _available_ram_bytes() -> int | None:
    memory_info = Path("/proc/meminfo")
    if not memory_info.is_file():
        return None
    match = re.search(r"^MemAvailable:\s+(\d+)\s+kB", memory_info.read_text(errors="replace"), re.MULTILINE)
    return int(match.group(1)) * 1024 if match else None


def _free_disk_bytes(path: Path) -> int | None:
    try:
        return int(shutil.disk_usage(path).free)
    except OSError:
        return None


def _raw_vector_bytes(n_vectors: int, dimension: int) -> int:
    return int(n_vectors) * int(dimension) * np.dtype(np.float32).itemsize


def capacity_check(
    spec: NgtSpec,
    *,
    n_vectors: int,
    dimension: int,
    available_ram_bytes: int | None,
    free_disk_bytes: int | None,
    max_ram_fraction: float,
    max_disk_fraction: float,
) -> CapacityCheck:
    """Conservatively gate the opaque, native NGT construction workspace."""

    raw = _raw_vector_bytes(n_vectors, dimension)
    graph_bytes = n_vectors * max(spec.edge_size_for_creation or 10, spec.edge_size_for_search or 40) * 16
    workspace = NATIVE_WORKSPACE_BYTES[spec.kind]
    if spec.kind == "NGT":
        # Production helper makes a contiguous database copy, then NGT stores
        # the objects and graph.  This is deliberately larger than final disk.
        estimated_ram = (2 * raw) + graph_bytes + workspace
        # The production build writes an index and then copies it to the
        # published database.  The benchmark exercises that same round trip,
        # so reserve both source and serialized copies.
        estimated_disk = (2 * (raw + graph_bytes)) + (512 * 1024**2)
    elif spec.kind == "QG":
        # QG first retains a regular NGT graph and then qbg's quantizer builds
        # a second representation.  Its native workspace is non-linear.
        estimated_ram = (3 * raw) + graph_bytes + workspace
        estimated_disk = (4 * raw) + (2 * graph_bytes) + (2 * 1024**3)
    else:
        # Pipeline QBG writes a text staging matrix through np.savetxt before
        # qbg builds.  Assume up to four bytes of text per raw binary byte.
        estimated_ram = (3 * raw) + workspace
        estimated_disk = (6 * raw) + (4 * 1024**3)
    safe_ram = int(available_ram_bytes * max_ram_fraction) if available_ram_bytes is not None else None
    safe_disk = int(free_disk_bytes * max_disk_fraction) if free_disk_bytes is not None else None
    if safe_ram is not None and estimated_ram > safe_ram:
        return CapacityCheck(
            False,
            estimated_ram,
            available_ram_bytes,
            safe_ram,
            estimated_disk,
            free_disk_bytes,
            safe_disk,
            "estimated NGT native build peak exceeds the configured safe MemAvailable fraction",
        )
    if safe_disk is not None and estimated_disk > safe_disk:
        return CapacityCheck(
            False,
            estimated_ram,
            available_ram_bytes,
            safe_ram,
            estimated_disk,
            free_disk_bytes,
            safe_disk,
            "estimated NGT temporary/index disk use exceeds the configured safe free-disk fraction",
        )
    return CapacityCheck(
        True,
        estimated_ram,
        available_ram_bytes,
        safe_ram,
        estimated_disk,
        free_disk_bytes,
        safe_disk,
    )


def _import_ngtpy() -> Any:
    try:
        import ngtpy
        from importlib.metadata import version
    except ImportError as exc:  # pragma: no cover - exercised in the pinned image
        raise RuntimeError("NGT is not installed; run this benchmark in the pinned NGT image") from exc
    installed = str(version("ngt"))
    if installed != "2.3.12":
        raise RuntimeError(f"this benchmark requires NGT 2.3.12, found {installed or 'unknown'}")
    if shutil.which("qbg") is None:
        raise RuntimeError("pinned NGT qbg executable is unavailable")
    return ngtpy


def _pipeline_ngt_adapter() -> Any:
    bin_dir = Path(__file__).resolve().parents[1] / "bin"
    if str(bin_dir) not in sys.path:
        sys.path.insert(0, str(bin_dir))
    return importlib.import_module("ngt_index")


def _directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _remove_owned_directory(path: Path, temporary_parent: Path) -> None:
    """Remove one tempfile-created runner directory, never a user path."""

    try:
        root = path.resolve()
        parent = temporary_parent.resolve()
    except OSError:
        return
    if root == parent or parent not in root.parents or not root.name.startswith("ginflow-ngt-"):
        raise RuntimeError(f"refusing to remove a non-runner NGT scratch path: {path}")
    shutil.rmtree(root, ignore_errors=True)


def _children_peak_rss_bytes() -> int | None:
    try:
        value = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    except (AttributeError, ValueError):
        return None
    if value <= 0:
        return None
    return value if sys.platform == "darwin" else value * 1024


@contextlib.contextmanager
def _temporary_directory_on(path: Path) -> Iterator[None]:
    """Make the production helper place its disposable build tree externally."""

    path.mkdir(parents=True, exist_ok=True)
    old_tempdir = tempfile.tempdir
    old_environment = {name: os.environ.get(name) for name in ("TMPDIR", "TEMP", "TMP")}
    tempfile.tempdir = str(path)
    for name in old_environment:
        os.environ[name] = str(path)
    try:
        yield
    finally:
        tempfile.tempdir = old_tempdir
        for name, value in old_environment.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def _run_qbg(*args: str) -> None:
    subprocess.run(["qbg", *args], check=True)


def _build_native_index(
    ngtpy: Any,
    adapter: Any,
    cache: CacheInputs,
    spec: NgtSpec,
    *,
    temporary_parent: Path,
) -> tuple[Any, dict[str, Any], Path]:
    """Build one explicitly experimental, native NGT variant on external disk."""

    root = Path(tempfile.mkdtemp(prefix="ginflow-ngt-native-", dir=temporary_parent))
    index_path = root / "index"
    vectors = np.ascontiguousarray(cache.vectors, dtype=np.float32)
    try:
        if spec.kind in {"NGT", "QG"}:
            ngtpy.create(
                str(index_path),
                cache.dimension,
                edge_size_for_creation=int(spec.edge_size_for_creation or 10),
                edge_size_for_search=int(spec.edge_size_for_search or 40),
                distance_type="Cosine",
                object_type="Float",
            )
            regular = ngtpy.Index(str(index_path))
            try:
                regular.batch_insert(vectors)
                regular.save()
            finally:
                regular.close()
            if spec.kind == "QG":
                create_args = ["create-qg"]
                if spec.qg_subvector_dimension is not None:
                    create_args.extend(["-Q", str(spec.qg_subvector_dimension)])
                create_args.append(str(index_path))
                _run_qbg(*create_args)
                build_args = ["build-qg", "-o", str(cache.dataset_window_count)]
                if spec.qg_max_edges is not None:
                    build_args.extend(["-E", str(spec.qg_max_edges)])
                if spec.qg_trials is not None:
                    build_args.extend(["-M", str(spec.qg_trials)])
                build_args.append(str(index_path))
                _run_qbg(*build_args)
        else:
            data_path = root / "vectors.tsv"
            np.savetxt(data_path, vectors, fmt="%.9g")
            create_args = ["create", "-d", str(cache.dimension), "-D", "2"]
            if spec.qbg_local_divisions is not None:
                create_args.extend(["-N", str(spec.qbg_local_divisions)])
            create_args.append(str(index_path))
            _run_qbg(*create_args)
            _run_qbg("append", str(index_path), str(data_path))
            build_args = ["build", "-o", str(cache.dataset_window_count)]
            if spec.qbg_max_edges is not None:
                build_args.extend(["-E", str(spec.qbg_max_edges)])
            if spec.qbg_trials is not None:
                build_args.extend(["-M", str(spec.qbg_trials)])
            build_args.append(str(index_path))
            _run_qbg(*build_args)
        index = adapter.NgtIndex(index_path, spec.kind, cache.dataset_window_count)
        details = {
            "backend": "ngt",
            "index_type": spec.kind,
            "index_class": spec.kind,
            "metric": index.metric,
            "gpu": False,
            "gpu_capable": False,
            "path": index_path,
            "benchmark_native_controls": True,
        }
        return index, details, root
    except Exception:
        shutil.rmtree(root, ignore_errors=True)
        raise
    finally:
        del vectors


def _build_index(
    ngtpy: Any,
    adapter: Any,
    cache: CacheInputs,
    spec: NgtSpec,
    *,
    temporary_parent: Path,
) -> BuildArtifact:
    """Build a faithful or labelled-native index and preserve it for searches."""

    temporary_parent.mkdir(parents=True, exist_ok=True)
    children_before = _children_peak_rss_bytes()
    source_index: Any | None = None
    loaded_index: Any | None = None
    serialized_root: Path | None = None
    with PeakResourceSampler() as sampler:
        try:
            started = time.perf_counter()
            if spec.mode == "pipeline":
                with _temporary_directory_on(temporary_parent):
                    source_index, details = adapter.build_populated_index(cache.vectors, spec.kind.lower())
            else:
                source_index, details, _ = _build_native_index(
                    ngtpy,
                    adapter,
                    cache,
                    spec,
                    temporary_parent=temporary_parent,
                )
            build_seconds = time.perf_counter() - started
            if int(source_index.ntotal) != cache.dataset_window_count:
                raise RuntimeError(f"NGT added {source_index.ntotal} vectors, expected {cache.dataset_window_count}")

            # BUILD_NGT_INDEX serializes this wrapper and SEARCH_NGT reopens it
            # in a distinct task.  Search the reopened representation here so
            # a benchmark point has the same artifact boundary as the pipeline.
            serialized_root = Path(tempfile.mkdtemp(prefix="ginflow-ngt-serialized-", dir=temporary_parent))
            serialized_path = serialized_root / "index"
            serialize_started = time.perf_counter()
            adapter.serialize_index(source_index, serialized_path)
            source_index = None
            serialize_seconds = time.perf_counter() - serialize_started
            meta = adapter.meta_from_details(dict(details))
            meta["n_windows"] = cache.dataset_window_count
            load_started = time.perf_counter()
            loaded_index = adapter.load_index(serialized_path, meta)
            load_seconds = time.perf_counter() - load_started
            if int(loaded_index.ntotal) != cache.dataset_window_count:
                raise RuntimeError(
                    f"reloaded NGT index has {loaded_index.ntotal} vectors, expected {cache.dataset_window_count}"
                )
            index_bytes = _directory_size(serialized_path)
        except Exception:
            if loaded_index is not None:
                close = getattr(loaded_index, "close", None)
                if close is not None:
                    close()
            if source_index is not None:
                close = getattr(source_index, "close", None)
                if close is not None:
                    close()
                try:
                    _remove_owned_directory(Path(source_index.path).parent, temporary_parent)
                except (OSError, RuntimeError):
                    pass
            if serialized_root is not None:
                try:
                    _remove_owned_directory(serialized_root, temporary_parent)
                except (OSError, RuntimeError):
                    pass
            raise
    children_after = _children_peak_rss_bytes()
    child_peak = children_after
    if children_before is not None and children_after is not None and children_after < children_before:
        child_peak = children_before
    return BuildArtifact(
        index=loaded_index,
        details=dict(details),
        build_seconds=build_seconds,
        serialize_seconds=serialize_seconds,
        load_seconds=load_seconds,
        index_bytes=index_bytes,
        peak_rss_bytes=sampler.peak_rss_bytes,
        peak_vram_bytes=sampler.peak_vram_bytes,
        child_max_rss_bytes=child_peak,
        temporary_root=serialized_root,
    )


def _native_search_set(index: Any, spec: NgtSpec) -> None:
    """Apply only documented ngtpy controls after native construction."""

    if spec.mode != "native":
        return
    raw_index = getattr(index, "_index", None)
    if raw_index is None:
        raise RuntimeError("NGT wrapper has no native search index")
    if spec.kind in {"NGT", "QG"}:
        values: dict[str, Any] = {}
        for name, value in (
            ("epsilon", spec.epsilon),
            ("search_radius", spec.search_radius),
            ("edge_size", spec.search_edge_size),
            ("expected_accuracy", spec.expected_accuracy),
            ("result_expansion", spec.result_expansion),
        ):
            if value is not None:
                values[name] = value
        if values:
            raw_index.set(**values)
        return
    values = {}
    for name, value in (
        ("epsilon", spec.epsilon),
        ("blob_epsilon", spec.blob_epsilon),
        ("result_expansion", spec.result_expansion),
        ("radius", spec.search_radius),
        ("edge_size", spec.search_edge_size),
        ("exploration_size", spec.exploration_size),
        ("exact_result_expansion", spec.exact_result_expansion),
        ("num_of_probes", spec.num_of_probes),
    ):
        if value is not None:
            values[name] = value
    if values:
        raw_index.set(**values)


def _cleanup_build(build: BuildArtifact, temporary_parent: Path) -> None:
    """Release only a directory created for this runner's external scratch."""

    try:
        close = getattr(build.index, "close", None)
        if close is not None:
            close()
    finally:
        _remove_owned_directory(build.temporary_root, temporary_parent)


def _run_id(
    cache: CacheInputs,
    spec: NgtSpec,
    parameters: Mapping[str, Any],
    hardware_id: str,
    image: str,
) -> str:
    return stable_id(
        "ngt-run",
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


def _result_path(output_dir: Path, label: str, run_id: str, repeat: int) -> Path:
    slug = re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")
    return output_dir / f"{slug}-{run_id}-repeat-{repeat}.json"


def _is_terminal_result(path: Path, run_id: str) -> bool:
    if not path.is_file():
        return False
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    return record.get("run_id") == run_id and not validate_result_record(record)


def _base_provenance(
    cache: CacheInputs,
    *,
    repository_root: Path,
    hardware: Mapping[str, Any],
    image: str,
    container: str,
    ngt_version: str | None,
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
        "ngt_version": ngt_version,
        "official_ngt_sources": list(OFFICIAL_NGT_SOURCES),
        "benchmark_build_method": "pipeline-helper-serialize-reload-v1",
        "production_build_difference": (
            "The benchmark supplies the production helper with the flattened float32 memmap, "
            "whereas BUILD_NGT_INDEX first concatenates all window shards. The runner does use "
            "the production build, serialize, and load code paths, but its build RSS is not a "
            "direct pipeline peak-RSS requirement."
        ),
        "production_build_rss_comparable": False,
        "metric_note": (
            "All cache vectors are L2-normalized. NGT/QG use NGT cosine; QBG uses the "
            "pipeline's L2 path, whose nearest-neighbour ranking is cosine-equivalent on unit vectors."
        ),
    }
    return provenance


def _unavailable_record(
    cache: CacheInputs,
    spec: NgtSpec,
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
        measurement["capacity_preflight"] = asdict(capacity)
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


def _measurement_payload(measured: Any, build: BuildArtifact, spec: NgtSpec) -> dict[str, Any]:
    payload = dict(getattr(measured, "measurement", {}) or {})
    payload.update(
        {
            "timing_scope": "warm in-process NGT search; index build and serialization excluded from QPS",
            "index_constructor": (
                "bin/ngt_index.py::build_populated_index"
                if spec.mode == "pipeline"
                else "benchmarks/run_ngt.py native NGT API/CLI experimental builder"
            ),
            "pipeline_faithful": spec.mode == "pipeline",
            "ngt_native_metric": build.details.get("metric"),
            "normalized_l2_cosine_equivalent": build.details.get("metric") == "l2",
            "serialized_index_bytes_method": "recursive byte sum of NGT index directory",
            "build_peak_rss_bytes": build.peak_rss_bytes,
            "build_peak_vram_bytes": build.peak_vram_bytes,
            "serialize_seconds": build.serialize_seconds,
            "load_seconds": build.load_seconds,
            "production_build_rss_comparable": False,
            "qbg_child_max_rss_bytes_cumulative": build.child_max_rss_bytes,
        }
    )
    return payload


def _max_optional(*values: int | None) -> int | None:
    present = [int(value) for value in values if value is not None]
    return max(present) if present else None


def _only_filter(specs: Sequence[NgtSpec], requested: str | None) -> list[NgtSpec]:
    if not requested:
        return list(specs)
    requested_types = {item.strip().upper() for item in requested.split(",") if item.strip()}
    unknown = sorted(requested_types.difference({"NGT", "QG", "QBG"}))
    if unknown:
        raise ValueError(f"unknown --only values: {', '.join(unknown)}")
    return [spec for spec in specs if spec.kind in requested_types]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    parser.add_argument("--dataset", required=True, help="prepared cache dataset id")
    parser.add_argument("--output-dir", type=Path, help="defaults to <cache>/<dataset>/results/ngt")
    parser.add_argument("--plan", choices=("frontier", "smoke"), default="frontier")
    parser.add_argument("--only", help="comma-separated NGT structures: ngt,qg,qbg")
    parser.add_argument(
        "--include-native-frontier",
        action="store_true",
        help="also benchmark documented NGT controls not exposed by GINflow; records are experimental",
    )
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--warmup-queries", type=int, default=DEFAULT_WARMUP_QUERIES)
    parser.add_argument("--query-batch-size", type=int, default=DEFAULT_QUERY_BATCH_SIZE)
    parser.add_argument("--max-ram-fraction", type=float, default=DEFAULT_MAX_RAM_FRACTION)
    parser.add_argument("--max-disk-fraction", type=float, default=DEFAULT_MAX_DISK_FRACTION)
    parser.add_argument("--image", default=NGT_IMAGE, help="pinned runtime image recorded in provenance")
    parser.add_argument("--container", default="docker", help="runtime type recorded in provenance")
    parser.add_argument("--scratch-dir", type=Path, help="external scratch root; defaults below the dataset cache")
    parser.add_argument("--no-resume", action="store_true", help="re-measure terminal result records")
    parser.add_argument("--dry-run", action="store_true", help="validate cache and print capacity plan without importing NGT")
    return parser.parse_args(argv)


def _validate_args(args: argparse.Namespace) -> None:
    for name in ("repeats", "warmup_queries", "query_batch_size"):
        if int(getattr(args, name)) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 1")
    for name in ("max_ram_fraction", "max_disk_fraction"):
        value = float(getattr(args, name))
        if not 0.0 < value <= 1.0:
            raise ValueError(f"--{name.replace('_', '-')} must be in (0, 1]")
    # Overnight results on the user-requested corpora must have >=3 repeats.
    if args.dataset in {"rouskin_6k", "rouskin_30k"} and int(args.repeats) < 3:
        raise ValueError("rouskin benchmark datasets require --repeats >= 3")


def run(args: argparse.Namespace) -> int:
    _validate_args(args)
    repository_root = Path(__file__).resolve().parents[1]
    cache_dir = cache_root_from_args(args.cache_dir)
    cache = load_cache_inputs(cache_dir, args.dataset)
    if cache.queries.shape[0] % args.query_batch_size:
        raise ValueError("prepared query count must be an exact multiple of --query-batch-size")
    if args.warmup_queries > cache.queries.shape[0]:
        raise ValueError("--warmup-queries must not exceed the prepared query count")
    if args.warmup_queries % args.query_batch_size:
        raise ValueError("--warmup-queries must be an exact multiple of --query-batch-size")

    output_dir = (args.output_dir or (cache.root / "results" / BACKEND)).expanduser().resolve()
    temporary_parent = (args.scratch_dir or (cache.root / "scratch" / BACKEND)).expanduser().resolve()
    hardware = hardware_snapshot()
    specs = _only_filter(
        frontier_specs(
            cache.dataset_window_count,
            cache.dimension,
            plan=args.plan,
            include_native_frontier=args.include_native_frontier,
        ),
        args.only,
    )
    if not specs:
        raise ValueError("the selected --only filter produced no benchmark settings")
    available_ram = _available_ram_bytes()
    free_disk = _free_disk_bytes(temporary_parent.parent if temporary_parent.parent.exists() else cache.root)
    capacities = {
        spec: capacity_check(
            spec,
            n_vectors=cache.dataset_window_count,
            dimension=cache.dimension,
            available_ram_bytes=available_ram,
            free_disk_bytes=free_disk,
            max_ram_fraction=args.max_ram_fraction,
            max_disk_fraction=args.max_disk_fraction,
        )
        for spec in specs
    }
    if args.dry_run:
        print(
            json.dumps(
                {
                    "schema_version": "ginflow-ngt-benchmark-plan-v1",
                    "dataset_id": cache.dataset_id,
                    "dataset_window_count": cache.dataset_window_count,
                    "dimension": cache.dimension,
                    "image": args.image,
                    "hardware_id": hardware["id"],
                    "settings": [
                        {
                            "label": spec.label(),
                            "spec": asdict(spec),
                            "capacity": asdict(capacities[spec]),
                        }
                        for spec in specs
                    ],
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    ngtpy = _import_ngtpy()
    from importlib.metadata import version

    adapter = _pipeline_ngt_adapter()
    base_provenance = _base_provenance(
        cache,
        repository_root=repository_root,
        hardware=hardware,
        image=args.image,
        container=args.container,
        ngt_version=str(version("ngt")),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for build_spec, search_specs in group_build_specs(specs):
        active: list[tuple[NgtSpec, dict[str, Any], str, list[Path]]] = []
        for spec in search_specs:
            parameters = spec.parameters(query_batch_size=args.query_batch_size)
            run_id = _run_id(cache, spec, parameters, str(hardware["id"]), args.image)
            paths = [_result_path(output_dir, spec.label(), run_id, repeat) for repeat in range(args.repeats)]
            terminal_paths = paths if capacities[spec].feasible else [paths[0]]
            if not args.no_resume and all(_is_terminal_result(path, run_id) for path in terminal_paths):
                continue
            capacity = capacities[spec]
            if not capacity.feasible:
                write_result(
                    paths[0],
                    _unavailable_record(
                        cache,
                        spec,
                        parameters=parameters,
                        run_id=run_id,
                        provenance=base_provenance,
                        status="skipped",
                        reason=str(capacity.reason),
                        capacity=capacity,
                    ),
                )
                continue
            active.append((spec, parameters, run_id, paths))
        if not active:
            continue
        try:
            build = _build_index(
                ngtpy,
                adapter,
                cache,
                build_spec,
                temporary_parent=temporary_parent,
            )
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            for spec, parameters, run_id, paths in active:
                write_result(
                    paths[0],
                    _unavailable_record(
                        cache,
                        spec,
                        parameters=parameters,
                        run_id=run_id,
                        provenance=base_provenance,
                        status="error",
                        reason=reason,
                        capacity=None,
                    ),
                )
            continue
        try:
            for spec, parameters, run_id, paths in active:
                try:
                    _native_search_set(build.index, spec)
                    for measured in measure_search_repeats(
                        build.index.search,
                        cache.queries,
                        cache.ground_truth_ids,
                        k=K,
                        warmup_queries=args.warmup_queries,
                        repeats=args.repeats,
                        query_batch_size=args.query_batch_size,
                    ):
                        repeat = int(measured.repeat)
                        if not args.no_resume and _is_terminal_result(paths[repeat], run_id):
                            continue
                        write_result(
                            paths[repeat],
                            make_result_record(
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
                                provenance={
                                    **base_provenance,
                                    "pipeline_faithful": spec.mode == "pipeline",
                                    "ngt_index_type": spec.kind,
                                },
                                qps=float(measured.qps),
                                latency_ms=measured.latency_ms,
                                recall=float(measured.recall_at_100),
                                index_bytes=build.index_bytes,
                                build_seconds=build.build_seconds,
                                peak_rss_bytes=_max_optional(build.peak_rss_bytes, measured.peak_rss_bytes),
                                peak_vram_bytes=_max_optional(build.peak_vram_bytes, measured.peak_vram_bytes),
                                measurement=_measurement_payload(measured, build, spec),
                            ),
                        )
                except Exception as exc:
                    reason = f"{type(exc).__name__}: {exc}"
                    write_result(
                        paths[0],
                        _unavailable_record(
                            cache,
                            spec,
                            parameters=parameters,
                            run_id=run_id,
                            provenance=base_provenance,
                            status="error",
                            reason=reason,
                            capacity=None,
                        ),
                    )
        finally:
            _cleanup_build(build, temporary_parent)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
