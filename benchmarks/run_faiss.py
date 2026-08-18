#!/usr/bin/env python3
"""Benchmark pinned FAISS index frontiers against a prepared GINflow cache.

The runner deliberately builds from the flattened ``vectors.npy`` memmap in
bounded add batches.  That makes the benchmark usable for the 30k corpus on a
workstation where materialising a second full database array would be unsafe.
It records every measured repeat with the backend-neutral
``ginflow-benchmark-v1`` contract and records capacity refusals as explicit
``skipped`` rows instead of risking an out-of-memory termination.

Run this script inside the repository's pinned FAISS CPU or GPU image.  For
example, the CPU invocation is:

```
docker run --rm -v "$PWD":/work -w /work \\
  -v /mnt/ssd_samsung/ginflow-benchmarks:/bench-cache \\
  community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e \\
  python3 benchmarks/run_faiss.py --cache-dir /bench-cache --dataset rouskin_6k
```
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import os
import re
import subprocess
import sys
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
    sha256_array,
    software_snapshot,
    stable_id,
    validate_result_record,
    write_result,
)


FAISS_CPU_IMAGE = "community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e"
FAISS_GPU_IMAGE = "community.wave.seqera.io/library/python_numpy_faiss-gpu:7382ed4d7c6e6258"
BACKEND = "faiss"
K = 100
DEFAULT_REPEATS = 3
DEFAULT_WARMUP_QUERIES = 64
DEFAULT_QUERY_BATCH_SIZE = 32
DEFAULT_ADD_BATCH_SIZE = 16_384
DEFAULT_RESOURCE_TEMP_MIB = 256
DEFAULT_MAX_RAM_FRACTION = 0.80
DEFAULT_MAX_VRAM_FRACTION = 0.75


@dataclass(frozen=True)
class CacheInputs:
    """Validated backend-neutral cache inputs needed for a FAISS run."""

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
    """One build and search setting; search fields are part of the frontier."""

    index_type: str
    device: str = "cpu"
    nlist: int | None = None
    nprobe: int | None = None
    pq_m: int | None = None
    pq_nbits: int | None = None
    hnsw_m: int | None = None
    hnsw_ef_construction: int | None = None
    hnsw_ef_search: int | None = None

    def build_key(self) -> tuple[Any, ...]:
        return (
            self.index_type,
            self.device,
            self.nlist,
            self.pq_m,
            self.pq_nbits,
            self.hnsw_m,
            self.hnsw_ef_construction,
        )

    def build_parameters(self) -> dict[str, Any]:
        return {
            "index_type": self.index_type,
            "device": self.device,
            "nlist": self.nlist,
            "pq_m": self.pq_m,
            "pq_nbits": self.pq_nbits,
            "hnsw_m": self.hnsw_m,
            "hnsw_ef_construction": self.hnsw_ef_construction,
        }

    def parameters(
        self,
        *,
        threads: int,
        add_batch_size: int,
        train_sample_size: int,
        query_batch_size: int,
        resource_temp_mib: int,
    ) -> dict[str, Any]:
        values = self.build_parameters()
        values.update(
            {
                "nprobe": self.nprobe,
                "hnsw_ef_search": self.hnsw_ef_search,
                "faiss_threads": int(threads),
                "add_batch_size": int(add_batch_size),
                "train_sample_size": int(train_sample_size),
                "query_batch_size": int(query_batch_size),
                "gpu_resource_temp_mib": int(resource_temp_mib) if self.device == "gpu" else None,
            }
        )
        return values

    def label(self) -> str:
        pieces = [self.device, self.index_type.lower()]
        if self.nlist is not None:
            pieces.append(f"nlist{self.nlist}")
        if self.pq_m is not None:
            pieces.append(f"m{self.pq_m}")
        if self.pq_nbits is not None:
            pieces.append(f"b{self.pq_nbits}")
        if self.hnsw_m is not None:
            pieces.append(f"m{self.hnsw_m}")
        if self.nprobe is not None:
            pieces.append(f"nprobe{self.nprobe}")
        if self.hnsw_ef_search is not None:
            pieces.append(f"ef{self.hnsw_ef_search}")
        return "-".join(pieces)


@dataclass
class BuildArtifact:
    """An in-memory index plus build measurements kept alive for its searches."""

    index: Any
    build_seconds: float
    serialize_seconds: float
    load_seconds: float
    index_bytes: int
    peak_rss_bytes: int | None
    peak_vram_bytes: int | None
    train_sample_size: int
    resources: Any | None = None


@dataclass(frozen=True)
class CapacityCheck:
    feasible: bool
    estimated_bytes: int
    available_bytes: int | None
    safe_limit_bytes: int | None
    reason: str | None = None


def cache_root_from_args(value: str | None) -> Path:
    """Resolve the shared external cache without ever creating it implicitly."""

    configured = value or os.environ.get("GINFLOW_BENCHMARK_CACHE")
    return Path(configured).expanduser().resolve() if configured else (Path.cwd() / ".benchmark-cache").resolve()


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _required_text(mapping: Mapping[str, Any], key: str, source: Path) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{source} is missing a non-empty {key!r}")
    return value


def load_cache_inputs(cache_dir: Path, dataset_id: str) -> CacheInputs:
    """Load and cross-check a prepared vector/query/exact-neighbour cache."""

    root = cache_dir / dataset_id
    flat_dir = root / "flat"
    query_dir = root / "queries"
    truth_dir = root / "ground-truth"
    flatten = _read_json(flat_dir / "flatten-manifest.json")
    selection = _read_json(query_dir / "query-selection.json")
    truth = _read_json(truth_dir / "ground-truth.json")
    vectors_path = flat_dir / str(flatten.get("vectors", {}).get("path") or "vectors.npy")
    query_path = query_dir / str(selection.get("query_vectors", {}).get("path") or "queries.npy")
    ground_truth_path = truth_dir / str(truth.get("ground_truth", {}).get("path") or "ground-truth.npz")
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
    if queries.ndim != 2 or queries.shape[0] < 1 or queries.shape[1] != vectors.shape[1]:
        raise ValueError("prepared query vectors do not match the flattened vector dimension")
    if ground_truth_ids.shape != (queries.shape[0], K):
        raise ValueError(
            f"exact ground truth shape {ground_truth_ids.shape} must be ({queries.shape[0]}, {K})"
        )
    if np.any(ground_truth_ids < 0) or np.any(ground_truth_ids >= vectors.shape[0]):
        raise ValueError("exact ground truth contains invalid vector identifiers")
    if not np.isfinite(queries).all() or not np.isfinite(vectors[: min(vectors.shape[0], 1024)]).all():
        raise ValueError("prepared cache contains non-finite vectors")

    embedding_cache_id = _required_text(flatten, "embedding_cache_id", flat_dir / "flatten-manifest.json")
    query_selection_id = _required_text(selection, "query_selection_id", query_dir / "query-selection.json")
    ground_truth_cache_id = _required_text(truth, "ground_truth_cache_id", truth_dir / "ground-truth.json")
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

    query_ids_sha256 = _required_text(selection, "query_ids_sha256", query_dir / "query-selection.json")
    if truth.get("query_ids_sha256") != query_ids_sha256:
        raise ValueError("exact ground truth and query selection have different query id digests")
    ground_truth_ids_sha256 = _required_text(
        truth.get("ground_truth", {}), "ids_sha256", truth_dir / "ground-truth.json"
    )
    if sha256_array(ground_truth_ids) != ground_truth_ids_sha256:
        raise ValueError("exact ground truth id digest does not match its manifest")
    if str(flatten.get("dataset_id")) != dataset_id or str(selection.get("dataset_id")) != dataset_id:
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


def _pq_m(dimension: int, preferred: int) -> int:
    for candidate in (preferred, 64, 32, 16, 8, 4, 2, 1):
        if candidate <= dimension and dimension % candidate == 0:
            return candidate
    raise ValueError(f"no product-quantizer subvector count divides dimension={dimension}")


def _frontier_nlists(n_vectors: int, *, plan: str, device: str) -> list[int]:
    if plan == "smoke":
        return _unique_positive((min(16, max(1, int(math.sqrt(n_vectors)))),), n_vectors)
    if device == "gpu":
        # The GPU frontier intentionally avoids a second large GPU build on an
        # 8 GiB host.  Its capacity check will still refuse unsafe data sizes.
        return _unique_positive((2048 if n_vectors < 1_500_000 else 4096,), n_vectors)
    candidates = (2048, 4096) if n_vectors < 1_500_000 else (4096, 8192)
    return _unique_positive(candidates, n_vectors)


def _probe_ladder(nlist: int, *, plan: str, n_vectors: int) -> list[int]:
    if plan == "smoke":
        return _unique_positive((1, min(2, nlist)), nlist)
    values = (8, 16, 32, 64, 128) if n_vectors < 1_500_000 else (16, 32, 64, 128, 256)
    return _unique_positive(values, nlist)


def frontier_specs(n_vectors: int, dimension: int, *, plan: str, device: str) -> list[IndexSpec]:
    """Return a compact, reproducible FAISS accuracy/throughput frontier."""

    if plan not in {"frontier", "smoke"}:
        raise ValueError("plan must be frontier or smoke")
    if device not in {"cpu", "gpu"}:
        raise ValueError("device must be cpu or gpu")
    nlists = _frontier_nlists(n_vectors, plan=plan, device=device)
    primary_nlist = nlists[0]
    probes = _probe_ladder(primary_nlist, plan=plan, n_vectors=n_vectors)
    specs: list[IndexSpec] = [IndexSpec("FlatIP", device=device)]
    for nlist in nlists:
        selected_probes = _probe_ladder(nlist, plan=plan, n_vectors=n_vectors)
        if plan == "frontier" and nlist != primary_nlist:
            # A second list count is a build sensitivity, not a duplicated
            # exhaustive ladder.
            selected_probes = _unique_positive((selected_probes[1], selected_probes[-2], selected_probes[-1]), nlist)
        specs.extend(IndexSpec("IVFFlat", device=device, nlist=nlist, nprobe=nprobe) for nprobe in selected_probes)

    pq_m = _pq_m(dimension, 32)
    if plan == "smoke":
        pq_m = _pq_m(dimension, min(8, dimension))
    # Classic FAISS GPU supports IVFPQ only with 8-bit codes.  The small CPU
    # smoke uses 4-bit codes to keep its training fixture compact.
    pq_nbits = 8 if device == "gpu" or plan == "frontier" else 4
    specs.extend(
        IndexSpec("IVFPQ", device=device, nlist=primary_nlist, nprobe=nprobe, pq_m=pq_m, pq_nbits=pq_nbits)
        for nprobe in probes
    )
    if plan == "frontier" and device == "cpu":
        compressed_m = _pq_m(dimension, 16)
        compression_probes = probes[max(1, len(probes) // 2) :]
        specs.extend(
            IndexSpec("IVFPQ", device=device, nlist=primary_nlist, nprobe=nprobe, pq_m=compressed_m, pq_nbits=8)
            for nprobe in compression_probes
        )
        low_bits_probes = probes[-2:]
        specs.extend(
            IndexSpec("IVFPQ", device=device, nlist=primary_nlist, nprobe=nprobe, pq_m=pq_m, pq_nbits=4)
            for nprobe in low_bits_probes
        )
    if device == "cpu":
        ef_values = (8, 16, 32) if plan == "smoke" else ((16, 32, 64, 128) if n_vectors < 1_500_000 else (32, 64, 128, 256))
        specs.extend(
            IndexSpec(
                "HNSW",
                device="cpu",
                hnsw_m=16 if plan == "smoke" else 32,
                hnsw_ef_construction=40 if plan == "smoke" else 80,
                hnsw_ef_search=int(value),
            )
            for value in ef_values
        )
    return specs


def group_build_specs(specs: Iterable[IndexSpec]) -> Iterator[tuple[IndexSpec, list[IndexSpec]]]:
    grouped: dict[tuple[Any, ...], list[IndexSpec]] = defaultdict(list)
    for spec in specs:
        grouped[spec.build_key()].append(spec)
    for values in grouped.values():
        yield values[0], values


def _requires_training(spec: IndexSpec) -> bool:
    return spec.index_type in {"IVFFlat", "IVFPQ"}


def train_sample_size(spec: IndexSpec, n_vectors: int, requested: int | None) -> int:
    if not _requires_training(spec):
        return 0
    if requested is not None:
        if requested < 1:
            raise ValueError("--train-sample-size must be >= 1")
        return min(int(requested), n_vectors)
    nlist = int(spec.nlist or 1)
    nbits = int(spec.pq_nbits or 0)
    # Forty examples per centroid/codeword is an intentionally modest but
    # defensible training floor for an overnight frontier, capped by the data.
    target = max(65_536, nlist * 40, (2**nbits) * 40 if nbits else 0)
    return min(target, n_vectors)


def _raw_vector_bytes(n_vectors: int, dimension: int) -> int:
    return int(n_vectors) * int(dimension) * np.dtype(np.float32).itemsize


def final_index_bytes_estimate(spec: IndexSpec, n_vectors: int, dimension: int) -> int:
    """Conservative final-storage estimate used only for capacity gating."""

    raw = _raw_vector_bytes(n_vectors, dimension)
    if spec.index_type == "FlatIP":
        return raw
    if spec.index_type == "IVFFlat":
        return raw + (n_vectors * 8) + ((spec.nlist or 1) * dimension * 4)
    if spec.index_type == "IVFPQ":
        code_bytes = math.ceil(int(spec.pq_m or 1) * int(spec.pq_nbits or 8) / 8)
        codebook_bytes = (2 ** int(spec.pq_nbits or 8)) * dimension * 4
        return (n_vectors * (code_bytes + 8)) + ((spec.nlist or 1) * dimension * 4) + codebook_bytes
    if spec.index_type == "HNSW":
        # FAISS documents 2 * 4 * M links on level zero; use 1.25x to leave
        # room for upper levels and allocator overhead.
        links = math.ceil(n_vectors * int(spec.hnsw_m or 32) * 2 * 4 * 1.25)
        return raw + links
    raise ValueError(f"unsupported FAISS index type {spec.index_type!r}")


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
    except (ValueError, IndexError):
        return None
    return max(0, int((total_mib - used_mib) * 1024 * 1024))


def cpu_capacity_check(
    spec: IndexSpec,
    *,
    n_vectors: int,
    dimension: int,
    training_rows: int,
    add_batch_size: int,
    available_bytes: int | None,
    max_fraction: float,
) -> CapacityCheck:
    final_bytes = final_index_bytes_estimate(spec, n_vectors, dimension)
    train_bytes = training_rows * dimension * 4
    batch_bytes = min(add_batch_size, n_vectors) * dimension * 4
    # Streaming avoids a full input duplicate. FAISS training nevertheless
    # allocates working buffers, and the production-equivalent artifact
    # roundtrip briefly keeps a source and a freshly reloaded CPU index alive.
    estimated = (2 * final_bytes) + (2 * train_bytes) + batch_bytes + (512 * 1024**2)
    if available_bytes is None:
        return CapacityCheck(True, estimated, None, None)
    limit = int(available_bytes * max_fraction)
    if estimated > limit:
        return CapacityCheck(
            False,
            estimated,
            available_bytes,
            limit,
            "estimated FAISS CPU peak exceeds the configured safe MemAvailable fraction",
        )
    return CapacityCheck(True, estimated, available_bytes, limit)


def gpu_capacity_check(
    spec: IndexSpec,
    *,
    n_vectors: int,
    dimension: int,
    training_rows: int,
    add_batch_size: int,
    resource_temp_mib: int,
    available_bytes: int | None,
    max_fraction: float,
) -> CapacityCheck:
    if spec.index_type not in {"FlatIP", "IVFFlat", "IVFPQ"}:
        return CapacityCheck(False, 0, available_bytes, None, "this FAISS index type is CPU-only")
    final_bytes = final_index_bytes_estimate(spec, n_vectors, dimension)
    train_bytes = training_rows * dimension * 4
    batch_bytes = min(add_batch_size, n_vectors) * dimension * 4
    # GPU IVF/PQ construction uses the final device index plus a train/add
    # buffer.  The 15% headroom covers allocator alignment and library work.
    estimated = int((final_bytes + train_bytes + batch_bytes + resource_temp_mib * 1024**2) * 1.15)
    if available_bytes is None:
        return CapacityCheck(False, estimated, None, None, "no visible NVIDIA GPU memory information")
    limit = int(available_bytes * max_fraction)
    if estimated > limit:
        return CapacityCheck(
            False,
            estimated,
            available_bytes,
            limit,
            "estimated FAISS GPU peak exceeds the configured free-VRAM fraction",
        )
    return CapacityCheck(True, estimated, available_bytes, limit)


def _import_faiss() -> Any:
    try:
        import faiss
    except ImportError as exc:  # pragma: no cover - tested by invoking the pinned image
        raise RuntimeError(
            "FAISS is not installed. Run this script in the pinned FAISS CPU or GPU container."
        ) from exc
    version = str(getattr(faiss, "__version__", ""))
    if not version.startswith("1.10"):
        raise RuntimeError(f"this benchmark requires pinned FAISS 1.10.x, found {version or 'unknown'}")
    return faiss


def _pipeline_faiss_adapter() -> Any:
    """Import GINflow's FAISS constructor shim without requiring a package."""

    bin_dir = Path(__file__).resolve().parents[1] / "bin"
    if str(bin_dir) not in sys.path:
        sys.path.insert(0, str(bin_dir))
    module = importlib.import_module("faiss_index")
    if getattr(module, "faiss", None) is None:
        raise RuntimeError("GINflow's faiss_index adapter could not import FAISS")
    return module


def _make_cpu_index(spec: IndexSpec, dimension: int, n_vectors: int) -> Any:
    """Reuse the production constructor and parameter-validation semantics."""

    adapter = _pipeline_faiss_adapter()
    options = adapter.IndexOptions(
        index_type=spec.index_type,
        nlist=spec.nlist,
        nprobe=spec.nprobe,
        pq_m=int(spec.pq_m or 16),
        pq_nbits=int(spec.pq_nbits or 8),
        hnsw_m=int(spec.hnsw_m or 32),
        hnsw_ef_construction=int(spec.hnsw_ef_construction or 40),
        hnsw_ef_search=int(spec.hnsw_ef_search or 16),
        gpu=False,
    )
    index, _ = adapter.make_cpu_index(dimension, n_vectors, options)
    return index


def _training_sample(vectors: np.ndarray, count: int) -> np.ndarray:
    if count < 1:
        raise ValueError("training sample count must be positive")
    if count >= vectors.shape[0]:
        return np.ascontiguousarray(vectors, dtype=np.float32)
    identifiers = np.linspace(0, vectors.shape[0] - 1, num=count, dtype=np.int64)
    return np.ascontiguousarray(vectors[identifiers], dtype=np.float32)


def _add_stream(index: Any, vectors: np.ndarray, *, batch_size: int) -> None:
    if batch_size < 1:
        raise ValueError("--add-batch-size must be >= 1")
    for start in range(0, vectors.shape[0], batch_size):
        stop = min(vectors.shape[0], start + batch_size)
        block = np.ascontiguousarray(vectors[start:stop], dtype=np.float32)
        index.add(block)
        del block


def _serialize_and_reload(
    faiss: Any,
    index: Any,
    *,
    scratch_dir: Path,
    build_id: str,
    gpu: bool,
    resources: Any | None,
    gpu_device: int,
) -> tuple[Any, int, float, float]:
    """Exercise the same FAISS artifact boundary used by BUILD/SEARCH.

    ``BUILD_FAISS_INDEX`` writes a CPU ``index.faiss`` and ``SEARCH_FAISS``
    reads it before optionally moving it to the GPU.  Keep the temporary file
    on the external cache volume and remove only that runner-owned file after
    loading it, so QPS comes from a warm reloaded index rather than the build
    object that happens to remain in memory.
    """

    scratch_dir.mkdir(parents=True, exist_ok=True)
    path = scratch_dir / f".{build_id}-{os.getpid()}-{time.time_ns()}.faiss"
    cpu_index: Any | None = None
    loaded_cpu: Any | None = None
    loaded_index: Any | None = None
    try:
        cpu_index = faiss.index_gpu_to_cpu(index) if gpu else index
        serialize_started = time.perf_counter()
        faiss.write_index(cpu_index, str(path))
        serialize_seconds = time.perf_counter() - serialize_started
        index_bytes = int(path.stat().st_size)
        load_started = time.perf_counter()
        loaded_cpu = faiss.read_index(str(path))
        if gpu:
            if resources is None:
                raise RuntimeError("FAISS GPU resources were not retained for the search reload")
            loaded_index = faiss.index_cpu_to_gpu(resources, gpu_device, loaded_cpu)
        else:
            loaded_index = loaded_cpu
        load_seconds = time.perf_counter() - load_started
        return loaded_index, index_bytes, float(serialize_seconds), float(load_seconds)
    except Exception:
        if loaded_index is not None:
            del loaded_index
        raise
    finally:
        if gpu and loaded_cpu is not None:
            del loaded_cpu
        if gpu and cpu_index is not None:
            del cpu_index
        if path.exists():
            path.unlink()


def build_index(
    faiss: Any,
    cache: CacheInputs,
    spec: IndexSpec,
    *,
    training_rows: int,
    add_batch_size: int,
    resource_temp_mib: int,
    gpu_device: int,
    scratch_dir: Path,
    build_id: str,
) -> BuildArtifact:
    """Build a CPU/GPU FAISS index from a memmap in bounded vector batches."""

    resources: Any | None = None
    gpu = spec.device == "gpu"
    with PeakResourceSampler(gpu_device=gpu_device if gpu else None) as sampler:
        cpu_index = _make_cpu_index(spec, cache.dimension, cache.dataset_window_count)
        if gpu:
            if not hasattr(faiss, "StandardGpuResources") or int(faiss.get_num_gpus()) <= gpu_device:
                raise RuntimeError("FAISS GPU runtime or requested CUDA device is unavailable")
            resources = faiss.StandardGpuResources()
            if hasattr(resources, "setTempMemory"):
                resources.setTempMemory(int(resource_temp_mib) * 1024 * 1024)
            index = faiss.index_cpu_to_gpu(resources, gpu_device, cpu_index)
            del cpu_index
        else:
            index = cpu_index
        started = time.perf_counter()
        if _requires_training(spec):
            sample = _training_sample(cache.vectors, training_rows)
            try:
                index.train(sample)
            finally:
                del sample
        _add_stream(index, cache.vectors, batch_size=add_batch_size)
        build_seconds = time.perf_counter() - started
        if int(index.ntotal) != cache.dataset_window_count:
            raise RuntimeError(f"FAISS added {index.ntotal} vectors, expected {cache.dataset_window_count}")
        reloaded_index, index_bytes, serialize_seconds, load_seconds = _serialize_and_reload(
            faiss,
            index,
            scratch_dir=scratch_dir,
            build_id=build_id,
            gpu=gpu,
            resources=resources,
            gpu_device=gpu_device,
        )
        del index
        index = reloaded_index
    return BuildArtifact(
        index=index,
        build_seconds=build_seconds,
        serialize_seconds=serialize_seconds,
        load_seconds=load_seconds,
        index_bytes=index_bytes,
        peak_rss_bytes=sampler.peak_rss_bytes,
        peak_vram_bytes=sampler.peak_vram_bytes,
        train_sample_size=training_rows,
        resources=resources,
    )


def _set_search_parameters(faiss: Any, index: Any, spec: IndexSpec) -> None:
    _ = faiss
    adapter = _pipeline_faiss_adapter()
    adapter.apply_search_params(index, nprobe=spec.nprobe, ef_search=spec.hnsw_ef_search)


def _faiss_metadata(faiss: Any) -> dict[str, Any]:
    compile_options = getattr(faiss, "get_compile_options", None)
    return {
        "faiss_version": str(getattr(faiss, "__version__", "unknown")),
        "faiss_compile_options": str(compile_options()) if compile_options is not None else None,
    }


def _run_id(cache: CacheInputs, spec: IndexSpec, parameters: Mapping[str, Any], hardware_id: str, image: str) -> str:
    return stable_id(
        "faiss-run",
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
    faiss: Any | None,
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
        "benchmark_build_method": "memmap-streaming-train-add-serialize-reload-v1",
        "production_build_difference": (
            "This benchmark streams normalized vectors from flat/vectors.npy and trains IVF/PQ "
            "on a bounded deterministic sample; GINflow build_faiss.py concatenates full shard "
            "arrays and trains on that full matrix. The runner does exercise the production FAISS "
            "write/read (and optional GPU reload) boundary, but build/RSS values remain harness "
            "measurements rather than production pipeline capacity requirements."
        ),
    }
    if faiss is not None:
        provenance.update(_faiss_metadata(faiss))
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


def _measurement_payload(measurement: Any, *, build: BuildArtifact, query_batch_size: int) -> dict[str, Any]:
    """Keep shared harness metadata while documenting the FAISS build boundary."""

    payload = dict(getattr(measurement, "measurement", {}) or {})
    payload.update(
        {
            "timing_scope": "warm serialized-and-reloaded FAISS search; build, serialization, and load excluded from QPS",
            "query_batch_size": int(query_batch_size),
            "build_peak_rss_bytes": build.peak_rss_bytes,
            "build_peak_vram_bytes": build.peak_vram_bytes,
            "serialized_index_bytes_method": "faiss.write_index CPU serialization to temporary external-cache file",
            "serialize_seconds": build.serialize_seconds,
            "load_seconds": build.load_seconds,
            "index_constructor": "bin/faiss_index.py::make_cpu_index",
            "build_method": "memmap-streaming-train-add-v1",
            "training_sample_method": "deterministic evenly spaced rows from the flattened memmap",
            "production_build_rss_comparable": False,
        }
    )
    return payload


def _gpu_synchronize_callback(gpu_device: int) -> tuple[Any | None, str]:
    """Return a CUDA synchronization hook when the pinned image exposes one."""

    try:
        import torch

        if torch.cuda.is_available():
            return (lambda: torch.cuda.synchronize(gpu_device)), "torch.cuda.synchronize"
    except (ImportError, RuntimeError):
        pass
    try:
        import cupy

        return (lambda: cupy.cuda.Device(gpu_device).synchronize()), "cupy.cuda.Device.synchronize"
    except (ImportError, RuntimeError):
        return None, "not_available"


def _measure_spec(
    faiss: Any,
    cache: CacheInputs,
    spec: IndexSpec,
    build: BuildArtifact,
    *,
    repeats: int,
    warmup_queries: int,
    query_batch_size: int,
    gpu_device: int,
) -> Iterator[Any]:
    _set_search_parameters(faiss, build.index, spec)
    synchronize, _ = _gpu_synchronize_callback(gpu_device) if spec.device == "gpu" else (None, "not_applicable")
    return measure_search_repeats(
        build.index.search,
        cache.queries,
        cache.ground_truth_ids,
        k=K,
        warmup_queries=warmup_queries,
        repeats=repeats,
        gpu_device=gpu_device if spec.device == "gpu" else None,
        query_batch_size=query_batch_size,
        synchronize=synchronize,
    )


def _only_filter(specs: Sequence[IndexSpec], requested: str | None) -> list[IndexSpec]:
    if not requested:
        return list(specs)
    allowed = {item.strip().lower() for item in requested.split(",") if item.strip()}
    valid = {"flatip", "ivfflat", "ivfpq", "hnsw"}
    unknown = sorted(allowed.difference(valid))
    if unknown:
        raise ValueError(f"unknown --only values: {', '.join(unknown)}")
    return [spec for spec in specs if spec.index_type.lower() in allowed]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    parser.add_argument("--dataset", required=True, help="prepared cache dataset id")
    parser.add_argument("--output-dir", type=Path, help="result directory; defaults to <cache>/<dataset>/results/faiss")
    parser.add_argument("--device", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--gpu-device", type=int, default=0)
    parser.add_argument("--plan", choices=("frontier", "smoke"), default="frontier")
    parser.add_argument("--only", help="comma-separated index types: flatip,ivfflat,ivfpq,hnsw")
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--warmup-queries", type=int, default=DEFAULT_WARMUP_QUERIES)
    parser.add_argument("--query-batch-size", type=int, default=DEFAULT_QUERY_BATCH_SIZE)
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--add-batch-size", type=int, default=DEFAULT_ADD_BATCH_SIZE)
    parser.add_argument("--train-sample-size", type=int)
    parser.add_argument("--resource-temp-mib", type=int, default=DEFAULT_RESOURCE_TEMP_MIB)
    parser.add_argument("--max-ram-fraction", type=float, default=DEFAULT_MAX_RAM_FRACTION)
    parser.add_argument("--max-vram-fraction", type=float, default=DEFAULT_MAX_VRAM_FRACTION)
    parser.add_argument("--image", help="pinned runtime image recorded in provenance")
    parser.add_argument("--container", default="docker", help="runtime type recorded in provenance")
    parser.add_argument("--no-resume", action="store_true", help="re-measure settings whose complete result records exist")
    parser.add_argument("--dry-run", action="store_true", help="validate cache and print capacity plan without importing FAISS")
    return parser.parse_args(argv)


def _validate_args(args: argparse.Namespace) -> None:
    for name in ("repeats", "warmup_queries", "query_batch_size", "threads", "add_batch_size", "resource_temp_mib"):
        if int(getattr(args, name)) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 1")
    for name in ("max_ram_fraction", "max_vram_fraction"):
        value = float(getattr(args, name))
        if not 0.0 < value <= 1.0:
            raise ValueError(f"--{name.replace('_', '-')} must be in (0, 1]")
    if args.dataset in {"rouskin_6k", "rouskin_30k"} and args.repeats < DEFAULT_REPEATS:
        raise ValueError("rouskin benchmark datasets require --repeats >= 3")


def run(args: argparse.Namespace) -> int:
    _validate_args(args)
    repository_root = Path(__file__).resolve().parents[1]
    cache_dir = cache_root_from_args(args.cache_dir)
    cache = load_cache_inputs(cache_dir, args.dataset)
    if cache.queries.shape[0] % args.query_batch_size:
        raise ValueError(
            "prepared query count must be an exact multiple of --query-batch-size for fixed-batch timing"
        )
    if args.warmup_queries % args.query_batch_size:
        raise ValueError("--warmup-queries must be an exact multiple of --query-batch-size")
    image = args.image or (FAISS_GPU_IMAGE if args.device == "gpu" else FAISS_CPU_IMAGE)
    output_dir = (args.output_dir or (cache.root / "results" / BACKEND)).expanduser().resolve()
    hardware = hardware_snapshot()
    specs = _only_filter(
        frontier_specs(cache.dataset_window_count, cache.dimension, plan=args.plan, device=args.device), args.only
    )
    if not specs:
        raise ValueError("the selected --only filter produced no benchmark settings")
    base_provenance = _base_provenance(
        cache,
        repository_root=repository_root,
        hardware=hardware,
        image=image,
        container=args.container,
        faiss=None,
    )
    ram_available = _available_ram_bytes()
    vram_available = _available_vram_bytes(args.gpu_device) if args.device == "gpu" else None
    plan_rows: list[dict[str, Any]] = []
    for spec in specs:
        training_rows = train_sample_size(spec, cache.dataset_window_count, args.train_sample_size)
        capacity = (
            gpu_capacity_check(
                spec,
                n_vectors=cache.dataset_window_count,
                dimension=cache.dimension,
                training_rows=training_rows,
                add_batch_size=args.add_batch_size,
                resource_temp_mib=args.resource_temp_mib,
                available_bytes=vram_available,
                max_fraction=args.max_vram_fraction,
            )
            if spec.device == "gpu"
            else cpu_capacity_check(
                spec,
                n_vectors=cache.dataset_window_count,
                dimension=cache.dimension,
                training_rows=training_rows,
                add_batch_size=args.add_batch_size,
                available_bytes=ram_available,
                max_fraction=args.max_ram_fraction,
            )
        )
        plan_rows.append(
            {
                "label": spec.label(),
                "spec": asdict(spec),
                "train_sample_size": training_rows,
                "capacity": asdict(capacity),
            }
        )
    if args.dry_run:
        print(
            json.dumps(
                {
                    "schema_version": "ginflow-faiss-benchmark-plan-v1",
                    "dataset_id": cache.dataset_id,
                    "dataset_window_count": cache.dataset_window_count,
                    "dimension": cache.dimension,
                    "image": image,
                    "hardware_id": hardware["id"],
                    "settings": plan_rows,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    faiss = _import_faiss()
    faiss.omp_set_num_threads(int(args.threads))
    base_provenance = _base_provenance(
        cache,
        repository_root=repository_root,
        hardware=hardware,
        image=image,
        container=args.container,
        faiss=faiss,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir = cache.root / "scratch" / BACKEND
    capacity_by_spec = {tuple(row["spec"].items()): row["capacity"] for row in plan_rows}
    train_by_spec = {tuple(row["spec"].items()): int(row["train_sample_size"]) for row in plan_rows}
    for build_spec, search_specs in group_build_specs(specs):
        active_specs: list[tuple[IndexSpec, dict[str, Any], str, list[Path]]] = []
        for spec in search_specs:
            key = tuple(asdict(spec).items())
            training_rows = train_by_spec[key]
            parameters = spec.parameters(
                threads=args.threads,
                add_batch_size=args.add_batch_size,
                train_sample_size=training_rows,
                query_batch_size=args.query_batch_size,
                resource_temp_mib=args.resource_temp_mib,
            )
            run_id = _run_id(cache, spec, parameters, str(hardware["id"]), image)
            paths = [_result_path(output_dir, spec.label(), run_id, repeat) for repeat in range(args.repeats)]
            if not args.no_resume and all(_is_complete_result(path, run_id) for path in paths):
                continue
            capacity_payload = capacity_by_spec[key]
            capacity = CapacityCheck(**capacity_payload)
            if not capacity.feasible:
                record = _unavailable_record(
                    cache,
                    spec,
                    parameters=parameters,
                    run_id=run_id,
                    provenance=base_provenance,
                    status="skipped",
                    reason=str(capacity.reason),
                    capacity=capacity,
                )
                write_result(paths[0], record)
                continue
            active_specs.append((spec, parameters, run_id, paths))
        if not active_specs:
            continue

        build_key = stable_id(
            "faiss-build",
            {
                "dataset": cache.dataset_id,
                "embedding_cache_id": cache.embedding_cache_id,
                "build_spec": asdict(build_spec),
                "image": image,
            },
        )
        build_training_rows = train_by_spec[tuple(asdict(build_spec).items())]
        try:
            build = build_index(
                faiss,
                cache,
                build_spec,
                training_rows=build_training_rows,
                add_batch_size=args.add_batch_size,
                resource_temp_mib=args.resource_temp_mib,
                gpu_device=args.gpu_device,
                scratch_dir=scratch_dir,
                build_id=build_key,
            )
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            for spec, parameters, run_id, paths in active_specs:
                record = _unavailable_record(
                    cache,
                    spec,
                    parameters=parameters,
                    run_id=run_id,
                    provenance=base_provenance,
                    status="error",
                    reason=reason,
                    capacity=None,
                )
                write_result(paths[0], record)
            continue

        try:
            for spec, parameters, run_id, paths in active_specs:
                try:
                    samples = _measure_spec(
                        faiss,
                        cache,
                        spec,
                        build,
                        repeats=args.repeats,
                        warmup_queries=args.warmup_queries,
                        query_batch_size=args.query_batch_size,
                        gpu_device=args.gpu_device,
                    )
                    for measured in samples:
                        repeat = int(measured.repeat)
                        if not args.no_resume and _is_complete_result(paths[repeat], run_id):
                            continue
                        record = make_result_record(
                            backend=BACKEND,
                            dataset_id=cache.dataset_id,
                            dataset_window_count=cache.dataset_window_count,
                            dimension=cache.dimension,
                            parameter_label=spec.label(),
                            parameters=parameters,
                            run_id=run_id,
                            repeat=repeat,
                            warmup_queries=min(args.warmup_queries, cache.queries.shape[0]),
                            timed_queries=cache.queries.shape[0],
                            query_ids_sha256=cache.query_ids_sha256,
                            ground_truth_ids_sha256=cache.ground_truth_ids_sha256,
                            provenance=base_provenance,
                            qps=float(measured.qps),
                            latency_ms=measured.latency_ms,
                            recall=float(measured.recall_at_100),
                            index_bytes=build.index_bytes,
                            build_seconds=build.build_seconds,
                            peak_rss_bytes=max(
                                value for value in (build.peak_rss_bytes, measured.peak_rss_bytes) if value is not None
                            )
                            if build.peak_rss_bytes is not None or measured.peak_rss_bytes is not None
                            else None,
                            peak_vram_bytes=max(
                                value for value in (build.peak_vram_bytes, measured.peak_vram_bytes) if value is not None
                            )
                            if build.peak_vram_bytes is not None or measured.peak_vram_bytes is not None
                            else None,
                            measurement=_measurement_payload(
                                measured, build=build, query_batch_size=args.query_batch_size
                            ),
                        )
                        write_result(paths[repeat], record)
                except Exception as exc:
                    reason = f"{type(exc).__name__}: {exc}"
                    record = _unavailable_record(
                        cache,
                        spec,
                        parameters=parameters,
                        run_id=run_id,
                        provenance=base_provenance,
                        status="error",
                        reason=reason,
                        capacity=None,
                    )
                    write_result(paths[0], record)
        finally:
            # GPU resources must outlive all searches but should be released
            # before a later large configuration begins.
            del build.index
            if build.resources is not None:
                del build.resources
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
