"""Shared, dependency-light utilities for reproducible index benchmarks.

The benchmark runners deliberately keep their result records independent from
the individual vector-index libraries.  This module owns the stable on-disk
schema, provenance IDs, recall calculation, and lightweight resource sampling
used by those runners.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import resource
import shutil
import socket
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping, Sequence

import numpy as np


SCHEMA_VERSION = "ginflow-benchmark-v1"
RUNNER_VERSION = "1.0"
MEASUREMENT_PROTOCOL = "fixed-query-batches-v1"
LATENCY_UNIT = "milliseconds_per_query_batch"
QPS_SCOPE = "queries_per_second_over_timed_batches"
RESULT_REQUIRED_FIELDS = frozenset(
    {
        "schema_version",
        "status",
        "backend",
        "dataset_id",
        "dataset_window_count",
        "dimension",
        "metric",
        "k",
        "parameter_label",
        "parameters",
        "run_id",
        "repeat",
        "warmup_queries",
        "timed_queries",
        "qps",
        "latency_ms",
        "recall_at_100",
        "query_ids_sha256",
        "ground_truth_ids_sha256",
        "index_bytes",
        "build_seconds",
        "peak_rss_bytes",
        "peak_vram_bytes",
        "provenance",
        "error",
    }
)


def canonical_json(value: Any) -> str:
    """Render JSON in the canonical form used by stable benchmark IDs."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("utf-8"))


def sha256_file(path: Path, chunk_bytes: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    header = canonical_json({"dtype": array.dtype.str, "shape": list(array.shape)}).encode("utf-8")
    return sha256_bytes(header + array.tobytes())


def stable_id(prefix: str, value: Any, length: int = 16) -> str:
    return f"{prefix}-{sha256_json(value)[:length]}"


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def atomic_json_dump(path: Path, value: Any) -> None:
    """Write JSON atomically so an interrupted overnight run leaves no fake cache."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}-{time.time_ns()}")
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_text_dump(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}-{time.time_ns()}")
    try:
        temporary.write_text(content)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def git_commit(repo_root: Path) -> str | None:
    try:
        output = subprocess.check_output(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return output.strip() or None


def _cpu_model() -> str | None:
    cpuinfo = Path("/proc/cpuinfo")
    if not cpuinfo.exists():
        return platform.processor() or None
    for line in cpuinfo.read_text(errors="replace").splitlines():
        if line.lower().startswith("model name") and ":" in line:
            return line.split(":", 1)[1].strip()
    return platform.processor() or None


def _ram_bytes() -> int | None:
    try:
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (AttributeError, OSError, ValueError):
        return None


def _nvidia_smi_rows(fields: Sequence[str]) -> list[list[str]]:
    command = shutil.which("nvidia-smi")
    if command is None:
        return []
    try:
        output = subprocess.check_output(
            [
                command,
                f"--query-gpu={','.join(fields)}",
                "--format=csv,noheader,nounits",
            ],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=10,
        )
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return []
    return [[item.strip() for item in row.split(",")] for row in output.splitlines() if row.strip()]


def hardware_snapshot() -> dict[str, Any]:
    """Capture stable enough host hardware evidence without extra Python packages."""
    gpus = []
    for row in _nvidia_smi_rows(("name", "memory.total", "driver_version", "compute_cap")):
        if len(row) != 4:
            continue
        name, vram_mib, driver, compute_cap = row
        try:
            vram_bytes: int | None = int(float(vram_mib) * 1024 * 1024)
        except ValueError:
            vram_bytes = None
        gpus.append(
            {
                "name": name,
                "vram_bytes": vram_bytes,
                "driver_version": driver,
                "compute_capability": compute_cap,
            }
        )
    machine = {
        "cpu_model": _cpu_model(),
        "logical_cpus": os.cpu_count(),
        "ram_bytes": _ram_bytes(),
        "gpus": gpus,
    }
    machine["id"] = stable_id("hardware", machine)
    machine["hostname"] = socket.gethostname()
    machine["platform"] = platform.platform()
    return machine


def software_snapshot(repo_root: Path | None = None) -> dict[str, Any]:
    """Capture interpreter and checkout identity; backend runners add image/package data."""
    snapshot: dict[str, Any] = {
        "runner_version": RUNNER_VERSION,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }
    if repo_root is not None:
        snapshot["git_commit"] = git_commit(repo_root)
    return snapshot


def percentile(values: Sequence[float], quantile: float) -> float | None:
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=np.float64), quantile))


def recall_at_k(labels: np.ndarray, expected: np.ndarray, k: int = 100) -> float:
    """Mean set recall@k, ignoring missing ``-1`` labels from an ANN backend."""
    got = np.asarray(labels, dtype=np.int64)
    truth = np.asarray(expected, dtype=np.int64)
    if got.ndim != 2 or truth.ndim != 2:
        raise ValueError("labels and expected neighbours must both be 2-D")
    if got.shape[0] != truth.shape[0]:
        raise ValueError("labels and expected neighbours have different query counts")
    width = min(int(k), got.shape[1], truth.shape[1])
    if width < 1:
        raise ValueError("k must leave at least one neighbour")
    recalls: list[float] = []
    for predicted, target in zip(got[:, :width], truth[:, :width]):
        target_set = {int(value) for value in target if int(value) >= 0}
        if not target_set:
            continue
        predicted_set = {int(value) for value in predicted if int(value) >= 0}
        recalls.append(len(predicted_set.intersection(target_set)) / len(target_set))
    if not recalls:
        return 0.0
    return float(np.mean(recalls))


def _current_rss_bytes() -> int | None:
    status = Path("/proc/self/status")
    if status.exists():
        match = re.search(r"^VmRSS:\s+(\d+)\s+kB", status.read_text(errors="replace"), re.MULTILINE)
        if match:
            return int(match.group(1)) * 1024
    try:
        # Linux ru_maxrss is KiB, macOS reports bytes. This fallback is only
        # used where /proc is absent, so distinguish Darwin explicitly.
        maximum = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return maximum if sys.platform == "darwin" else maximum * 1024
    except (AttributeError, ValueError):
        return None


def _current_vram_bytes(device: int | None = None) -> int | None:
    rows = _nvidia_smi_rows(("memory.used",))
    if not rows:
        return None
    selected = rows if device is None else rows[device : device + 1]
    values = []
    for row in selected:
        try:
            values.append(int(float(row[0])) * 1024 * 1024)
        except (IndexError, ValueError):
            continue
    return max(values) if values else None


class PeakResourceSampler:
    """Sample process RSS and GPU-wide VRAM during an in-process benchmark.

    VRAM is the NVIDIA-device total reported by ``nvidia-smi``. It is a useful
    upper bound on a single-user benchmark host, rather than a per-process
    attribution; that distinction is preserved in the result provenance.
    """

    def __init__(self, interval_seconds: float = 0.2, gpu_device: int | None = None) -> None:
        if interval_seconds <= 0:
            raise ValueError("interval_seconds must be > 0")
        self.interval_seconds = float(interval_seconds)
        self.gpu_device = gpu_device
        self.peak_rss_bytes: int | None = None
        self.peak_vram_bytes: int | None = None
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _sample(self) -> None:
        rss = _current_rss_bytes()
        vram = _current_vram_bytes(self.gpu_device)
        if rss is not None:
            self.peak_rss_bytes = max(self.peak_rss_bytes or 0, rss)
        if vram is not None:
            self.peak_vram_bytes = max(self.peak_vram_bytes or 0, vram)

    def _run(self) -> None:
        while not self._stop.wait(self.interval_seconds):
            self._sample()

    def __enter__(self) -> "PeakResourceSampler":
        self._sample()
        self._thread = threading.Thread(target=self._run, name="benchmark-resource-sampler", daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_: Any) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval_seconds * 3)
        self._sample()


@dataclass(frozen=True)
class SearchMeasurement:
    repeat: int
    qps: float
    latency_ms: dict[str, float]
    recall_at_100: float
    peak_rss_bytes: int | None
    peak_vram_bytes: int | None
    measurement: dict[str, Any]


def _labels_from_search_result(value: Any) -> np.ndarray:
    if isinstance(value, tuple) and len(value) == 2:
        # GINflow wrappers return (distances, labels).
        value = value[1]
    labels = np.asarray(value, dtype=np.int64)
    if labels.ndim != 2:
        raise ValueError("the benchmark search callable must return a 2-D label array or (scores, labels)")
    return labels


def _validate_search_labels(labels: np.ndarray, *, query_count: int, k: int) -> None:
    """Reject incomplete search output before it can be called recall@k."""
    if labels.shape[0] != query_count:
        raise ValueError("search result has a different query count")
    if labels.shape[1] < k:
        raise ValueError(f"search result has {labels.shape[1]} neighbours, but k={k} was requested")


def measure_search_repeats(
    search: Callable[[np.ndarray, int], Any],
    queries: np.ndarray,
    ground_truth_ids: np.ndarray,
    *,
    k: int = 100,
    warmup_queries: int = 32,
    repeats: int = 3,
    query_batch_size: int = 32,
    gpu_device: int | None = None,
    synchronize: Callable[[], None] | None = None,
) -> Iterator[SearchMeasurement]:
    """Measure identical fixed query batches after a whole-batch warm-up.

    Each latency quantile is calculated from actual search-call durations for
    batches of ``query_batch_size`` queries. QPS is calculated from precisely
    those same calls, never from an unrelated whole-matrix average. The query
    and warm-up counts must therefore be positive multiples of the batch size.
    GPU runners may supply ``synchronize`` so asynchronous kernels are included
    in the per-batch elapsed time.
    """
    query_matrix = np.ascontiguousarray(queries, dtype=np.float32)
    if query_matrix.ndim != 2 or query_matrix.shape[0] < 1:
        raise ValueError("queries must be a non-empty 2-D float32 matrix")
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    if k < 1:
        raise ValueError("k must be >= 1")
    batch_size = int(query_batch_size)
    if batch_size < 1:
        raise ValueError("query_batch_size must be >= 1")
    if batch_size > query_matrix.shape[0]:
        raise ValueError("query_batch_size must not exceed the timed query count")
    if query_matrix.shape[0] % batch_size:
        raise ValueError("timed query count must be divisible by query_batch_size")
    warmup = int(warmup_queries)
    if warmup < 1:
        raise ValueError("warmup_queries must be >= 1")
    if warmup > query_matrix.shape[0]:
        raise ValueError("warmup_queries must not exceed the timed query count")
    if warmup % batch_size:
        raise ValueError("warmup_queries must be divisible by query_batch_size")

    def run_search(batch: np.ndarray) -> np.ndarray:
        if synchronize is not None:
            synchronize()
        labels = _labels_from_search_result(search(batch, k))
        if synchronize is not None:
            synchronize()
        _validate_search_labels(labels, query_count=batch.shape[0], k=k)
        return labels

    for warmup_start in range(0, warmup, batch_size):
        run_search(query_matrix[warmup_start : warmup_start + batch_size])

    timed_batch_count = query_matrix.shape[0] // batch_size
    measurement_template = {
        "protocol": MEASUREMENT_PROTOCOL,
        "query_batch_size": batch_size,
        "timed_batch_count": timed_batch_count,
        "timed_queries": int(query_matrix.shape[0]),
        "latency_unit": LATENCY_UNIT,
        "qps_scope": QPS_SCOPE,
    }

    for repeat in range(repeats):
        batch_latencies_seconds: list[float] = []
        label_batches: list[np.ndarray] = []
        with PeakResourceSampler(gpu_device=gpu_device) as sampler:
            for batch_start in range(0, query_matrix.shape[0], batch_size):
                batch = query_matrix[batch_start : batch_start + batch_size]
                started = time.perf_counter()
                label_batches.append(run_search(batch))
                batch_latencies_seconds.append(time.perf_counter() - started)
        elapsed_seconds = float(sum(batch_latencies_seconds))
        if elapsed_seconds <= 0:
            raise RuntimeError("search timing did not advance")
        labels = np.concatenate(label_batches, axis=0)
        _validate_search_labels(labels, query_count=query_matrix.shape[0], k=k)
        batch_latencies_ms = [value * 1000.0 for value in batch_latencies_seconds]
        yield SearchMeasurement(
            repeat=repeat,
            qps=float(query_matrix.shape[0] / elapsed_seconds) if elapsed_seconds > 0 else float("inf"),
            latency_ms={
                "mean": float(np.mean(batch_latencies_ms)),
                "p50": percentile(batch_latencies_ms, 50.0) or 0.0,
                "p95": percentile(batch_latencies_ms, 95.0) or 0.0,
            },
            recall_at_100=recall_at_k(labels, ground_truth_ids, k=k),
            peak_rss_bytes=sampler.peak_rss_bytes,
            peak_vram_bytes=sampler.peak_vram_bytes,
            measurement={**measurement_template, "timed_seconds": elapsed_seconds},
        )


def make_result_record(
    *,
    backend: str,
    dataset_id: str,
    dataset_window_count: int,
    dimension: int,
    parameter_label: str,
    parameters: Mapping[str, Any],
    run_id: str,
    repeat: int,
    warmup_queries: int,
    timed_queries: int,
    query_ids_sha256: str,
    ground_truth_ids_sha256: str,
    provenance: Mapping[str, Any],
    status: str = "ok",
    metric: str = "cosine",
    k: int = 100,
    qps: float | None = None,
    latency_ms: Mapping[str, float] | None = None,
    recall: float | None = None,
    index_bytes: int | None = None,
    build_seconds: float | None = None,
    peak_rss_bytes: int | None = None,
    peak_vram_bytes: int | None = None,
    error: str | None = None,
    measurement: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Create a flat, report-friendly result object for one timed repeat."""
    if status not in {"ok", "skipped", "error"}:
        raise ValueError("status must be one of ok, skipped, error")
    result_measurement = dict(measurement or {})
    result_provenance = dict(provenance)
    if status == "ok" and result_measurement:
        result_provenance.setdefault("measurement_protocol", result_measurement.get("protocol"))
        result_provenance.setdefault("query_batch_size", result_measurement.get("query_batch_size"))
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": status,
        "backend": backend,
        "dataset_id": dataset_id,
        "dataset_window_count": int(dataset_window_count),
        "dimension": int(dimension),
        "metric": metric,
        "k": int(k),
        "parameter_label": parameter_label,
        "parameters": dict(parameters),
        "run_id": run_id,
        "repeat": int(repeat),
        "warmup_queries": int(warmup_queries),
        "timed_queries": int(timed_queries),
        "qps": qps,
        "latency_ms": None if latency_ms is None else dict(latency_ms),
        "recall_at_100": recall,
        "query_ids_sha256": query_ids_sha256,
        "ground_truth_ids_sha256": ground_truth_ids_sha256,
        "index_bytes": index_bytes,
        "build_seconds": build_seconds,
        "peak_rss_bytes": peak_rss_bytes,
        "peak_vram_bytes": peak_vram_bytes,
        "provenance": result_provenance,
        "error": error,
        "measurement": result_measurement,
        "created_at": utc_now(),
    }
    errors = validate_result_record(record)
    if errors:
        raise ValueError("invalid benchmark result: " + "; ".join(errors))
    return record


def validate_result_record(record: Mapping[str, Any]) -> list[str]:
    """Validate the stable schema without imposing a jsonschema dependency."""
    errors: list[str] = []
    missing = sorted(RESULT_REQUIRED_FIELDS.difference(record))
    if missing:
        errors.append("missing fields: " + ", ".join(missing))
        return errors
    if record.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must be {SCHEMA_VERSION}")
    status = record.get("status")
    if status not in {"ok", "skipped", "error"}:
        errors.append("status must be ok, skipped, or error")
    if int(record.get("k") or 0) != 100:
        errors.append("this benchmark schema is defined for k=100")
    if record.get("metric") != "cosine":
        errors.append("metric must be cosine")
    if int(record.get("dataset_window_count") or 0) < 1:
        errors.append("dataset_window_count must be >= 1")
    if int(record.get("dimension") or 0) < 1:
        errors.append("dimension must be >= 1")
    if int(record.get("timed_queries") or 0) < 0:
        errors.append("timed_queries must be >= 0")
    if int(record.get("warmup_queries") or 0) < 0:
        errors.append("warmup_queries must be >= 0")
    provenance = record.get("provenance")
    if not isinstance(provenance, Mapping):
        errors.append("provenance must be an object")
    elif status == "ok":
        for key in (
            "git_commit",
            "runner_version",
            "hardware_id",
            "embedding_cache_id",
            "query_selection_id",
            "ground_truth_cache_id",
            "measurement_protocol",
            "query_batch_size",
        ):
            if not provenance.get(key):
                errors.append(f"provenance.{key} is required for successful results")
    if status == "ok":
        for key in ("qps", "recall_at_100", "index_bytes", "build_seconds"):
            if record.get(key) is None:
                errors.append(f"{key} is required for successful results")
        recall = record.get("recall_at_100")
        if recall is not None and not (0.0 <= float(recall) <= 1.0):
            errors.append("recall_at_100 must be between 0 and 1")
        if int(record.get("warmup_queries") or 0) < 1:
            errors.append("warmup_queries must be >= 1 for successful results")
        if int(record.get("timed_queries") or 0) < 1:
            errors.append("timed_queries must be >= 1 for successful results")
        latency = record.get("latency_ms")
        if not isinstance(latency, Mapping):
            errors.append("latency_ms must be an object for successful results")
        else:
            values: dict[str, float] = {}
            for key in ("mean", "p50", "p95"):
                value = latency.get(key)
                if not isinstance(value, (int, float)) or isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                    errors.append(f"latency_ms.{key} must be a positive finite number")
                else:
                    values[key] = float(value)
            if values.get("p50", 0.0) > values.get("p95", float("inf")):
                errors.append("latency_ms.p50 must not exceed latency_ms.p95")
        benchmark_measurement = record.get("measurement")
        if not isinstance(benchmark_measurement, Mapping):
            errors.append("measurement must be an object for successful results")
        else:
            if benchmark_measurement.get("protocol") != MEASUREMENT_PROTOCOL:
                errors.append(f"measurement.protocol must be {MEASUREMENT_PROTOCOL}")
            if benchmark_measurement.get("latency_unit") != LATENCY_UNIT:
                errors.append(f"measurement.latency_unit must be {LATENCY_UNIT}")
            if benchmark_measurement.get("qps_scope") != QPS_SCOPE:
                errors.append(f"measurement.qps_scope must be {QPS_SCOPE}")
            batch_size = benchmark_measurement.get("query_batch_size")
            batch_count = benchmark_measurement.get("timed_batch_count")
            timed = int(record.get("timed_queries") or 0)
            if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size < 1:
                errors.append("measurement.query_batch_size must be a positive integer")
            elif timed and (batch_size > timed or timed % batch_size):
                errors.append("timed_queries must be a positive multiple of measurement.query_batch_size")
            if not isinstance(batch_count, int) or isinstance(batch_count, bool) or batch_count < 1:
                errors.append("measurement.timed_batch_count must be a positive integer")
            elif isinstance(batch_size, int) and batch_size > 0 and timed and batch_count != timed // batch_size:
                errors.append("measurement.timed_batch_count must equal timed_queries / query_batch_size")
            if benchmark_measurement.get("timed_queries") != timed:
                errors.append("measurement.timed_queries must equal timed_queries")
            if provenance and provenance.get("measurement_protocol") != benchmark_measurement.get("protocol"):
                errors.append("provenance.measurement_protocol must equal measurement.protocol")
            if provenance and provenance.get("query_batch_size") != batch_size:
                errors.append("provenance.query_batch_size must equal measurement.query_batch_size")
    if status == "error" and not record.get("error"):
        errors.append("error is required when status=error")
    return errors


def write_result(path: Path, record: Mapping[str, Any]) -> None:
    errors = validate_result_record(record)
    if errors:
        raise ValueError("invalid benchmark result: " + "; ".join(errors))
    atomic_json_dump(path, dict(record))
