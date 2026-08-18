#!/usr/bin/env python3
"""Create reproducible vector, query, and exact-neighbour benchmark caches.

The large graph/embed/window step is delegated to ``cache_windows.nf``.  The
remaining operations deliberately stream one window shard or one bounded vector
block at a time, so the 30k corpus never needs a second whole-database array in
RAM (or GPU VRAM).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Sequence

import numpy as np

from benchmark_utils import (
    atomic_json_dump,
    atomic_text_dump,
    hardware_snapshot,
    sha256_array,
    sha256_file,
    sha256_json,
    software_snapshot,
    stable_id,
    utc_now,
)


DEFAULT_QUERY_COUNT = 512
DEFAULT_QUERY_SEED = "ginflow-index-benchmark-v1"
DEFAULT_DATABASE_CHUNK = 32_768
DEFAULT_QUERY_CHUNK = 16
WINDOW_SUFFIX = ".windows.npz"
WINDOW_MANIFEST_SUFFIX = ".windows.manifest.json"
WINDOW_BINARY_PATH = "bin/slice_windows.py"
PREPARATION_IMPLEMENTATION_PATHS = (
    "benchmarks/cache_windows.nf",
    "workflows/prepare_windows.nf",
    "modules/build_rna_graphs/main.nf",
    "modules/build_rna_graphs/environment.yml",
    "modules/embed_rna_graphs/main.nf",
    "modules/embed_rna_graphs/environment.yml",
    "modules/embed_rna_graphs/environment.gpu.yml",
    "modules/generate_windows/main.nf",
    "modules/generate_windows/environment.yml",
    WINDOW_BINARY_PATH,
    "nextflow.config",
    "conf/base.config",
)


def cache_root_from_args(value: str | None) -> Path:
    if value:
        return Path(value).expanduser().resolve()
    configured = os.environ.get("GINFLOW_BENCHMARK_CACHE")
    if configured:
        return Path(configured).expanduser().resolve()
    return (Path.cwd() / ".benchmark-cache").resolve()


def dataset_root(cache_dir: Path, dataset_id: str) -> Path:
    if not dataset_id or any(char in dataset_id for char in "/\\"):
        raise ValueError("dataset id must be a non-empty simple name")
    return cache_dir / dataset_id


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def _read_json_list(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, list) or not all(isinstance(item, dict) for item in payload):
        raise ValueError(f"{path} is not a list of JSON objects")
    return payload


def _relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _cache_provenance_path(root: Path) -> Path:
    return root / "provenance.json"


def update_cache_provenance(root: Path, updates: dict[str, Any]) -> dict[str, Any]:
    """Merge a phase's evidence into a cache-wide provenance document."""
    destination = _cache_provenance_path(root)
    existing = _read_json(destination) if destination.exists() else {}
    existing.update(updates)
    existing.setdefault("schema_version", "ginflow-benchmark-cache-v1")
    existing.setdefault("created_at", utc_now())
    existing["updated_at"] = utc_now()
    atomic_json_dump(destination, existing)
    return existing


def _profile_names(profile: str) -> list[str]:
    names = [item.strip() for item in profile.split(",") if item.strip()]
    if not names:
        raise ValueError("profile must name at least one Nextflow profile")
    return names


def window_binary_integrity(repo_root: Path) -> dict[str, str]:
    """Validate the canonical executable in the repository ``bin/`` directory."""
    entrypoint = repo_root / WINDOW_BINARY_PATH
    if not entrypoint.is_file() or not os.access(entrypoint, os.X_OK):
        raise ValueError(f"window binary is missing or not executable: {entrypoint}")
    if entrypoint.is_symlink():
        raise ValueError(f"window binary must be a regular file in bin/: {entrypoint}")
    return {
        "entrypoint": WINDOW_BINARY_PATH,
        "sha256": sha256_file(entrypoint),
    }


def preparation_implementation_files(repo_root: Path, profile: str) -> list[Path]:
    """Return every live source/config file that can change cached windows.

    The list follows the benchmark entry workflow through its three processes,
    including their environment definitions, the called window script, common
    configuration, and any built-in profile config named by ``profile``. Paths
    are repository-relative in the fingerprint, so moving a checkout does not
    change the implementation identity.
    """
    relative_paths = set(PREPARATION_IMPLEMENTATION_PATHS)
    for profile_name in _profile_names(profile):
        candidate = f"conf/{profile_name}.config"
        if (repo_root / candidate).is_file():
            relative_paths.add(candidate)
    files = [repo_root / relative for relative in sorted(relative_paths)]
    missing = [str(path.relative_to(repo_root)) for path in files if not path.is_file()]
    if missing:
        raise ValueError("window-preparation implementation file(s) missing: " + ", ".join(missing))
    window_binary_integrity(repo_root)
    return files


def preparation_implementation_fingerprint(repo_root: Path, profile: str) -> dict[str, Any]:
    """Hash live preparation code/configuration without relying on Git state."""
    files = preparation_implementation_files(repo_root, profile)
    entries = [
        {"path": path.relative_to(repo_root).as_posix(), "sha256": sha256_file(path)}
        for path in files
    ]
    identity = {
        "schema_version": "ginflow-window-preparation-implementation-v2",
        "files": entries,
        "window_binary": window_binary_integrity(repo_root),
    }
    return {**identity, "id": stable_id("window-preparation", identity)}


def _window_cache_request_identity(request: dict[str, Any]) -> dict[str, Any] | None:
    """Read the semantic request identity, refusing old path/timestamp-only caches."""
    identity = request.get("identity")
    return dict(identity) if isinstance(identity, dict) else None


@contextmanager
def _cache_lock(root: Path) -> Iterator[None]:
    """Serialize a dataset cache so two overnight invocations cannot interleave."""
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / ".window-cache.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        except ImportError:  # pragma: no cover - the benchmark host is POSIX
            pass
        try:
            yield
        finally:
            try:
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except ImportError:  # pragma: no cover - the benchmark host is POSIX
                pass


def cache_windows_request(
    repo_root: Path,
    input_path: Path,
    dataset_id: str,
    *,
    profile: str,
    shard_size: int,
    window_size: int,
    window_stride: int,
) -> dict[str, Any]:
    input_sha256 = sha256_file(input_path)
    implementation = preparation_implementation_fingerprint(repo_root, profile)
    identity = {
        "schema_version": "ginflow-window-cache-request-identity-v2",
        "dataset_id": dataset_id,
        "input_sha256": input_sha256,
        "profile": profile,
        "shard_size": int(shard_size),
        "window_size": int(window_size),
        "window_stride": int(window_stride),
        "preparation_implementation_id": implementation["id"],
    }
    return {
        "schema_version": "ginflow-window-cache-request-v2",
        "identity": identity,
        "request_id": stable_id("window-cache-request", identity),
        "dataset_id": dataset_id,
        "input_path": str(input_path.resolve()),
        "input_sha256": input_sha256,
        "profile": profile,
        "shard_size": int(shard_size),
        "window_size": int(window_size),
        "window_stride": int(window_stride),
        "preparation_implementation": implementation,
        "created_at": utc_now(),
    }


def _artifacts_ready(artifacts: Path) -> bool:
    try:
        graph_paths = sorted((artifacts / "graphs").rglob("*.safetensors"))
        embedding_paths = sorted((artifacts / "embeddings").rglob("*.npz"))
        window_pairs = find_window_shards(artifacts / "windows")
    except (OSError, ValueError):
        return False
    if not graph_paths or not embedding_paths or not window_pairs:
        return False
    if not all(path.with_suffix(".json").is_file() for path in graph_paths):
        return False
    if not all(path.with_suffix(".manifest.json").is_file() for path in embedding_paths):
        return False
    # PREPARE_WINDOWS emits one graph, embedding, and window shard per input
    # shard. A count mismatch means publication was interrupted or mixed.
    return len(graph_paths) == len(embedding_paths) == len(window_pairs)


def run_window_cache(
    *,
    repo_root: Path,
    cache_dir: Path,
    dataset_id: str,
    input_path: Path,
    profile: str,
    shard_size: int,
    window_size: int,
    window_stride: int,
    nextflow: str,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Materialize graph/embed/window artifacts once using Nextflow ``-resume``.

    Existing artifacts are reused only when their immutable semantic request
    matches exactly. A failed attempt leaves its request as ``in_progress`` or
    ``failed`` so the same invocation can safely resume; only a complete request
    with a complete artifact set is reused.
    """
    if not input_path.is_file():
        raise ValueError(f"dataset TSV does not exist: {input_path}")
    root = dataset_root(cache_dir, dataset_id)
    artifacts = root / "artifacts"
    request = cache_windows_request(
        repo_root,
        input_path,
        dataset_id,
        profile=profile,
        shard_size=shard_size,
        window_size=window_size,
        window_stride=window_stride,
    )
    request_path = root / "window-cache-request.json"
    work_dir = root / "nextflow-work"
    command = [
        nextflow,
        "run",
        str(repo_root / "benchmarks" / "cache_windows.nf"),
        "-c",
        str(repo_root / "nextflow.config"),
        "-profile",
        profile,
        "-resume",
        "-work-dir",
        str(work_dir),
        "-output-dir",
        str(artifacts),
        "--input",
        str(input_path.resolve()),
        "--benchmark_prefix",
        dataset_id,
        "--shard_size",
        str(shard_size),
        "--window_size",
        str(window_size),
        "--window_stride",
        str(window_stride),
    ]
    execution = {
        "command": command,
        "started_at": utc_now(),
        "hardware": hardware_snapshot(),
        "software": software_snapshot(repo_root),
    }
    if dry_run:
        return {"status": "planned", "root": str(root), "request": request, "execution": execution}
    with _cache_lock(root):
        if request_path.exists():
            prior = _read_json(request_path)
            prior_identity = _window_cache_request_identity(prior)
            if prior_identity != request["identity"]:
                raise ValueError(
                    "cache request differs from existing artifacts; choose a new --cache-dir "
                    f"or dataset id. Existing={prior_identity!r}, requested={request['identity']!r}"
                )
            if prior.get("state") == "complete" and _artifacts_ready(artifacts):
                return {"status": "reused", "root": str(root), "request": prior}
            request["attempt"] = int(prior.get("attempt") or 0) + 1
        elif artifacts.exists() and any(artifacts.iterdir()):
            raise ValueError(
                f"{artifacts} already contains files without a cache request; refusing to mix provenance"
            )
        else:
            request["attempt"] = 1

        request["state"] = "in_progress"
        request["execution"] = execution
        atomic_json_dump(request_path, request)
        try:
            subprocess.run(command, cwd=repo_root, check=True)
            if not _artifacts_ready(artifacts):
                raise RuntimeError("Nextflow completed but the required graph/embed/window artifacts were not published")
        except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
            request["state"] = "failed"
            request["failed_at"] = utc_now()
            request["error"] = str(exc)
            atomic_json_dump(request_path, request)
            raise

        request["state"] = "complete"
        request["completed_at"] = utc_now()
        atomic_json_dump(request_path, request)
        update_cache_provenance(
            root,
            {
                "dataset_id": dataset_id,
                "source": {"path": str(input_path.resolve()), "sha256": request["input_sha256"]},
                "window_cache_request_id": request["request_id"],
                "window_cache": {"request": _relative_or_absolute(request_path, root), "execution": execution},
            },
        )
        return {"status": "created", "root": str(root), "request": request, "execution": execution}


def find_window_shards(windows_dir: Path) -> list[tuple[Path, Path]]:
    if not windows_dir.is_dir():
        raise ValueError(f"windows directory does not exist: {windows_dir}")
    pairs: list[tuple[Path, Path]] = []
    for window_path in sorted(windows_dir.rglob(f"*{WINDOW_SUFFIX}")):
        stem = window_path.name[: -len(WINDOW_SUFFIX)]
        manifest_path = window_path.with_name(stem + WINDOW_MANIFEST_SUFFIX)
        if not manifest_path.is_file():
            raise ValueError(f"missing manifest for {window_path}: expected {manifest_path.name}")
        pairs.append((window_path, manifest_path))
    if not pairs:
        raise ValueError(f"no {WINDOW_SUFFIX} shards found in {windows_dir}")
    return pairs


def inspect_window_shards(windows_dir: Path) -> tuple[list[dict[str, Any]], int, int, dict[str, Any]]:
    """First pass: validate manifests and determine the exact memmap shape."""
    shard_entries: list[dict[str, Any]] = []
    total_windows = 0
    window_dim: int | None = None
    reference: dict[str, Any] | None = None
    seen_identifiers: set[str] = set()
    for window_path, manifest_path in find_window_shards(windows_dir):
        manifest = _read_json(manifest_path)
        if not manifest.get("l2_normalized"):
            raise ValueError(f"{manifest_path} does not declare l2_normalized=true")
        dim = int(manifest.get("window_dim") or 0)
        if dim < 1:
            raise ValueError(f"{manifest_path} has invalid window_dim")
        if window_dim is None:
            window_dim = dim
            reference = manifest
        elif dim != window_dim:
            raise ValueError(f"{manifest_path} has window_dim={dim}, expected {window_dim}")
        records = manifest.get("records")
        if not isinstance(records, list):
            raise ValueError(f"{manifest_path} has no records list")
        count = 0
        for item in records:
            if not isinstance(item, dict):
                raise ValueError(f"{manifest_path} contains a non-object record")
            identifier = str(item.get("identifier") or "")
            n_windows = int(item.get("n_windows") or 0)
            if not identifier or n_windows < 1:
                raise ValueError(f"{manifest_path} has an invalid record entry")
            if identifier in seen_identifiers:
                raise ValueError(f"duplicate record identifier {identifier!r} across window shards")
            seen_identifiers.add(identifier)
            count += n_windows
        shard_entries.append(
            {
                "windows": str(window_path),
                "manifest": str(manifest_path),
                "manifest_sha256": sha256_file(manifest_path),
                "n_windows": count,
                "n_records": len(records),
            }
        )
        total_windows += count
    if window_dim is None or reference is None or total_windows < 1:
        raise ValueError("window manifests contain no indexable vectors")
    return shard_entries, total_windows, window_dim, reference


def _normalise_block(block: np.ndarray) -> np.ndarray:
    output = np.ascontiguousarray(block, dtype=np.float32)
    if output.ndim != 2:
        raise ValueError(f"window block must be 2-D, got {output.shape}")
    if not np.isfinite(output).all():
        raise ValueError("window vectors contain non-finite values")
    norms = np.linalg.norm(output, axis=1, keepdims=True)
    if np.any(norms <= np.float32(1e-12)):
        raise ValueError("window vectors contain zero-norm rows")
    output /= norms
    return output


def flatten_window_shards(
    windows_dir: Path,
    flat_dir: Path,
    *,
    dataset_id: str,
    row_chunk: int = 4096,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Flatten sharded compressed windows to a normalized ``.npy`` memmap.

    ``row_chunk`` is the largest independent row block held in addition to one
    compressed shard.  The output itself is an mmap, not a concatenated RAM
    matrix, which keeps the 30k corpus streamable on ordinary workstations.
    """
    if row_chunk < 1:
        raise ValueError("row_chunk must be >= 1")
    shard_entries, total_windows, window_dim, reference = inspect_window_shards(windows_dir)
    flat_dir.mkdir(parents=True, exist_ok=True)
    vectors_path = flat_dir / "vectors.npy"
    records_path = flat_dir / "records.json"
    windows_path = flat_dir / "windows.tsv"
    manifest_path = flat_dir / "flatten-manifest.json"
    outputs = (vectors_path, records_path, windows_path, manifest_path)
    if any(path.exists() for path in outputs) and not overwrite:
        raise ValueError(f"flattened output already exists in {flat_dir}; use --overwrite only after review")

    temporary_vectors = flat_dir / ".vectors.partial.npy"
    temporary_records = flat_dir / ".records.partial.json"
    temporary_windows = flat_dir / ".windows.partial.tsv"
    if any(path.exists() for path in (temporary_vectors, temporary_records, temporary_windows)):
        raise ValueError(f"partial flatten output exists in {flat_dir}; inspect it before retrying")

    vector_digest = hashlib.sha256()
    records: list[dict[str, Any]] = []
    cursor = 0
    memmap = np.lib.format.open_memmap(
        temporary_vectors,
        mode="w+",
        dtype=np.float32,
        shape=(total_windows, window_dim),
    )
    try:
        with temporary_windows.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["vector_id", "record_ordinal", "transcript_id", "window_offset"])
            for entry in shard_entries:
                window_file = Path(entry["windows"])
                manifest_file = Path(entry["manifest"])
                manifest = _read_json(manifest_file)
                stride = int(manifest.get("stride") or 0)
                if stride < 1:
                    raise ValueError(f"{manifest_file} has invalid stride")
                with np.load(window_file, allow_pickle=False) as archive:
                    for source_record in manifest["records"]:
                        identifier = str(source_record["identifier"])
                        expected_rows = int(source_record["n_windows"])
                        if identifier not in archive.files:
                            raise ValueError(f"{identifier} is listed in {manifest_file} but absent from {window_file}")
                        array = np.asarray(archive[identifier], dtype=np.float32)
                        if array.shape != (expected_rows, window_dim):
                            raise ValueError(
                                f"{identifier} shape {array.shape} disagrees with manifest "
                                f"({expected_rows}, {window_dim})"
                            )
                        start = cursor
                        for chunk_start in range(0, expected_rows, row_chunk):
                            chunk_end = min(expected_rows, chunk_start + row_chunk)
                            normalised = _normalise_block(array[chunk_start:chunk_end])
                            memmap[cursor : cursor + normalised.shape[0]] = normalised
                            vector_digest.update(normalised.tobytes())
                            cursor += normalised.shape[0]
                        record_ordinal = len(records)
                        records.append(
                            {
                                "record_ordinal": record_ordinal,
                                "identifier": identifier,
                                "vector_start": start,
                                "n_windows": expected_rows,
                                "window_stride": stride,
                            }
                        )
                        writer.writerows(
                            (start + offset, record_ordinal, identifier, offset * stride)
                            for offset in range(expected_rows)
                        )
        if cursor != total_windows:
            raise RuntimeError(f"flatten cursor={cursor}, expected {total_windows}")
        memmap.flush()
    finally:
        # Release the mmap before atomically publishing it on all supported OSes.
        del memmap

    atomic_json_dump(temporary_records, records)
    records_sha256 = sha256_json(records)
    windows_sha256 = sha256_file(temporary_windows)
    embedding_identity = {
        "schema_version": "ginflow-embedding-cache-identity-v1",
        "dataset_id": dataset_id,
        "vectors": {
            "shape": [total_windows, window_dim],
            "dtype": "float32",
            "l2_normalized": True,
            "payload_sha256": vector_digest.hexdigest(),
        },
        "records_sha256": records_sha256,
        "windows_sha256": windows_sha256,
        "window_parameters": {
            "window_size": int(reference.get("window_size") or 0),
            "stride": int(reference.get("stride") or 0),
            "embedding_dim": int(reference.get("embedding_dim") or 0),
            "window_dim": window_dim,
        },
        "embedding_provenance": {
            "ginfinity_version": reference.get("ginfinity_version"),
            "model_version": reference.get("model_version"),
            "checkpoint_sha256": reference.get("checkpoint_sha256"),
        },
    }
    flatten_manifest = {
        "schema_version": "ginflow-flatten-cache-v1",
        "dataset_id": dataset_id,
        "created_at": utc_now(),
        "vectors": {
            "path": vectors_path.name,
            "shape": [total_windows, window_dim],
            "dtype": "float32",
            "l2_normalized": True,
            "payload_sha256": vector_digest.hexdigest(),
        },
        "records": {"path": records_path.name, "count": len(records)},
        "windows": {"path": windows_path.name, "count": total_windows, "sha256": windows_sha256},
        "window_parameters": {
            "window_size": int(reference.get("window_size") or 0),
            "stride": int(reference.get("stride") or 0),
            "embedding_dim": int(reference.get("embedding_dim") or 0),
            "window_dim": window_dim,
        },
        "embedding_provenance": {
            "ginfinity_version": reference.get("ginfinity_version"),
            "model_version": reference.get("model_version"),
            "checkpoint_sha256": reference.get("checkpoint_sha256"),
        },
        "source_shards": shard_entries,
        "identity": embedding_identity,
    }
    flatten_manifest["embedding_cache_id"] = stable_id("embedding-cache", embedding_identity)
    temporary_manifest = flat_dir / ".flatten-manifest.partial.json"
    atomic_json_dump(temporary_manifest, flatten_manifest)
    os.replace(temporary_vectors, vectors_path)
    os.replace(temporary_records, records_path)
    os.replace(temporary_windows, windows_path)
    os.replace(temporary_manifest, manifest_path)
    return flatten_manifest


def _hash_to_int(*parts: str) -> int:
    digest = hashlib.sha256("\x00".join(parts).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False)


def select_record_balanced_queries(
    records: Sequence[dict[str, Any]],
    vectors: np.ndarray,
    *,
    query_count: int = DEFAULT_QUERY_COUNT,
    seed: str = DEFAULT_QUERY_SEED,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Select deterministic queries with equal record weight, not length weight."""
    if query_count < 1:
        raise ValueError("query_count must be >= 1")
    if not records:
        raise ValueError("cannot select queries from zero records")
    sorted_records = sorted(
        records,
        key=lambda item: (_hash_to_int(seed, str(item["identifier"]), "record"), str(item["identifier"])),
    )
    chosen: list[dict[str, Any]] = []
    cycle = 0
    while len(chosen) < query_count:
        for record in sorted_records:
            if len(chosen) >= query_count:
                break
            n_windows = int(record["n_windows"])
            if n_windows < 1:
                continue
            identifier = str(record["identifier"])
            offset = _hash_to_int(seed, identifier, "window", str(cycle)) % n_windows
            flat_id = int(record["vector_start"]) + offset
            chosen.append(
                {
                    "query_ordinal": len(chosen),
                    "vector_id": flat_id,
                    "record_ordinal": int(record["record_ordinal"]),
                    "transcript_id": identifier,
                    "window_offset": offset * int(record["window_stride"]),
                    "selection_cycle": cycle,
                }
            )
        cycle += 1
    ids = np.asarray([item["vector_id"] for item in chosen], dtype=np.int64)
    query_vectors = np.ascontiguousarray(vectors[ids], dtype=np.float32)
    norms = np.linalg.norm(query_vectors, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(np.abs(norms - 1.0) > 2e-5):
        raise ValueError("flattened cache does not contain unit-normalized query vectors")
    return query_vectors, chosen


def make_query_cache(
    flat_dir: Path,
    query_dir: Path,
    *,
    dataset_id: str,
    query_count: int = DEFAULT_QUERY_COUNT,
    seed: str = DEFAULT_QUERY_SEED,
    overwrite: bool = False,
) -> dict[str, Any]:
    flatten_manifest = _read_json(flat_dir / "flatten-manifest.json")
    records = _read_json_list(flat_dir / "records.json")
    vectors = np.load(flat_dir / "vectors.npy", mmap_mode="r", allow_pickle=False)
    query_dir.mkdir(parents=True, exist_ok=True)
    vectors_path = query_dir / "queries.npy"
    table_path = query_dir / "queries.tsv"
    manifest_path = query_dir / "query-selection.json"
    if any(path.exists() for path in (vectors_path, table_path, manifest_path)) and not overwrite:
        raise ValueError(f"query cache already exists in {query_dir}; use --overwrite only after review")
    partials = (
        query_dir / ".queries.partial.npy",
        query_dir / ".queries.partial.tsv",
        query_dir / ".query-selection.partial.json",
    )
    if any(path.exists() for path in partials):
        raise ValueError(f"partial query cache exists in {query_dir}; inspect it before retrying")
    queries, selected = select_record_balanced_queries(
        records,
        vectors,
        query_count=query_count,
        seed=seed,
    )
    temporary_vectors = query_dir / ".queries.partial.npy"
    np.save(temporary_vectors, queries, allow_pickle=False)
    table_rows = ["query_ordinal\tvector_id\trecord_ordinal\ttranscript_id\twindow_offset\tselection_cycle\n"]
    table_rows.extend(
        "{query_ordinal}\t{vector_id}\t{record_ordinal}\t{transcript_id}\t{window_offset}\t{selection_cycle}\n".format(**row)
        for row in selected
    )
    temporary_table = query_dir / ".queries.partial.tsv"
    atomic_text_dump(temporary_table, "".join(table_rows))
    vector_ids = np.asarray([row["vector_id"] for row in selected], dtype=np.int64)
    selected_sha256 = sha256_json(selected)
    query_identity = {
        "schema_version": "ginflow-query-selection-identity-v1",
        "dataset_id": dataset_id,
        "selection_method": "sha256-ranked-records-one-window-per-record-cycle",
        "seed": seed,
        "query_count": int(query_count),
        "dimension": int(queries.shape[1]),
        "embedding_cache_id": flatten_manifest["embedding_cache_id"],
        "query_vectors_sha256": sha256_array(queries),
        "query_ids_sha256": sha256_array(vector_ids),
        "selection_rows_sha256": selected_sha256,
        "queries_table_sha256": sha256_file(temporary_table),
    }
    manifest = {
        "schema_version": "ginflow-query-selection-v1",
        "dataset_id": dataset_id,
        "created_at": utc_now(),
        "selection_method": "sha256-ranked-records-one-window-per-record-cycle",
        "seed": seed,
        "query_count": int(query_count),
        "dimension": int(queries.shape[1]),
        "embedding_cache_id": flatten_manifest["embedding_cache_id"],
        "query_vectors": {"path": vectors_path.name, "sha256": query_identity["query_vectors_sha256"]},
        "query_ids_sha256": query_identity["query_ids_sha256"],
        "queries_table": {"path": table_path.name, "sha256": query_identity["queries_table_sha256"]},
        "identity": query_identity,
    }
    manifest["query_selection_id"] = stable_id("query-selection", query_identity)
    temporary_manifest = query_dir / ".query-selection.partial.json"
    atomic_json_dump(temporary_manifest, manifest)
    os.replace(temporary_vectors, vectors_path)
    os.replace(temporary_table, table_path)
    os.replace(temporary_manifest, manifest_path)
    return manifest


def _merge_topk(
    best_scores: np.ndarray,
    best_ids: np.ndarray,
    scores: np.ndarray,
    block_start: int,
    k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge a score block into a deterministic score-desc/id-asc top-k table."""
    if scores.ndim != 2:
        raise ValueError("score block must be 2-D")
    take = min(k, scores.shape[1])
    if take < 1:
        return best_scores, best_ids
    # ``argpartition`` does not define which equal-valued rows cross its
    # boundary. Resolve that boundary explicitly by lowest vector id before
    # combining the block with the retained candidates. This makes ties stable
    # even when an equal-score group spans database chunks.
    block_ids_full = np.arange(scores.shape[1], dtype=np.int64) + int(block_start)
    selected_positions: list[np.ndarray] = []
    for row in scores:
        partition = np.argpartition(-row, kth=take - 1)[:take]
        threshold = row[partition].min()
        strict = np.flatnonzero(row > threshold)
        needed = take - strict.size
        ties = np.flatnonzero(row == threshold)
        tie_ids = block_ids_full[ties]
        tie_order = np.argsort(tie_ids, kind="stable")[:needed]
        selected_positions.append(np.concatenate((strict, ties[tie_order])))
    positions = np.stack(selected_positions, axis=0)
    block_scores = np.take_along_axis(scores, positions, axis=1)
    block_ids = positions.astype(np.int64, copy=False) + int(block_start)
    candidate_scores = np.concatenate((best_scores, block_scores), axis=1)
    candidate_ids = np.concatenate((best_ids, block_ids), axis=1)
    order = np.lexsort((candidate_ids, -candidate_scores), axis=1)[:, :k]
    return (
        np.take_along_axis(candidate_scores, order, axis=1),
        np.take_along_axis(candidate_ids, order, axis=1),
    )


def _cupy_module(engine: str, gpu_device: int) -> tuple[Any | None, str]:
    if engine not in {"auto", "cpu", "gpu"}:
        raise ValueError("engine must be auto, cpu, or gpu")
    if engine == "cpu":
        return None, "cpu"
    try:
        import cupy as cp

        if int(cp.cuda.runtime.getDeviceCount()) <= gpu_device:
            raise RuntimeError(f"CUDA device {gpu_device} is not visible")
        return cp, "gpu"
    except Exception as exc:
        if engine == "gpu":
            raise ValueError(f"GPU exact ground truth was requested but is unavailable: {exc}") from exc
        return None, "cpu"


def compute_exact_topk(
    vectors: np.ndarray,
    queries: np.ndarray,
    *,
    k: int = 100,
    database_chunk: int = DEFAULT_DATABASE_CHUNK,
    query_chunk: int = DEFAULT_QUERY_CHUNK,
    engine: str = "auto",
    gpu_device: int = 0,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Compute exact cosine top-k with a bounded database block on CPU or GPU.

    The database block is outermost.  GPU mode therefore transfers each block
    once, then compares all query chunks before releasing it; it does not copy a
    30k database into 8 GiB VRAM or retransmit it for every query batch.
    """
    if vectors.ndim != 2 or queries.ndim != 2:
        raise ValueError("vectors and queries must be 2-D")
    if vectors.shape[1] != queries.shape[1]:
        raise ValueError("vectors and queries have incompatible dimensions")
    if vectors.shape[0] < 1 or queries.shape[0] < 1:
        raise ValueError("vectors and queries must be non-empty")
    if not (1 <= k <= vectors.shape[0]):
        raise ValueError("k must be between 1 and the database size")
    if database_chunk < 1 or query_chunk < 1:
        raise ValueError("database_chunk and query_chunk must be >= 1")
    cp, resolved_engine = _cupy_module(engine, gpu_device)
    query_matrix = np.ascontiguousarray(queries, dtype=np.float32)
    n_queries = query_matrix.shape[0]
    best_scores = np.full((n_queries, k), -np.inf, dtype=np.float32)
    best_ids = np.full((n_queries, k), -1, dtype=np.int64)
    query_gpu = None
    if cp is not None:
        with cp.cuda.Device(gpu_device):
            query_gpu = cp.asarray(query_matrix)

    for database_start in range(0, vectors.shape[0], database_chunk):
        database_end = min(vectors.shape[0], database_start + database_chunk)
        # At most one database block is materialized beyond the source memmap.
        database_block = np.ascontiguousarray(vectors[database_start:database_end], dtype=np.float32)
        if cp is not None:
            with cp.cuda.Device(gpu_device):
                database_gpu = cp.asarray(database_block)
                for query_start in range(0, n_queries, query_chunk):
                    query_end = min(n_queries, query_start + query_chunk)
                    scores_gpu = query_gpu[query_start:query_end] @ database_gpu.T
                    # Make the transfer boundary explicit: the returned host
                    # block is a completed GPU GEMM result, not a pending CUDA
                    # operation that could be mixed with the next block.
                    cp.cuda.get_current_stream().synchronize()
                    scores = cp.asnumpy(scores_gpu)
                    merged_scores, merged_ids = _merge_topk(
                        best_scores[query_start:query_end],
                        best_ids[query_start:query_end],
                        np.asarray(scores, dtype=np.float32),
                        database_start,
                        k,
                    )
                    best_scores[query_start:query_end] = merged_scores
                    best_ids[query_start:query_end] = merged_ids
                del database_gpu
        else:
            for query_start in range(0, n_queries, query_chunk):
                query_end = min(n_queries, query_start + query_chunk)
                scores = query_matrix[query_start:query_end] @ database_block.T
                merged_scores, merged_ids = _merge_topk(
                    best_scores[query_start:query_end],
                    best_ids[query_start:query_end],
                    np.asarray(scores, dtype=np.float32),
                    database_start,
                    k,
                )
                best_scores[query_start:query_end] = merged_scores
                best_ids[query_start:query_end] = merged_ids
    return best_scores, best_ids, resolved_engine


def make_ground_truth_cache(
    flat_dir: Path,
    query_dir: Path,
    ground_truth_dir: Path,
    *,
    dataset_id: str,
    k: int = 100,
    database_chunk: int = DEFAULT_DATABASE_CHUNK,
    query_chunk: int = DEFAULT_QUERY_CHUNK,
    engine: str = "auto",
    gpu_device: int = 0,
    overwrite: bool = False,
) -> dict[str, Any]:
    flatten_manifest = _read_json(flat_dir / "flatten-manifest.json")
    selection_manifest = _read_json(query_dir / "query-selection.json")
    vectors = np.load(flat_dir / "vectors.npy", mmap_mode="r", allow_pickle=False)
    queries = np.load(query_dir / "queries.npy", allow_pickle=False)
    ground_truth_dir.mkdir(parents=True, exist_ok=True)
    values_path = ground_truth_dir / "ground-truth.npz"
    manifest_path = ground_truth_dir / "ground-truth.json"
    if any(path.exists() for path in (values_path, manifest_path)) and not overwrite:
        raise ValueError(f"ground-truth cache already exists in {ground_truth_dir}; use --overwrite only after review")
    partials = (
        ground_truth_dir / ".ground-truth.partial.npz",
        ground_truth_dir / ".ground-truth.partial.json",
    )
    if any(path.exists() for path in partials):
        raise ValueError(f"partial ground-truth cache exists in {ground_truth_dir}; inspect it before retrying")
    started = time.perf_counter()
    scores, ids, resolved_engine = compute_exact_topk(
        vectors,
        queries,
        k=k,
        database_chunk=database_chunk,
        query_chunk=query_chunk,
        engine=engine,
        gpu_device=gpu_device,
    )
    elapsed = time.perf_counter() - started
    temporary_values = ground_truth_dir / ".ground-truth.partial.npz"
    np.savez_compressed(temporary_values, ids=ids, scores=scores)
    ids_sha256 = sha256_array(ids)
    ground_truth_identity = {
        "schema_version": "ginflow-exact-ground-truth-identity-v1",
        "dataset_id": dataset_id,
        "k": int(k),
        "metric": "cosine",
        "database_window_count": int(vectors.shape[0]),
        "dimension": int(vectors.shape[1]),
        "embedding_cache_id": flatten_manifest["embedding_cache_id"],
        "query_selection_id": selection_manifest["query_selection_id"],
        "query_ids_sha256": selection_manifest["query_ids_sha256"],
        # Exact neighbours are the semantic ground truth. Scores remain
        # provenance because CPU/GPU BLAS may differ by harmless last bits.
        "ground_truth_ids_sha256": ids_sha256,
    }
    manifest = {
        "schema_version": "ginflow-exact-ground-truth-v1",
        "dataset_id": dataset_id,
        "created_at": utc_now(),
        "k": int(k),
        "metric": "cosine",
        "database_window_count": int(vectors.shape[0]),
        "dimension": int(vectors.shape[1]),
        "embedding_cache_id": flatten_manifest["embedding_cache_id"],
        "query_selection_id": selection_manifest["query_selection_id"],
        "query_ids_sha256": selection_manifest["query_ids_sha256"],
        "ground_truth": {
            "path": values_path.name,
            "ids_sha256": ids_sha256,
            "scores_sha256": sha256_array(scores),
        },
        "engine": {
            "requested": engine,
            "used": resolved_engine,
            "gpu_device": gpu_device if resolved_engine == "gpu" else None,
            "database_chunk": int(database_chunk),
            "query_chunk": int(query_chunk),
            "elapsed_seconds": elapsed,
        },
        "identity": ground_truth_identity,
    }
    manifest["ground_truth_cache_id"] = stable_id("ground-truth", ground_truth_identity)
    temporary_manifest = ground_truth_dir / ".ground-truth.partial.json"
    atomic_json_dump(temporary_manifest, manifest)
    os.replace(temporary_values, values_path)
    os.replace(temporary_manifest, manifest_path)
    return manifest


def _default_dataset_id(input_path: Path | None) -> str:
    if input_path is None:
        raise ValueError("--dataset is required when no input path is supplied")
    return input_path.stem


def _require_flat_paths(root: Path) -> tuple[Path, Path, Path]:
    return root / "flat", root / "queries", root / "ground-truth"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)

    cache = subcommands.add_parser("cache-windows", help="run graph/embed/window Nextflow cache once")
    cache.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    cache.add_argument("--dataset", help="stable cache dataset id; defaults to input filename stem")
    cache.add_argument("--input", type=Path, required=True)
    cache.add_argument("--profile", default="docker,gpu")
    cache.add_argument("--shard-size", type=int, default=50)
    cache.add_argument("--window-size", type=int, default=11)
    cache.add_argument("--window-stride", type=int, default=1)
    cache.add_argument("--nextflow", default="nextflow")
    cache.add_argument("--dry-run", action="store_true")

    flatten = subcommands.add_parser("flatten", help="stream windows into a normalized memmap")
    flatten.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    flatten.add_argument("--dataset", required=True)
    flatten.add_argument("--windows-dir", type=Path)
    flatten.add_argument("--row-chunk", type=int, default=4096)
    flatten.add_argument("--overwrite", action="store_true")

    queries = subcommands.add_parser("select-queries", help="make deterministic record-balanced query cache")
    queries.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    queries.add_argument("--dataset", required=True)
    queries.add_argument("--query-count", type=int, default=DEFAULT_QUERY_COUNT)
    queries.add_argument("--seed", default=DEFAULT_QUERY_SEED)
    queries.add_argument("--overwrite", action="store_true")

    truth = subcommands.add_parser("ground-truth", help="compute bounded-memory exact cosine top-100")
    truth.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    truth.add_argument("--dataset", required=True)
    truth.add_argument("--k", type=int, default=100)
    truth.add_argument("--database-chunk", type=int, default=DEFAULT_DATABASE_CHUNK)
    truth.add_argument("--query-chunk", type=int, default=DEFAULT_QUERY_CHUNK)
    truth.add_argument("--engine", choices=("auto", "cpu", "gpu"), default="auto")
    truth.add_argument("--gpu-device", type=int, default=0)
    truth.add_argument("--overwrite", action="store_true")

    inspect = subcommands.add_parser("inspect", help="print cache provenance and artifact identity")
    inspect.add_argument("--cache-dir", help="external cache root (or $GINFLOW_BENCHMARK_CACHE)")
    inspect.add_argument("--dataset", required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    repo_root = Path(__file__).resolve().parents[1]
    cache_dir = cache_root_from_args(getattr(args, "cache_dir", None))
    try:
        if args.command == "cache-windows":
            dataset_id = args.dataset or _default_dataset_id(args.input)
            result = run_window_cache(
                repo_root=repo_root,
                cache_dir=cache_dir,
                dataset_id=dataset_id,
                input_path=args.input.resolve(),
                profile=args.profile,
                shard_size=args.shard_size,
                window_size=args.window_size,
                window_stride=args.window_stride,
                nextflow=args.nextflow,
                dry_run=args.dry_run,
            )
        else:
            root = dataset_root(cache_dir, args.dataset)
            flat_dir, query_dir, ground_truth_dir = _require_flat_paths(root)
            if args.command == "inspect":
                path = _cache_provenance_path(root)
                if not path.is_file():
                    raise ValueError(f"no provenance file: {path}")
                result = _read_json(path)
            else:
                # A manifest is published last, but the lock also prevents a
                # reader/writer race during an intentional --overwrite.
                with _cache_lock(root):
                    if args.command == "flatten":
                        windows_dir = args.windows_dir or root / "artifacts" / "windows"
                        manifest = flatten_window_shards(
                            windows_dir.resolve(),
                            flat_dir,
                            dataset_id=args.dataset,
                            row_chunk=args.row_chunk,
                            overwrite=args.overwrite,
                        )
                        update_cache_provenance(
                            root,
                            {
                                "dataset_id": args.dataset,
                                "embedding_cache_id": manifest["embedding_cache_id"],
                                "flatten": manifest,
                                "hardware": hardware_snapshot(),
                                "software": software_snapshot(repo_root),
                            },
                        )
                        result = manifest
                    elif args.command == "select-queries":
                        manifest = make_query_cache(
                            flat_dir,
                            query_dir,
                            dataset_id=args.dataset,
                            query_count=args.query_count,
                            seed=args.seed,
                            overwrite=args.overwrite,
                        )
                        update_cache_provenance(root, {"query_selection": manifest})
                        result = manifest
                    elif args.command == "ground-truth":
                        manifest = make_ground_truth_cache(
                            flat_dir,
                            query_dir,
                            ground_truth_dir,
                            dataset_id=args.dataset,
                            k=args.k,
                            database_chunk=args.database_chunk,
                            query_chunk=args.query_chunk,
                            engine=args.engine,
                            gpu_device=args.gpu_device,
                            overwrite=args.overwrite,
                        )
                        update_cache_provenance(root, {"ground_truth": manifest})
                        result = manifest
                    else:  # pragma: no cover - argparse enforces subcommands
                        raise AssertionError(args.command)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
