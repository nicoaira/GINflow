#!/usr/bin/env python3
"""GPU CAGRA companion search for the HNSWLIB workflow.

The CPU HNSWLIB path indexes centroid-code windows with a registered custom
distance.  cuVS 24.10 cannot register that distance or persist float16 CAGRA
datasets, so this optional accelerator indexes the normalized original window
vectors after one global int8 scale.  CAGRA supplies candidates and the final
scores are recomputed from the preserved residue embeddings.
"""
from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hnswlib_index import quantization_dir
from record_pack import pack_records, write_records
from search_compact_hnswlib import load_embeddings, rerank_candidates


INDEX_TYPE = "HNSWLIB_GPU_CAGRA"
WINDOW_SUFFIX = ".windows.npz"
MANIFEST_SUFFIX = ".windows.manifest.json"
DEFAULT_INT8_SCALE = 850.0
COMPAT_KEYS = (
    ("window_size", "window_size"),
    ("stride", "window_stride"),
    ("window_dim", "window_dim"),
    ("checkpoint_sha256", "checkpoint_sha256"),
)


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def pair_window_shards(
    windows_dir: Path, manifests_dir: Path | None = None
) -> list[tuple[Path, Path]]:
    manifests_dir = manifests_dir or windows_dir
    windows = {path.name[: -len(WINDOW_SUFFIX)]: path for path in windows_dir.glob(f"*{WINDOW_SUFFIX}")}
    manifests = {
        path.name[: -len(MANIFEST_SUFFIX)]: path
        for path in manifests_dir.glob(f"*{MANIFEST_SUFFIX}")
    }
    only_windows = sorted(set(windows) - set(manifests))
    only_manifests = sorted(set(manifests) - set(windows))
    if only_windows or only_manifests:
        raise ValueError(
            "window / manifest shard names do not match: "
            f"only-windows={only_windows} only-manifests={only_manifests}"
        )
    if not windows:
        raise ValueError(f"no original window shards found in {windows_dir}")
    return [(windows[key], manifests[key]) for key in sorted(windows)]


def compatible_tuple(manifest: dict[str, Any]) -> tuple[Any, ...]:
    return tuple(manifest.get(key) for key, _db_key in COMPAT_KEYS)


def window_count(manifest: dict[str, Any]) -> int:
    return int(sum(int(record.get("n_windows", 0)) for record in manifest.get("records", [])))


def read_window_rows(
    window_path: Path, manifest: dict[str, Any]
) -> tuple[np.ndarray, list[tuple[str, int, int]]]:
    window_size = int(manifest["window_size"])
    window_dim = int(manifest["window_dim"])
    stride = int(manifest["stride"])
    values: list[np.ndarray] = []
    mapping: list[tuple[str, int, int]] = []
    with np.load(window_path) as arrays:
        for record in manifest.get("records", []):
            identifier = str(record["identifier"])
            if identifier not in arrays.files:
                raise KeyError(f"{identifier} is in {manifest} but missing from {window_path}")
            rows = np.asarray(arrays[identifier], dtype=np.float32)
            expected = int(record.get("n_windows", rows.shape[0]))
            if rows.ndim != 2 or rows.shape[1] != window_dim or rows.shape[0] != expected:
                raise ValueError(
                    f"{identifier} in {window_path} has shape {rows.shape}; "
                    f"expected ({expected}, {window_dim})"
                )
            if rows.size and not np.isfinite(rows).all():
                raise ValueError(f"{identifier} in {window_path} contains non-finite values")
            values.append(np.ascontiguousarray(rows, dtype=np.float32))
            mapping.extend(
                (identifier, offset * stride, offset * stride + window_size)
                for offset in range(rows.shape[0])
            )
    if not values:
        return np.empty((0, window_dim), dtype=np.float32), mapping
    return np.ascontiguousarray(np.concatenate(values, axis=0), dtype=np.float32), mapping


def copy_quantization(source: Path, destination: Path) -> None:
    source_dir = quantization_dir(source)
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("centroids.npy", "similarity.npy", "quantization.json"):
        source_file = source_dir / name
        if not source_file.exists():
            raise FileNotFoundError(f"missing quantization artifact {source_file}")
        shutil.copy2(source_file, destination / name)


def write_mapping(path: Path, mapping: Iterable[tuple[str, int, int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        writer.writerows(
            (row_id, identifier, start, end)
            for row_id, (identifier, start, end) in enumerate(mapping)
        )


def resolve_scale(max_abs: float, requested: float | None) -> float:
    if not np.isfinite(max_abs) or max_abs <= 0:
        raise ValueError(f"window values have an invalid maximum absolute value: {max_abs}")
    scale = (
        float(requested)
        if requested is not None
        else min(DEFAULT_INT8_SCALE, 126.0 / max_abs)
    )
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("int8 scale must be finite and positive")
    if scale * max_abs > 127.0 + 1e-5:
        raise ValueError(
            f"int8 scale {scale:g} clips a window maximum of {max_abs:g}; "
            "use a smaller scale"
        )
    return scale


def quantize_windows(values: np.ndarray, scale: float) -> np.ndarray:
    rows = np.asarray(values, dtype=np.float32)
    if rows.ndim != 2:
        raise ValueError(f"windows must be a 2-D array, got {rows.shape}")
    scaled = np.rint(rows * np.float32(scale))
    if scaled.size and (scaled.min() < -127 or scaled.max() > 127):
        raise ValueError("scaled windows exceed the signed int8 range")
    return np.ascontiguousarray(scaled, dtype=np.int8)


def require_gpu():
    try:
        import cupy as cp
        from cuvs.neighbors import cagra
    except ImportError as exc:  # pragma: no cover - GPU environment specific
        raise ValueError(
            "hnswlib GPU mode requires CuPy and cuVS 24.10; use -profile gpu"
        ) from exc
    if int(cp.cuda.runtime.getDeviceCount()) < 1:
        raise ValueError("hnswlib GPU mode requires a visible NVIDIA GPU")
    return cp, cagra


def build_database(
    windows_dir: Path,
    quantization: Path,
    outdir: Path,
    *,
    intermediate_graph_degree: int,
    graph_degree: int,
    build_algo: str,
    int8_scale: float | None,
    itopk_size: int,
    search_batch_size: int,
    embeddings: list[Path],
    graph_metadata: list[Path],
    manifests_dir: Path | None = None,
) -> dict[str, Any]:
    if intermediate_graph_degree < 1 or graph_degree < 1:
        raise ValueError("CAGRA graph degrees must be positive")
    if graph_degree > intermediate_graph_degree:
        raise ValueError("CAGRA graph degree cannot exceed intermediate graph degree")
    if itopk_size < 1 or search_batch_size < 1:
        raise ValueError("GPU search parameters must be positive")
    shards = pair_window_shards(windows_dir, manifests_dir)
    manifests = [load_json(manifest_path) for _window_path, manifest_path in shards]
    reference = manifests[0]
    expected = compatible_tuple(reference)
    for manifest_path, manifest in zip((item[1] for item in shards), manifests):
        if compatible_tuple(manifest) != expected:
            raise ValueError(f"{manifest_path} is incompatible with the first window manifest")
    total_windows = sum(window_count(manifest) for manifest in manifests)
    if total_windows < 1:
        raise ValueError("no original windows to index")
    window_dim = int(reference["window_dim"])

    max_abs = 0.0
    for window_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        values, _mapping = read_window_rows(window_path, manifest)
        if values.size:
            max_abs = max(max_abs, float(np.max(np.abs(values))))
    scale = resolve_scale(max_abs, int8_scale)

    cp, cagra = require_gpu()
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="hnswlib-gpu-", dir=outdir) as work:
        data_path = Path(work) / "windows.int8.npy"
        data = np.lib.format.open_memmap(
            data_path,
            mode="w+",
            dtype=np.int8,
            shape=(total_windows, window_dim),
        )
        mapping: list[tuple[str, int, int]] = []
        base = 0
        for window_path, manifest_path in shards:
            manifest = load_json(manifest_path)
            values, local_mapping = read_window_rows(window_path, manifest)
            quantized = quantize_windows(values, scale)
            stop = base + quantized.shape[0]
            data[base:stop] = quantized
            mapping.extend(local_mapping)
            base = stop
        data.flush()
        del data
        host_data = np.load(data_path, mmap_mode="r")
        device_data = cp.asarray(host_data)
        params = cagra.IndexParams(
            metric="sqeuclidean",
            intermediate_graph_degree=int(intermediate_graph_degree),
            graph_degree=int(graph_degree),
            build_algo=str(build_algo),
        )
        index = cagra.build(params, device_data)
        cp.cuda.Stream.null.synchronize()
        index_dir = outdir / "cagra"
        index_dir.mkdir(parents=True, exist_ok=True)
        cagra.save(str(index_dir / "index.bin"), index, include_dataset=True)
        cp.cuda.Stream.null.synchronize()

    write_mapping(outdir / "windows.tsv", mapping)
    copy_quantization(quantization, outdir / "quantization")
    packed_embeddings, records = pack_records(embeddings, graph_metadata)
    write_records(outdir, packed_embeddings, records)
    first = next(iter(packed_embeddings.values()), None)
    meta: dict[str, Any] = {
        "backend": "hnswlib",
        "index_type": INDEX_TYPE,
        "index_class": "cuvs.neighbors.cagra.Index",
        "space": "sqeuclidean_int8",
        "metric": "original_window_cosine_after_exact_rerank",
        "candidate_representation": "original_normalized_window_int8",
        "quantized_nodes": False,
        "node_quantization_available": True,
        "node_quantization_candidate_only": False,
        "quantization": "quantization",
        "gpu": True,
        "gpu_index": "cagra",
        "gpu_index_path": "cagra/index.bin",
        "gpu_metric": "sqeuclidean",
        "gpu_data_dtype": "int8",
        "gpu_int8_scale": scale,
        "gpu_intermediate_graph_degree": int(intermediate_graph_degree),
        "gpu_graph_degree": int(graph_degree),
        "gpu_build_algo": str(build_algo),
        "gpu_itopk_size": int(itopk_size),
        "gpu_search_batch_size": int(search_batch_size),
        "cuvs_version": importlib.metadata.version("cuvs"),
        "window_size": int(reference["window_size"]),
        "window_stride": int(reference["stride"]),
        "embedding_dim": int(reference["embedding_dim"]),
        "window_dim": window_dim,
        "l2_normalized": True,
        "score_scale": 1.0,
        "n_records": int(sum(len(manifest.get("records", [])) + len(manifest.get("skipped_short", [])) for manifest in manifests)),
        "n_windows": int(total_windows),
        "n_skipped_short": int(sum(len(manifest.get("skipped_short", [])) for manifest in manifests)),
        "ginfinity_version": reference.get("ginfinity_version"),
        "model_version": reference.get("model_version"),
        "checkpoint_sha256": reference.get("checkpoint_sha256"),
        "original_embeddings_preserved": True,
        "original_embedding_dim": int(reference["embedding_dim"]),
        "original_embedding_dtype": str(first.dtype) if first is not None else None,
        "has_residue_embeddings": True,
        "n_packed_records": len(records),
        "quantization_metadata": load_json(outdir / "quantization" / "quantization.json"),
    }
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def load_targets(path: Path) -> list[tuple[str, int, int]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"faiss_id", "transcript_id", "start", "end"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must have columns {sorted(required)}")
        rows: list[tuple[str, int, int] | None] = []
        for record in reader:
            faiss_id = int(record["faiss_id"])
            while len(rows) <= faiss_id:
                rows.append(None)
            rows[faiss_id] = (
                record["transcript_id"],
                int(record["start"]),
                int(record["end"]),
            )
    if any(row is None for row in rows):
        raise ValueError("windows.tsv has gaps in faiss_id")
    return [row for row in rows if row is not None]


def check_compatible(manifest: dict[str, Any], meta: dict[str, Any]) -> None:
    mismatches = [
        f"{query_key}={manifest.get(query_key)!r} vs {db_key}={meta.get(db_key)!r}"
        for query_key, db_key in COMPAT_KEYS
        if manifest.get(query_key) != meta.get(db_key)
    ]
    if mismatches:
        raise ValueError("query windows are incompatible with the GPU HNSWLIB database: " + "; ".join(mismatches))


def search_database(
    query_windows_dir: Path,
    database: Path,
    *,
    k: int,
    candidate_k: int | None,
    query_embeddings_paths: list[Path],
    min_similarity: float,
    itopk_size: int | None,
    search_batch_size: int | None,
    manifests_dir: Path | None = None,
) -> list[tuple]:
    if k < 1:
        raise ValueError("--k must be >= 1")
    meta = load_json(database / "meta.json")
    if str(meta.get("index_type", "")).upper() != INDEX_TYPE:
        raise ValueError(f"{database} is not a GPU CAGRA HNSWLIB database")
    targets = load_targets(database / "windows.tsv")
    requested_k = int(candidate_k or k)
    if requested_k < k:
        raise ValueError("--candidate-k must be >= --k")
    search_k = min(requested_k, len(targets))
    beam = min(
        max(search_k, int(itopk_size or meta.get("gpu_itopk_size", 256))),
        len(targets),
    )
    batch_size = int(search_batch_size or meta.get("gpu_search_batch_size", 512))
    if batch_size < 1:
        raise ValueError("GPU search batch size must be positive")
    scale = float(meta["gpu_int8_scale"])
    database_embeddings = load_embeddings([database / "embeddings.npz"])
    query_embeddings = load_embeddings(query_embeddings_paths)
    cp, cagra = require_gpu()
    index = cagra.load(str(database / "cagra" / "index.bin"))
    cp.cuda.Stream.null.synchronize()
    search_params = cagra.SearchParams(itopk_size=beam, max_queries=batch_size)
    hits: list[tuple] = []

    for window_path, manifest_path in pair_window_shards(query_windows_dir, manifests_dir):
        manifest = load_json(manifest_path)
        check_compatible(manifest, meta)
        values, mapping = read_window_rows(window_path, manifest)
        if values.shape[0] == 0:
            continue
        label_parts: list[np.ndarray] = []
        for start in range(0, values.shape[0], batch_size):
            stop = min(start + batch_size, values.shape[0])
            query_codes = cp.asarray(quantize_windows(values[start:stop], scale))
            _distances, labels = cagra.search(
                search_params,
                index,
                query_codes,
                search_k,
            )
            cp.cuda.Stream.null.synchronize()
            label_parts.append(np.asarray(cp.asnumpy(labels), dtype=np.uint64))
        labels = np.concatenate(label_parts, axis=0)
        final_labels, final_scores = rerank_candidates(
            labels,
            mapping,
            targets,
            query_embeddings,
            database_embeddings,
            k,
        )
        stride = int(manifest["stride"])
        for (identifier, query_start, query_end), row_labels, row_scores in zip(
            mapping, final_labels, final_scores
        ):
            rank = 0
            for target_label, score in zip(row_labels, row_scores):
                if int(target_label) == np.iinfo(np.uint64).max:
                    continue
                value = float(score)
                if value < min_similarity:
                    continue
                rank += 1
                target_id, target_start, target_end = targets[int(target_label)]
                hits.append(
                    (
                        identifier,
                        query_start,
                        query_end,
                        target_id,
                        target_start,
                        target_end,
                        value,
                        rank,
                    )
                )
                if rank >= k:
                    break
    return hits


def write_seeds(path: Path, hits: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "query_id",
                "query_start",
                "query_end",
                "target_id",
                "target_start",
                "target_end",
                "score",
                "rank",
            ]
        )
        for row in hits:
            writer.writerow([row[0], row[1], row[2], row[3], row[4], row[5], f"{row[6]:.6f}", row[7]])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build")
    build.add_argument("--windows-dir", type=Path, required=True)
    build.add_argument("--quantization", type=Path, required=True)
    build.add_argument("--embeddings", type=Path, nargs="*", required=True)
    build.add_argument("--graph-metadata", type=Path, nargs="*", required=True)
    build.add_argument("--outdir", type=Path, required=True)
    build.add_argument("--intermediate-graph-degree", type=int, default=128)
    build.add_argument("--graph-degree", type=int, default=64)
    build.add_argument("--build-algo", default="nn_descent")
    build.add_argument("--int8-scale", type=float)
    build.add_argument("--itopk-size", type=int, default=256)
    build.add_argument("--search-batch-size", type=int, default=512)
    build.add_argument("--manifests-dir", type=Path)
    search = subparsers.add_parser("search")
    search.add_argument("--windows-dir", type=Path, required=True)
    search.add_argument("--manifests-dir", type=Path)
    search.add_argument("--database", type=Path, required=True)
    search.add_argument("--query-embeddings", type=Path, nargs="*", required=True)
    search.add_argument("--output", type=Path, required=True)
    search.add_argument("--k", type=int, default=50)
    search.add_argument("--candidate-k", type=int)
    search.add_argument("--min-similarity", type=float, default=0.8)
    search.add_argument("--itopk-size", type=int)
    search.add_argument("--search-batch-size", type=int)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "build":
            meta = build_database(
                args.windows_dir,
                args.quantization,
                args.outdir,
                intermediate_graph_degree=args.intermediate_graph_degree,
                graph_degree=args.graph_degree,
                build_algo=args.build_algo,
                int8_scale=args.int8_scale,
                itopk_size=args.itopk_size,
                search_batch_size=args.search_batch_size,
                embeddings=args.embeddings,
                graph_metadata=args.graph_metadata,
                manifests_dir=args.manifests_dir,
            )
            print(json.dumps({"outdir": str(args.outdir), "n_windows": meta["n_windows"], "index_type": meta["index_type"]}))
        else:
            hits = search_database(
                args.windows_dir,
                args.database,
                k=args.k,
                candidate_k=args.candidate_k,
                query_embeddings_paths=args.query_embeddings,
                min_similarity=args.min_similarity,
                itopk_size=args.itopk_size,
                search_batch_size=args.search_batch_size,
                manifests_dir=args.manifests_dir,
            )
            write_seeds(args.output, hits)
            print(json.dumps({"output": str(args.output), "n_seeds": len(hits)}))
    except (OSError, KeyError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
