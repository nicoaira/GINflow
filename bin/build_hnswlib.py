#!/usr/bin/env python3
"""Build an HNSWLIB index over centroid-code windows."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hnswlib_index import create_index, encode_code_windows, load_quantization, quantization_dir
from record_pack import pack_records, write_records


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def pair_window_shards(windows_dir: Path) -> list[tuple[Path, Path]]:
    window_paths = sorted(windows_dir.glob("*.windows.npz"))
    manifest_paths = sorted(windows_dir.glob("*.windows.manifest.json"))
    window_map = {path.name[: -len(".windows.npz")]: path for path in window_paths}
    manifest_map = {path.name[: -len(".windows.manifest.json")]: path for path in manifest_paths}
    only_windows = sorted(set(window_map) - set(manifest_map))
    only_manifests = sorted(set(manifest_map) - set(window_map))
    if only_windows or only_manifests:
        raise ValueError(
            "window / manifest shard names do not match: "
            f"only-windows={only_windows} only-manifests={only_manifests}"
        )
    if not window_map:
        raise ValueError(f"no quantized windows found in {windows_dir}")
    return [(window_map[key], manifest_map[key]) for key in sorted(window_map)]


def compatible_tuple(manifest: dict) -> tuple:
    return (
        manifest.get("window_size"),
        manifest.get("stride"),
        manifest.get("window_dim"),
        manifest.get("checkpoint_sha256"),
        manifest.get("n_centroids"),
        manifest.get("embedding_dim"),
    )


def window_count(manifest: dict) -> int:
    return int(sum(int(record.get("n_windows", 0)) for record in manifest.get("records", [])))


def read_code_rows(
    path: Path, manifest: dict
) -> tuple[np.ndarray, list[tuple[str, int, int]]]:
    window_size = int(manifest["window_size"])
    codes: list[np.ndarray] = []
    mapping: list[tuple[str, int, int]] = []
    stride = int(manifest["stride"])
    with np.load(path) as arrays:
        for record in manifest.get("records", []):
            identifier = str(record["identifier"])
            if identifier not in arrays.files:
                raise KeyError(f"{identifier} is in {manifest} but missing from {path}")
            values = np.asarray(arrays[identifier])
            if values.ndim != 2 or values.shape[1] != window_size:
                raise ValueError(f"{identifier} in {path} has invalid code-window shape {values.shape}")
            if values.shape[0] != int(record.get("n_windows", values.shape[0])):
                raise ValueError(f"{identifier} in {path} disagrees with its manifest window count")
            codes.append(np.ascontiguousarray(values))
            for offset in range(values.shape[0]):
                start = offset * stride
                mapping.append((identifier, start, start + window_size))
    if not codes:
        return np.empty((0, window_size), dtype=np.uint16), []
    return np.ascontiguousarray(np.concatenate(codes, axis=0)), mapping


def copy_quantization(source: Path, destination: Path) -> None:
    source_dir = quantization_dir(source)
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("centroids.npy", "similarity.npy", "quantization.json"):
        source_file = source_dir / name
        if not source_file.exists():
            raise FileNotFoundError(f"missing quantization artifact {source_file}")
        shutil.copy2(source_file, destination / name)


def write_mapping(path: Path, mapping: list[tuple[int, str, int, int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        writer.writerows(mapping)


def build_database(
    windows_dir: Path,
    quantization: Path,
    outdir: Path,
    *,
    m: int,
    ef_construction: int,
    ef_search: int,
    random_seed: int,
    num_threads: int,
    embeddings: list[Path] | None = None,
    graph_metadata: list[Path] | None = None,
) -> dict:
    shards = pair_window_shards(windows_dir)
    manifests = [load_json(manifest_path) for _window_path, manifest_path in shards]
    reference = manifests[0]
    expected = compatible_tuple(reference)
    if not reference.get("quantized_windows"):
        raise ValueError("HNSWLIB requires quantized window shards")
    for manifest_path, manifest in zip((item[1] for item in shards), manifests):
        if compatible_tuple(manifest) != expected:
            raise ValueError(f"{manifest_path} is incompatible with the first quantized window manifest")
    centroids, similarity, quantization_metadata = load_quantization(quantization)
    if int(reference["n_centroids"]) != centroids.shape[0]:
        raise ValueError("quantized window n_centroids does not match the centroid file")
    if int(reference["embedding_dim"]) != centroids.shape[1]:
        raise ValueError("quantized window embedding_dim does not match the centroid file")

    total_windows = sum(window_count(manifest) for manifest in manifests)
    if total_windows < 1:
        raise ValueError("no quantized windows to index")
    dimension = int(reference["window_size"]) * int(reference["embedding_dim"])
    index = create_index(
        dimension,
        total_windows,
        m,
        ef_construction,
        ef_search,
        random_seed,
        num_threads,
    )
    mapping: list[tuple[int, str, int, int]] = []
    base = 0
    for window_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        codes, local_mapping = read_code_rows(window_path, manifest)
        if codes.shape[0] == 0:
            continue
        vectors = encode_code_windows(codes, centroids)
        ids = np.arange(base, base + vectors.shape[0], dtype=np.int64)
        index.add_items(vectors, ids)
        mapping.extend((int(row_id), identifier, start, end) for row_id, (identifier, start, end) in zip(ids, local_mapping))
        base += int(vectors.shape[0])
    if base != total_windows or len(mapping) != total_windows:
        raise ValueError(f"indexed {base} windows but expected {total_windows}")

    outdir.mkdir(parents=True, exist_ok=True)
    index.save_index(str(outdir / "index.bin"))
    write_mapping(outdir / "windows.tsv", mapping)
    copy_quantization(quantization, outdir / "quantization")

    meta: dict[str, Any] = {
        "backend": "hnswlib",
        "index_type": "HNSWLIB",
        "index_class": "hnswlib.Index",
        "space": "ip",
        "metric": "sum_positionwise_centroid_similarity",
        "raw_metric": "sum_positionwise_centroid_similarity",
        "candidate_representation": "centroid_codes",
        "quantized_nodes": True,
        "quantization": "quantization",
        "n_centroids": int(centroids.shape[0]),
        "centroid_dim": int(centroids.shape[1]),
        "centroid_dtype": str(centroids.dtype),
        "similarity_dtype": str(similarity.dtype),
        "window_size": int(reference["window_size"]),
        "window_stride": int(reference["stride"]),
        "embedding_dim": int(reference["embedding_dim"]),
        "window_dim": dimension,
        "l2_normalized": False,
        "score_scale": 1.0,
        "n_records": int(sum(len(manifest.get("records", [])) + len(manifest.get("skipped_short", [])) for manifest in manifests)),
        "n_windows": int(total_windows),
        "n_skipped_short": int(sum(len(manifest.get("skipped_short", [])) for manifest in manifests)),
        "ginfinity_version": reference.get("ginfinity_version"),
        "model_version": reference.get("model_version"),
        "checkpoint_sha256": reference.get("checkpoint_sha256"),
        "hnsw_m": int(m),
        "hnsw_ef_construction": int(ef_construction),
        "hnsw_ef_search": int(ef_search),
        "hnsw_random_seed": int(random_seed),
        "hnsw_num_threads": int(num_threads),
        "original_embeddings_preserved": bool(embeddings and graph_metadata),
        "original_embedding_dim": int(reference["embedding_dim"]),
        "original_embedding_dtype": None,
    }
    if embeddings or graph_metadata:
        if not embeddings or not graph_metadata:
            raise ValueError("--embeddings and --graph-metadata must be passed together")
        packed_embeddings, records = pack_records(embeddings, graph_metadata)
        write_records(outdir, packed_embeddings, records)
        meta["has_residue_embeddings"] = True
        meta["n_packed_records"] = len(records)
        first = next(iter(packed_embeddings.values()), None)
        meta["original_embedding_dtype"] = str(first.dtype) if first is not None else None
    else:
        meta["has_residue_embeddings"] = False
    meta["quantization_metadata"] = quantization_metadata
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--quantization", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, nargs="*")
    parser.add_argument("--graph-metadata", type=Path, nargs="*")
    parser.add_argument("--m", type=int, default=32)
    parser.add_argument("--ef-construction", type=int, default=200)
    parser.add_argument("--ef-search", type=int, default=100)
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--num-threads", type=int, default=0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        meta = build_database(
            args.windows_dir,
            args.quantization,
            args.outdir,
            m=args.m,
            ef_construction=args.ef_construction,
            ef_search=args.ef_search,
            random_seed=args.random_seed,
            num_threads=args.num_threads,
            embeddings=args.embeddings,
            graph_metadata=args.graph_metadata,
        )
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"outdir": str(args.outdir), **{key: meta[key] for key in ("n_windows", "n_centroids", "hnsw_m", "hnsw_ef_search")}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
