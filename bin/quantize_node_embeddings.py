#!/usr/bin/env python3
"""Fit/apply cosine node quantization and persist centroid similarity data."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterator

import numpy as np


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def shard_prefix(path: Path, suffix: str) -> str:
    if not path.name.endswith(suffix):
        raise ValueError(f"{path} does not end with {suffix}")
    return path.name[: -len(suffix)]


def pair_shards(embeddings: list[Path], manifests: list[Path]) -> list[tuple[Path, Path]]:
    embedding_map = {shard_prefix(path, ".npz"): path for path in embeddings}
    manifest_map = {shard_prefix(path, ".manifest.json"): path for path in manifests}
    missing_embeddings = sorted(set(manifest_map) - set(embedding_map))
    missing_manifests = sorted(set(embedding_map) - set(manifest_map))
    if missing_embeddings or missing_manifests:
        raise ValueError(
            "embedding / manifest shard names do not match: "
            f"only-embeddings={missing_manifests} only-manifests={missing_embeddings}"
        )
    return [(embedding_map[key], manifest_map[key]) for key in sorted(embedding_map)]


def normalize_rows(vectors: np.ndarray) -> np.ndarray:
    values = np.ascontiguousarray(vectors, dtype=np.float32)
    if values.ndim != 2:
        raise ValueError(f"node embeddings must be 2-D, got {values.shape}")
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.maximum(norms, np.float32(1e-12))


def iter_nodes(shards: list[tuple[Path, Path]]) -> Iterator[tuple[Path, dict, str, np.ndarray]]:
    for embedding_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        with np.load(embedding_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                if identifier not in arrays.files:
                    raise KeyError(f"{identifier} is in {manifest_path} but missing from {embedding_path}")
                values = normalize_rows(np.asarray(arrays[identifier]))
                yield embedding_path, manifest, identifier, values


def sample_training_nodes(
    shards: list[tuple[Path, Path]], sample_size: int, seed: int
) -> tuple[np.ndarray, int]:
    if sample_size < 1:
        raise ValueError("--sample-size must be >= 1")
    rng = np.random.default_rng(seed)
    reservoir: np.ndarray | None = None
    seen = 0
    dimension: int | None = None
    for _path, _manifest, _identifier, values in iter_nodes(shards):
        if dimension is None:
            dimension = int(values.shape[1])
            reservoir = np.empty((sample_size, dimension), dtype=np.float32)
        elif values.shape[1] != dimension:
            raise ValueError(f"inconsistent node embedding dimensions: {values.shape[1]} vs {dimension}")
        assert reservoir is not None
        for vector in values:
            seen += 1
            if seen <= sample_size:
                reservoir[seen - 1] = vector
                continue
            replacement = int(rng.integers(0, seen))
            if replacement < sample_size:
                reservoir[replacement] = vector
    if seen == 0 or reservoir is None or dimension is None:
        raise ValueError("no node embeddings were provided")
    return np.ascontiguousarray(reservoir[: min(seen, sample_size)]), seen


def _numpy_kmeans(vectors: np.ndarray, k: int, niter: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    initial = rng.choice(vectors.shape[0], size=k, replace=False)
    centroids = np.array(vectors[initial], dtype=np.float32, copy=True)
    for _ in range(niter):
        labels = np.argmax(vectors @ centroids.T, axis=1)
        sums = np.zeros_like(centroids)
        np.add.at(sums, labels, vectors)
        counts = np.bincount(labels, minlength=k)
        empty = counts == 0
        nonempty = ~empty
        centroids[nonempty] = sums[nonempty] / counts[nonempty, None]
        if np.any(empty):
            centroids[empty] = vectors[rng.choice(vectors.shape[0], size=int(empty.sum()), replace=False)]
        centroids = normalize_rows(centroids)
    return centroids


def fit_centroids(vectors: np.ndarray, requested_k: int, niter: int, seed: int) -> tuple[np.ndarray, int]:
    if requested_k < 1:
        raise ValueError("--k must be >= 1")
    if niter < 1:
        raise ValueError("--niter must be >= 1")
    if vectors.ndim != 2 or vectors.shape[0] == 0:
        raise ValueError("k-means needs a non-empty 2-D training matrix")
    k = min(requested_k, int(vectors.shape[0]))
    training = normalize_rows(vectors)
    try:
        import faiss

        kmeans = faiss.Kmeans(
            int(training.shape[1]),
            k,
            niter=niter,
            nredo=1,
            verbose=False,
            spherical=True,
            seed=seed,
        )
        kmeans.train(np.ascontiguousarray(training, dtype=np.float32))
        centroids = np.asarray(kmeans.centroids, dtype=np.float32).reshape(k, training.shape[1])
        centroids = normalize_rows(centroids)
    except (ImportError, RuntimeError, ValueError):
        centroids = _numpy_kmeans(training, k, niter, seed)
    return np.ascontiguousarray(centroids, dtype=np.float16), k


def centroid_similarity(centroids: np.ndarray) -> np.ndarray:
    values = np.ascontiguousarray(centroids, dtype=np.float32)
    if values.ndim != 2:
        raise ValueError(f"centroids must be 2-D, got {values.shape}")
    similarity = values @ values.T
    similarity = (similarity + similarity.T) / np.float32(2.0)
    return np.clip(similarity, np.float32(-1.0), np.float32(1.0)).astype(np.float32)


def assign_nodes(values: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    normalized = normalize_rows(values)
    codebook = np.ascontiguousarray(centroids, dtype=np.float32)
    if normalized.shape[1] != codebook.shape[1]:
        raise ValueError(f"node dimension {normalized.shape[1]} does not match centroid dimension {codebook.shape[1]}")
    labels = np.argmax(normalized @ codebook.T, axis=1)
    dtype = np.uint16 if codebook.shape[0] <= np.iinfo(np.uint16).max else np.uint32
    return np.asarray(labels, dtype=dtype)


def _quantization_metadata(
    centroids: np.ndarray,
    requested_k: int,
    source_shards: list[tuple[Path, Path]],
    *,
    sample_size: int | None = None,
    niter: int | None = None,
    seed: int | None = None,
    total_nodes: int | None = None,
    fitted: bool,
) -> dict:
    reference_manifest = load_json(source_shards[0][1]) if source_shards else {}
    return {
        "format_version": 1,
        "quantized_nodes": True,
        "metric": "cosine",
        "normalization": "per_node_l2",
        "requested_k": int(requested_k),
        "n_centroids": int(centroids.shape[0]),
        "embedding_dim": int(centroids.shape[1]),
        "centroid_dtype": "float16",
        "similarity_dtype": "float32",
        "similarity_file": "similarity.npy",
        "centroid_file": "centroids.npy",
        "fitted_here": fitted,
        "sample_size": sample_size,
        "niter": niter,
        "seed": seed,
        "total_nodes_seen": total_nodes,
        "source_shards": [path.name for path, _manifest in source_shards],
        "ginfinity_version": reference_manifest.get("ginfinity_version"),
        "model_version": reference_manifest.get("model_version"),
        "checkpoint_sha256": reference_manifest.get("checkpoint_sha256"),
    }


def write_quantization(
    shards: list[tuple[Path, Path]],
    outdir: Path,
    centroids: np.ndarray,
    metadata: dict,
) -> dict:
    metadata = {
        "format_version": 1,
        "quantized_nodes": True,
        "metric": "cosine",
        "normalization": "per_node_l2",
        "n_centroids": int(centroids.shape[0]),
        "embedding_dim": int(centroids.shape[1]),
        "centroid_dtype": "float16",
        "similarity_dtype": "float32",
        "similarity_file": "similarity.npy",
        "centroid_file": "centroids.npy",
        **metadata,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    nodes_dir = outdir / "nodes"
    nodes_dir.mkdir(parents=True, exist_ok=True)
    np.save(outdir / "centroids.npy", np.ascontiguousarray(centroids, dtype=np.float16))
    np.save(outdir / "similarity.npy", centroid_similarity(centroids))

    records_written = 0
    for embedding_path, manifest_path in shards:
        source_manifest = load_json(manifest_path)
        quantized: dict[str, np.ndarray] = {}
        records = []
        with np.load(embedding_path) as arrays:
            for record in source_manifest.get("records", []):
                identifier = str(record["identifier"])
                if identifier not in arrays.files:
                    raise KeyError(f"{identifier} is in {manifest_path} but missing from {embedding_path}")
                values = np.asarray(arrays[identifier])
                codes = assign_nodes(values, centroids)
                quantized[identifier] = codes
                records.append(
                    {
                        "identifier": identifier,
                        "length": int(record.get("core_length", values.shape[0])),
                        "n_nodes": int(values.shape[0]),
                        "shape": [int(codes.shape[0])],
                    }
                )
                records_written += 1

        prefix = shard_prefix(embedding_path, ".npz")
        output_npz = nodes_dir / f"{prefix}.quantized.npz"
        output_manifest = nodes_dir / f"{prefix}.quantized.manifest.json"
        np.savez_compressed(output_npz, **quantized)
        output_manifest.write_text(
            json.dumps(
                {
                    **metadata,
                    "status": "complete",
                    "input": embedding_path.name,
                    "input_manifest": manifest_path.name,
                    "records": records,
                },
                indent=2,
            )
            + "\n"
        )

    metadata = {**metadata, "status": "complete", "n_shards": len(shards), "n_records": records_written}
    (outdir / "quantization.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def load_existing_quantization(path: Path) -> tuple[np.ndarray, dict]:
    source = path
    if source.is_file():
        source = source.parent
    centroids_path = source / "centroids.npy"
    metadata_path = source / "quantization.json"
    if not centroids_path.exists():
        raise FileNotFoundError(f"missing centroid file: {centroids_path}")
    centroids = np.load(centroids_path)
    if centroids.dtype != np.float16:
        raise ValueError(f"{centroids_path} must contain float16 centroids, got {centroids.dtype}")
    if centroids.ndim != 2 or centroids.shape[0] == 0 or centroids.shape[1] == 0:
        raise ValueError(f"invalid centroid shape {centroids.shape}")
    metadata = load_json(metadata_path) if metadata_path.exists() else {}
    if int(metadata.get("n_centroids", centroids.shape[0])) != centroids.shape[0]:
        raise ValueError("quantization metadata n_centroids does not match centroids.npy")
    return np.ascontiguousarray(centroids), metadata


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, nargs="+", required=True)
    parser.add_argument("--manifests", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--centroids-dir", type=Path)
    parser.add_argument("--k", type=int, default=2048)
    parser.add_argument("--sample-size", type=int, default=500_000)
    parser.add_argument("--niter", type=int, default=25)
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        shards = pair_shards(args.embeddings, args.manifests)
        if args.centroids_dir:
            centroids, source_metadata = load_existing_quantization(args.centroids_dir)
            metadata = _quantization_metadata(
                centroids,
                int(source_metadata.get("requested_k", centroids.shape[0])),
                shards,
                fitted=False,
            )
            metadata.update(
                {
                    "source_quantization": str(args.centroids_dir),
                    "source_fit_seed": source_metadata.get("seed"),
                    "source_fit_niter": source_metadata.get("niter"),
                }
            )
        else:
            training, total_nodes = sample_training_nodes(shards, args.sample_size, args.seed)
            centroids, _effective_k = fit_centroids(training, args.k, args.niter, args.seed)
            metadata = _quantization_metadata(
                centroids,
                args.k,
                shards,
                sample_size=min(args.sample_size, total_nodes),
                niter=args.niter,
                seed=args.seed,
                total_nodes=total_nodes,
                fitted=True,
            )
        result = write_quantization(shards, args.outdir, centroids, metadata)
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"outdir": str(args.outdir), **result}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
