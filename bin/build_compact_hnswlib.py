#!/usr/bin/env python3
"""Build a compact custom-distance hnswlib index over quantized windows."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compact_hnswlib_io import copy_quantization, pack_code_windows, write_mapping, write_similarity_binary
from hnswlib_index import load_quantization, quantization_dir
from record_pack import pack_records, write_records


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def build_database(
    windows_dir: Path,
    quantization: Path,
    outdir: Path,
    executable: Path,
    *,
    m: int,
    ef_construction: int,
    ef_search: int,
    random_seed: int,
    num_threads: int,
    embeddings: list[Path] | None = None,
    graph_metadata: list[Path] | None = None,
) -> dict[str, Any]:
    quant_dir = quantization_dir(quantization)
    centroids, similarity, quantization_metadata = load_quantization(quant_dir)
    if centroids.shape[0] > np.iinfo(np.uint16).max + 1:
        raise ValueError("compact hnswlib supports at most 65536 centroids")
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="compact-hnsw-", dir=outdir) as work:
        workdir = Path(work)
        codes_path = workdir / "codes.bin"
        n_windows, window_size, mapping, manifests = pack_code_windows(windows_dir, codes_path)
        similarity_path = workdir / "similarity.bin"
        write_similarity_binary(quant_dir, similarity_path)
        index_path = outdir / "index.bin"
        executable_path = executable if executable.is_absolute() else Path.cwd() / executable
        command = [
            str(executable_path),
            "build",
            "--codes", str(codes_path),
            "--similarity", str(similarity_path),
            "--count", str(n_windows),
            "--window-size", str(window_size),
            "--n-centroids", str(centroids.shape[0]),
            "--index", str(index_path),
            "--m", str(m),
            "--ef-construction", str(ef_construction),
            "--ef-search", str(ef_search),
            "--random-seed", str(random_seed),
            "--threads", str(num_threads),
        ]
        subprocess.run(command, check=True)

    write_mapping(outdir / "windows.tsv", mapping)
    copy_quantization(quant_dir, outdir / "quantization")
    reference = manifests[0]
    meta: dict[str, Any] = {
        "backend": "hnswlib",
        "index_type": "HNSWLIB_COMPACT",
        "index_class": "hnswlib::HierarchicalNSW<float>",
        "space": "custom_code_sum",
        "metric": "sum_positionwise_centroid_similarity",
        "raw_metric": "sum_positionwise_centroid_similarity",
        "distance_function": "negative_sum_positionwise_centroid_similarity",
        "candidate_representation": "uint16_centroid_codes",
        "quantized_nodes": True,
        "quantization": "quantization",
        "n_centroids": int(centroids.shape[0]),
        "centroid_dim": int(centroids.shape[1]),
        "centroid_dtype": str(centroids.dtype),
        "similarity_dtype": str(similarity.dtype),
        "window_size": int(reference["window_size"]),
        "window_stride": int(reference["stride"]),
        "embedding_dim": int(reference["embedding_dim"]),
        "window_dim": int(reference["window_size"]) * int(reference["embedding_dim"]),
        "index_data_dtype": "uint16",
        "index_data_bytes_per_element": int(reference["window_size"]) * np.dtype(np.uint16).itemsize,
        "l2_normalized": False,
        "score_scale": 1.0,
        "n_records": int(sum(len(manifest.get("records", [])) + len(manifest.get("skipped_short", [])) for manifest in manifests)),
        "n_windows": int(n_windows),
        "n_skipped_short": int(sum(len(manifest.get("skipped_short", [])) for manifest in manifests)),
        "ginfinity_version": reference.get("ginfinity_version"),
        "model_version": reference.get("model_version"),
        "checkpoint_sha256": reference.get("checkpoint_sha256"),
        "hnsw_m": int(m),
        "hnsw_ef_construction": int(ef_construction),
        "hnsw_ef_search": int(ef_search),
        "hnsw_random_seed": int(random_seed),
        "hnsw_num_threads": int(num_threads),
        "hnswlib_source_version": "0.8.0",
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
    parser.add_argument("--executable", type=Path, required=True)
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
            args.executable,
            m=args.m,
            ef_construction=args.ef_construction,
            ef_search=args.ef_search,
            random_seed=args.random_seed,
            num_threads=args.num_threads,
            embeddings=args.embeddings,
            graph_metadata=args.graph_metadata,
        )
    except (OSError, KeyError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"outdir": str(args.outdir), **{key: meta[key] for key in ("n_windows", "n_centroids", "hnsw_m", "hnsw_ef_search", "index_type")}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
