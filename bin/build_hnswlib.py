#!/usr/bin/env python3
"""Build a custom-distance HNSWLIB index over node-PQ code windows."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from record_pack import pack_records, write_records
from node_quantization import load_json, node_code_bytes


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


def compile_driver(bundle: Path, output: Path) -> None:
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++11",
            "-fopenmp",
            "-I",
            str(bundle),
            str(bundle / "pq_hnswlib.cpp"),
            "-o",
            str(output),
        ],
        check=True,
    )


def load_packed_windows(windows_dir: Path) -> tuple[np.ndarray, list[tuple[str, int, int]], dict[str, Any]]:
    shards = pair_window_shards(windows_dir)
    mapping: list[tuple[str, int, int]] = []
    blocks: list[np.ndarray] = []
    reference: dict[str, Any] | None = None
    for window_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        if reference is None:
            reference = manifest
        window_size = int(manifest["window_size"])
        stride = int(manifest["stride"])
        with np.load(window_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                values = np.asarray(arrays[identifier])
                if values.ndim != 2:
                    raise ValueError(f"{identifier} has invalid packed-window shape {values.shape}")
                blocks.append(np.ascontiguousarray(values))
                for offset in range(values.shape[0]):
                    start = offset * stride
                    mapping.append((identifier, start, start + window_size))
    if not blocks or reference is None:
        raise ValueError("no packed PQ windows found")
    return np.concatenate(blocks, axis=0), mapping, reference


def build_database(
    windows_dir: Path,
    quantization: Path,
    hnsw_bundle: Path,
    outdir: Path,
    *,
    m: int,
    ef_construction: int,
    ef_search: int,
    random_seed: int,
    num_threads: int,
    embeddings: list[Path] | None = None,
    graph_metadata: list[Path] | None = None,
) -> dict[str, Any]:
    quantizer = load_json(quantization / "quantizer.json")
    if quantizer.get("mode") not in {"pq", "opq"}:
        raise ValueError("HNSWLIB custom distance requires --quantize pq or opq")
    codes, mapping, window_meta = load_packed_windows(windows_dir)
    pq_m = int(quantizer["pq_m"])
    nbits = int(quantizer["pq_nbits"])
    window_size = int(window_meta["window_size"])
    expected_width = window_size * node_code_bytes(pq_m, nbits)
    if codes.shape[1] != expected_width:
        raise ValueError(f"packed window width {codes.shape[1]} != {expected_width}")
    outdir.mkdir(parents=True, exist_ok=True)
    quant_out = outdir / "quantization"
    if quant_out.exists():
        shutil.rmtree(quant_out)
    shutil.copytree(quantization, quant_out)
    codes_path = outdir / "window_codes.bin"
    codes.tofile(codes_path)
    similarity = np.load(quantization / "sdc_lut.npy")
    similarity_path = outdir / "sdc_lut.bin"
    np.ascontiguousarray(similarity, dtype=np.float32).tofile(similarity_path)
    executable = outdir / "pq_hnswlib"
    compile_driver(hnsw_bundle, executable)
    index_path = outdir / "index.bin"
    command = [
        str(executable),
        "build",
        "--codes",
        str(codes_path),
        "--similarity",
        str(similarity_path),
        "--index",
        str(index_path),
        "--count",
        str(codes.shape[0]),
        "--window-size",
        str(window_size),
        "--pq-m",
        str(pq_m),
        "--nbits",
        str(nbits),
        "--m",
        str(m),
        "--ef-construction",
        str(ef_construction),
        "--ef-search",
        str(ef_search),
        "--random-seed",
        str(random_seed),
        "--threads",
        str(num_threads),
    ]
    subprocess.run(command, check=True)
    with (outdir / "windows.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        for faiss_id, (transcript_id, start, end) in enumerate(mapping):
            writer.writerow([faiss_id, transcript_id, start, end])
    if embeddings and graph_metadata:
        packed_embeddings, packed_records = pack_records(embeddings, graph_metadata)
        write_records(outdir, packed_embeddings, packed_records)
    meta = {
        "backend": "hnswlib",
        "index_type": "HNSWLIB_PQ",
        "quantize": quantizer.get("mode"),
        "candidate_representation": "node_pq_codes",
        "distance": "sdc_build_adc_search",
        "n_windows": int(codes.shape[0]),
        "window_size": window_size,
        "window_stride": int(window_meta.get("stride", 1)),
        "pq_m": pq_m,
        "pq_nbits": nbits,
        "embedding_dim": int(quantizer.get("embedding_dim", 128)),
        "hnsw_m": m,
        "hnsw_ef_construction": ef_construction,
        "hnsw_ef_search": ef_search,
        "random_seed": random_seed,
    }
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return meta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--quantization", type=Path, required=True)
    parser.add_argument("--hnsw-bundle", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--m", type=int, default=32)
    parser.add_argument("--ef-construction", type=int, default=200)
    parser.add_argument("--ef-search", type=int, default=200)
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--num-threads", type=int, default=0)
    parser.add_argument("--embeddings", nargs="*", type=Path, default=None)
    parser.add_argument("--graph-metadata", nargs="*", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_database(
        args.windows_dir,
        args.quantization,
        args.hnsw_bundle,
        args.outdir,
        m=args.m,
        ef_construction=args.ef_construction,
        ef_search=args.ef_search,
        random_seed=args.random_seed,
        num_threads=args.num_threads,
        embeddings=args.embeddings,
        graph_metadata=args.graph_metadata,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
