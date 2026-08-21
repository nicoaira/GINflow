#!/usr/bin/env python3
"""Build a node-PQ CAGRA graph (GPU) over pre-quantized windows."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    import pq_cagra_adc as pq
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "pq_cagra_adc is not installed. Use -profile conda or -profile docker "
        "so BUILD_PQ_CAGRA_INDEX gets nicolas.aira::pq-cagra-adc."
    ) from exc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from node_quantization import load_json, node_code_bytes, unpack_pq_codes
from record_pack import pack_records, write_records


def pair_window_shards(windows_dir: Path) -> list[tuple[Path, Path]]:
    window_paths = sorted(windows_dir.glob("*.windows.npz"))
    manifest_paths = sorted(windows_dir.glob("*.windows.manifest.json"))
    window_map = {path.name[: -len(".windows.npz")]: path for path in window_paths}
    manifest_map = {path.name[: -len(".windows.manifest.json")]: path for path in manifest_paths}
    if set(window_map) != set(manifest_map) or not window_map:
        raise ValueError(f"quantized window shards are missing in {windows_dir}")
    return [(window_map[key], manifest_map[key]) for key in sorted(window_map)]


def load_unpacked_windows(windows_dir: Path, pq_m: int, nbits: int) -> tuple[np.ndarray, list[tuple[str, int, int]], dict[str, Any]]:
    blocks: list[np.ndarray] = []
    mapping: list[tuple[str, int, int]] = []
    reference: dict[str, Any] | None = None
    width = None
    for window_path, manifest_path in pair_window_shards(windows_dir):
        manifest = load_json(manifest_path)
        if reference is None:
            reference = manifest
        window_size = int(manifest["window_size"])
        stride = int(manifest["stride"])
        packed_width = window_size * node_code_bytes(pq_m, nbits)
        with np.load(window_path) as arrays:
            for record in manifest.get("records", []):
                identifier = str(record["identifier"])
                packed = np.asarray(arrays[identifier], dtype=np.uint8)
                if packed.ndim != 2 or packed.shape[1] != packed_width:
                    raise ValueError(f"{identifier} packed width {packed.shape} != {packed_width}")
                unpacked = unpack_pq_codes(
                    packed.reshape(packed.shape[0] * window_size, node_code_bytes(pq_m, nbits)),
                    pq_m,
                    nbits,
                ).reshape(packed.shape[0], window_size, pq_m)
                blocks.append(np.ascontiguousarray(unpacked))
                for offset in range(unpacked.shape[0]):
                    start = offset * stride
                    mapping.append((identifier, start, start + window_size))
        width = packed_width
    if not blocks or reference is None:
        raise ValueError("no PQ windows to index")
    _ = width
    return np.concatenate(blocks, axis=0), mapping, reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows-dir", type=Path, required=True)
    parser.add_argument("--quantization", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--graph-degree", type=int, default=32)
    parser.add_argument("--intermediate-graph-degree", type=int, default=64)
    parser.add_argument("--nndescent-iterations", type=int, default=6)
    parser.add_argument("--itopk-size", type=int, default=1024)
    parser.add_argument("--embeddings", nargs="*", type=Path, default=None)
    parser.add_argument("--graph-metadata", nargs="*", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    quantizer = load_json(args.quantization / "quantizer.json")
    if quantizer.get("mode") not in {"pq", "opq"}:
        raise SystemExit("PQ-CAGRA requires --quantize pq or opq")
    pq_m = int(quantizer["pq_m"])
    nbits = int(quantizer["pq_nbits"])
    codebook = np.load(args.quantization / "codebook.npy")
    rotation_path = args.quantization / "rotation.npy"
    codes, mapping, window_meta = load_unpacked_windows(args.windows_dir, pq_m, nbits)
    index = pq.build_index(
        codes,
        codebook,
        graph_degree=args.graph_degree,
        intermediate_graph_degree=args.intermediate_graph_degree,
        nndescent_iterations=args.nndescent_iterations,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    pq.save_index(index, args.outdir / "cagra.index")
    shutil.copytree(args.quantization, args.outdir / "quantization", dirs_exist_ok=True)
    with (args.outdir / "windows.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        for faiss_id, (transcript_id, start, end) in enumerate(mapping):
            writer.writerow([faiss_id, transcript_id, start, end])
    if args.embeddings and args.graph_metadata:
        packed_embeddings, packed_records = pack_records(args.embeddings, args.graph_metadata)
        write_records(args.outdir, packed_embeddings, packed_records)
    meta = {
        "backend": "cagra",
        "index_type": "PQ_CAGRA",
        "quantize": quantizer.get("mode"),
        "n_windows": int(codes.shape[0]),
        "window_size": int(window_meta["window_size"]),
        "window_stride": int(window_meta.get("stride", 1)),
        "pq_m": pq_m,
        "pq_nbits": nbits,
        "embedding_dim": int(quantizer.get("embedding_dim", 128)),
        "graph_degree": args.graph_degree,
        "intermediate_graph_degree": args.intermediate_graph_degree,
        "nndescent_iterations": args.nndescent_iterations,
        "itopk_size": args.itopk_size,
        "has_rotation": rotation_path.exists(),
        "gpu_build": True,
        "cpu_search": True,
    }
    (args.outdir / "meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"n_windows": meta["n_windows"], "index_type": "PQ_CAGRA"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
