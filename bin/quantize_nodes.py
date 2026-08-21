#!/usr/bin/env python3
"""Fit or apply node-level SQ / PQ / OPQ and persist codebook + SDC artifacts."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from node_quantization import (
    collect_nodes,
    encode_pq,
    load_json,
    pair_shards,
    rotate_nodes,
    sdc_lookup_table,
    train_node_opq,
    train_node_pq,
    train_node_sq,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", nargs="+", type=Path, required=True)
    parser.add_argument("--manifests", nargs="+", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mode", choices=("sq", "pq", "opq"), required=True)
    parser.add_argument("--quantizer-dir", type=Path, default=None, help="Apply a published quantizer instead of fitting.")
    parser.add_argument("--pq-m", type=int, default=16)
    parser.add_argument("--pq-nbits", type=int, default=4)
    parser.add_argument("--sample-size", type=int, default=250000)
    parser.add_argument("--niter", type=int, default=12)
    parser.add_argument("--opq-iters", type=int, default=10)
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def save_node_shards(
    shards: list[tuple[Path, Path]],
    codes: np.ndarray,
    outdir: Path,
    code_width: int,
) -> list[dict[str, Any]]:
    nodes_dir = outdir / "nodes"
    nodes_dir.mkdir(parents=True, exist_ok=True)
    offset = 0
    summaries: list[dict[str, Any]] = []
    for embedding_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        arrays: dict[str, np.ndarray] = {}
        records: list[dict[str, Any]] = []
        for record in manifest.get("records", []):
            identifier = str(record["identifier"])
            length = int(record.get("length", record.get("core_length", 0)))
            if length < 1:
                with np.load(embedding_path) as source:
                    length = int(np.asarray(source[identifier]).shape[0])
            node_codes = codes[offset : offset + length]
            if node_codes.shape[0] != length:
                raise ValueError(f"{identifier}: expected {length} node codes, got {node_codes.shape[0]}")
            arrays[identifier] = np.ascontiguousarray(node_codes)
            records.append({"identifier": identifier, "length": length})
            offset += length
        prefix = embedding_path.name[: -len(".npz")] if embedding_path.name.endswith(".npz") else embedding_path.stem
        np.savez_compressed(nodes_dir / f"{prefix}.quantized.npz", **arrays)
        write_json(
            nodes_dir / f"{prefix}.quantized.manifest.json",
            {
                "records": records,
                "code_width": code_width,
                "source_embeddings": embedding_path.name,
            },
        )
        summaries.extend(records)
    if offset != codes.shape[0]:
        raise ValueError(f"wrote {offset} node codes but encoded {codes.shape[0]}")
    return summaries


def fit_or_apply(args: argparse.Namespace) -> dict[str, Any]:
    shards = pair_shards(args.embeddings, args.manifests)
    nodes, _records = collect_nodes(shards)
    dimension = int(nodes.shape[1])
    rotation = None
    codebook = None
    scale = None
    zero = None
    if args.quantizer_dir:
        quantizer = load_json(args.quantizer_dir / "quantizer.json")
        mode = str(quantizer["mode"])
        if mode != args.mode:
            raise ValueError(f"published quantizer mode {mode} does not match --mode {args.mode}")
        if mode == "sq":
            scale = np.load(args.quantizer_dir / "sq_scale.npy")
            zero = np.load(args.quantizer_dir / "sq_zero.npy")
            codes = np.clip(np.rint((nodes - zero) / np.maximum(scale, 1e-12)), 0, 255).astype(np.uint8)
            code_width = dimension
        else:
            codebook = np.load(args.quantizer_dir / "codebook.npy")
            if mode == "opq":
                rotation = np.load(args.quantizer_dir / "rotation.npy")
                encoded = rotate_nodes(nodes, rotation)
            else:
                encoded = nodes
            codes = encode_pq(encoded, codebook)
            code_width = int(quantizer["pq_m"])
        meta = dict(quantizer)
        meta["n_nodes"] = int(nodes.shape[0])
        meta["applied"] = True
        return {
            "meta": meta,
            "codes": codes,
            "codebook": codebook,
            "rotation": rotation,
            "scale": scale,
            "zero": zero,
            "code_width": code_width,
        }

    if args.mode == "sq":
        scale, zero, codes = train_node_sq(nodes)
        meta = {
            "mode": "sq",
            "embedding_dim": dimension,
            "n_nodes": int(nodes.shape[0]),
            "sq_type": "8bit",
        }
        code_width = dimension
    elif args.mode == "pq":
        codebook, codes = train_node_pq(
            nodes, args.pq_m, args.pq_nbits, args.sample_size, args.niter, args.seed
        )
        meta = {
            "mode": "pq",
            "embedding_dim": dimension,
            "n_nodes": int(nodes.shape[0]),
            "pq_m": args.pq_m,
            "pq_nbits": args.pq_nbits,
            "sample_size": args.sample_size,
            "niter": args.niter,
            "seed": args.seed,
        }
        code_width = args.pq_m
    else:
        rotation, codebook, codes = train_node_opq(
            nodes,
            args.pq_m,
            args.pq_nbits,
            args.sample_size,
            args.niter,
            args.opq_iters,
            args.seed,
        )
        meta = {
            "mode": "opq",
            "embedding_dim": dimension,
            "n_nodes": int(nodes.shape[0]),
            "pq_m": args.pq_m,
            "pq_nbits": args.pq_nbits,
            "sample_size": args.sample_size,
            "niter": args.niter,
            "opq_iters": args.opq_iters,
            "seed": args.seed,
        }
        code_width = args.pq_m
    return {
        "meta": meta,
        "codes": codes,
        "codebook": codebook,
        "rotation": rotation,
        "scale": scale,
        "zero": zero,
        "code_width": code_width,
    }


def main() -> int:
    args = parse_args()
    if args.outdir.exists():
        shutil.rmtree(args.outdir)
    args.outdir.mkdir(parents=True)
    shards = pair_shards(args.embeddings, args.manifests)
    fitted = fit_or_apply(args)
    meta = fitted["meta"]
    codes = fitted["codes"]
    if fitted["codebook"] is not None:
        np.save(args.outdir / "codebook.npy", fitted["codebook"])
        np.save(args.outdir / "sdc_lut.npy", sdc_lookup_table(fitted["codebook"]))
    if fitted["rotation"] is not None:
        np.save(args.outdir / "rotation.npy", fitted["rotation"])
    if fitted["scale"] is not None and fitted["zero"] is not None:
        np.save(args.outdir / "sq_scale.npy", fitted["scale"])
        np.save(args.outdir / "sq_zero.npy", fitted["zero"])
    save_node_shards(shards, codes, args.outdir, fitted["code_width"])
    write_json(args.outdir / "quantizer.json", meta)
    print(
        f"wrote {args.mode} quantizer: nodes={meta['n_nodes']} dim={meta['embedding_dim']}"
        + (f" pq_m={meta.get('pq_m')} nbits={meta.get('pq_nbits')}" if args.mode != "sq" else "")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
