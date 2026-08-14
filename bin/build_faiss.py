#!/usr/bin/env python3
"""Build a FAISS IndexFlatIP database from window embedding shards."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import faiss
import numpy as np


COMPAT_KEYS = ("window_size", "stride", "window_dim", "checkpoint_sha256")


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def shard_prefix(path: Path, suffix: str) -> str:
    name = path.name
    if not name.endswith(suffix):
        raise ValueError(f"{path} does not end with {suffix}")
    return name[: -len(suffix)]


def pair_shards(windows: list[Path], manifests: list[Path]) -> list[tuple[Path, Path]]:
    win_map = {shard_prefix(path, ".windows.npz"): path for path in windows}
    man_map = {shard_prefix(path, ".windows.manifest.json"): path for path in manifests}
    extra_windows = sorted(set(win_map) - set(man_map))
    extra_manifests = sorted(set(man_map) - set(win_map))
    if extra_windows or extra_manifests:
        raise ValueError(
            "window / manifest shard names do not match: "
            f"only-windows={extra_windows} only-manifests={extra_manifests}"
        )
    return [(win_map[key], man_map[key]) for key in sorted(win_map)]


def compat_tuple(manifest: dict) -> tuple:
    return tuple(manifest.get(key) for key in COMPAT_KEYS)


def load_shard(npz_path: Path, manifest: dict) -> tuple[np.ndarray, list[tuple[str, int, int]]]:
    arrays = np.load(npz_path)
    window_size = int(manifest["window_size"])
    vectors = []
    rows = []
    for record in manifest.get("records", []):
        identifier = record["identifier"]
        if identifier not in arrays.files:
            raise KeyError(f"{identifier} is in {npz_path.name} manifest but missing from the NPZ")
        windows = np.asarray(arrays[identifier], dtype=np.float32)
        if windows.ndim != 2:
            raise ValueError(f"{identifier} windows must be 2-D, got {windows.shape}")
        for offset, vector in enumerate(windows):
            start = offset * int(manifest["stride"])
            rows.append((identifier, start, start + window_size))
            vectors.append(vector)
    if not vectors:
        return np.empty((0, int(manifest["window_dim"])), dtype=np.float32), []
    stacked = np.ascontiguousarray(np.stack(vectors, axis=0), dtype=np.float32)
    return stacked, rows


def build_index(shards: list[tuple[Path, Path]]) -> tuple[faiss.Index, list[tuple[int, str, int, int]], dict]:
    if not shards:
        raise ValueError("no window shards were provided")

    all_vectors = []
    mapping = []
    expected = None
    n_records = 0
    n_skipped = 0
    reference_manifest = None

    for npz_path, manifest_path in shards:
        manifest = load_json(manifest_path)
        if expected is None:
            expected = compat_tuple(manifest)
            reference_manifest = manifest
        elif compat_tuple(manifest) != expected:
            raise ValueError(
                f"{manifest_path} is incompatible with the first shard "
                f"({dict(zip(COMPAT_KEYS, expected))} vs "
                f"{dict(zip(COMPAT_KEYS, compat_tuple(manifest)))})"
            )
        vectors, rows = load_shard(npz_path, manifest)
        n_records += len(manifest.get("records", [])) + len(manifest.get("skipped_short", []))
        n_skipped += len(manifest.get("skipped_short", []))
        base = sum(item.shape[0] for item in all_vectors)
        for local_id, (identifier, start, end) in enumerate(rows):
            mapping.append((base + local_id, identifier, start, end))
        if vectors.shape[0]:
            all_vectors.append(vectors)

    if not all_vectors:
        raise ValueError("no windows to index (every sequence was shorter than --window-size)")

    xb = np.ascontiguousarray(np.concatenate(all_vectors, axis=0), dtype=np.float32)
    index = faiss.IndexFlatIP(xb.shape[1])
    index.add(xb)

    assert reference_manifest is not None
    meta = {
        "window_size": int(reference_manifest["window_size"]),
        "window_stride": int(reference_manifest["stride"]),
        "embedding_dim": int(reference_manifest["embedding_dim"]),
        "window_dim": int(reference_manifest["window_dim"]),
        "metric": "inner_product",
        "l2_normalized": True,
        "ginfinity_version": reference_manifest.get("ginfinity_version"),
        "model_version": reference_manifest.get("model_version"),
        "checkpoint_sha256": reference_manifest.get("checkpoint_sha256"),
        "n_records": n_records,
        "n_windows": int(xb.shape[0]),
        "n_skipped_short": n_skipped,
    }
    return index, mapping, meta


def write_database(outdir: Path, index: faiss.Index, mapping: list[tuple[int, str, int, int]], meta: dict) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    faiss.write_index(index, str(outdir / "index.faiss"))
    with (outdir / "windows.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        writer.writerows(mapping)
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, nargs="+", required=True)
    parser.add_argument("--manifests", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        shards = pair_shards(args.windows, args.manifests)
        index, mapping, meta = build_index(shards)
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    write_database(args.outdir, index, mapping, meta)
    print(json.dumps({"outdir": str(args.outdir), **{k: meta[k] for k in ("n_windows", "n_records", "window_dim")}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
