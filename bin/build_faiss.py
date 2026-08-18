#!/usr/bin/env python3
"""Build a reusable FAISS window database from embedding shards."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from faiss_index import IndexOptions, build_populated_index, meta_from_details
from cuvs_index import build_populated_index as build_cuvs_index, meta_from_details as cuvs_meta_from_details, serialize_index as serialize_cuvs_index
from ngt_index import build_populated_index as build_ngt_index, meta_from_details as ngt_meta_from_details, serialize_index as serialize_ngt_index
from scann_index import build_populated_searcher, is_scann_type, serialize_index

SLICE_ID_RE = re.compile(r"^(?P<base>.+):(?P<start>\d+)-(?P<end>\d+)$")


COMPAT_KEYS = ("window_size", "stride", "window_dim", "checkpoint_sha256")


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def parse_slice_id(identifier: str) -> tuple[int, int] | None:
    match = SLICE_ID_RE.fullmatch(identifier)
    if not match:
        return None
    return int(match["start"]), int(match["end"])


def pair_table(structure: str) -> list[int]:
    partners = [-1] * len(structure)
    stack: list[int] = []
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                continue
            opening = stack.pop()
            partners[opening] = index
            partners[index] = opening
    return partners


def close_window_structure(structure: str, start: int, end: int) -> str:
    """Keep only pairs whose both ends lie inside ``[start, end)``."""
    partners = pair_table(structure)
    out: list[str] = []
    for index in range(start, end):
        partner = partners[index]
        if partner < start or partner >= end:
            out.append(".")
        else:
            out.append("(" if partner > index else ")")
    return "".join(out)


def subject_sequence(identifier: str, sequence: str, structure: str) -> tuple[str, str]:
    """Return the independent subject (core window, or the full molecule)."""
    if len(sequence) != len(structure):
        raise ValueError(f"{identifier} sequence/structure length mismatch")
    span = parse_slice_id(identifier)
    if span is None:
        return sequence, structure
    start, end = span
    if end - start == len(sequence):
        return sequence, close_window_structure(structure, 0, len(structure))
    if not (0 <= start < end <= len(sequence)):
        raise ValueError(
            f"{identifier} slice [{start}, {end}) is outside a {len(sequence)} nt sequence"
        )
    return sequence[start:end], close_window_structure(structure, start, end)


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


def build_index(
    shards: list[tuple[Path, Path]],
    options: IndexOptions | None = None,
    backend: str = "faiss",
    ngt_index_type: str = "ngt",
    ngt_options: dict[str, Any] | None = None,
    cuvs_index_type: str = "cagra",
    cuvs_options: dict[str, Any] | None = None,
) -> tuple[Any, list[tuple[int, str, int, int]], dict]:
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
    chosen = options or IndexOptions()
    if backend == "ngt":
        index, details = build_ngt_index(xb, ngt_index_type, ngt_options)
    elif backend == "cuvs":
        index, details = build_cuvs_index(xb, cuvs_index_type, **(cuvs_options or {}))
    elif is_scann_type(chosen.index_type):
        index, details = build_populated_searcher(xb, chosen)
    else:
        index, details = build_populated_index(xb, chosen)

    assert reference_manifest is not None
    meta = {
        "window_size": int(reference_manifest["window_size"]),
        "window_stride": int(reference_manifest["stride"]),
        "embedding_dim": int(reference_manifest["embedding_dim"]),
        "window_dim": int(reference_manifest["window_dim"]),
        "l2_normalized": True,
        "ginfinity_version": reference_manifest.get("ginfinity_version"),
        "model_version": reference_manifest.get("model_version"),
        "checkpoint_sha256": reference_manifest.get("checkpoint_sha256"),
        "n_records": n_records,
        "n_windows": int(xb.shape[0]),
        "n_skipped_short": n_skipped,
    }
    if backend == "ngt":
        meta.update(ngt_meta_from_details(details))
    elif backend == "cuvs":
        meta.update(cuvs_meta_from_details(details))
    else:
        meta.update(meta_from_details(details))
    return index, mapping, meta


def pack_records(embedding_paths: list[Path], metadata_paths: list[Path]) -> tuple[dict[str, np.ndarray], list[tuple[str, str, str]]]:
    embeddings: dict[str, np.ndarray] = {}
    for path in embedding_paths:
        with np.load(path) as archive:
            for key in archive.files:
                if key in embeddings:
                    raise ValueError(f"duplicate embedding id {key!r} in {path}")
                embeddings[key] = np.asarray(archive[key])

    records: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for path in metadata_paths:
        payload = load_json(path)
        identifiers = payload.get("identifiers")
        sequences = payload.get("sequences")
        structures = payload.get("structures")
        if not (identifiers and sequences and structures):
            raise ValueError(f"{path} is not a GINFINITY graph sidecar")
        if not (len(identifiers) == len(sequences) == len(structures)):
            raise ValueError(f"{path} identifier/sequence/structure length mismatch")
        for identifier, sequence, structure in zip(identifiers, sequences, structures):
            if identifier in seen:
                raise ValueError(f"duplicate record id {identifier!r} in {path}")
            if identifier not in embeddings:
                raise ValueError(f"{identifier} is in {path} but missing from residue embeddings")
            sequence, structure = subject_sequence(identifier, sequence, structure)
            if len(sequence) != embeddings[identifier].shape[0]:
                raise ValueError(
                    f"{identifier} sequence length {len(sequence)} does not match "
                    f"embedding rows {embeddings[identifier].shape[0]}"
                )
            seen.add(identifier)
            records.append((identifier, sequence, structure))

    missing = sorted(set(embeddings) - seen)
    if missing:
        raise ValueError(f"embeddings without graph metadata: {missing[:8]}")
    return embeddings, records


def write_database(
    outdir: Path,
    index: Any,
    mapping: list[tuple[int, str, int, int]],
    meta: dict,
    packed: tuple[dict[str, np.ndarray], list[tuple[str, str, str]]] | None = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if str(meta.get("backend") or "").lower() == "ngt":
        serialize_ngt_index(index, outdir / "ngt")
    elif str(meta.get("backend") or "").lower() == "cuvs":
        serialize_cuvs_index(index, outdir / "cuvs")
    elif is_scann_type(str(meta.get("index_type") or "")):
        serialize_index(index, outdir / "scann")
    else:
        import faiss

        faiss.write_index(index, str(outdir / "index.faiss"))
    with (outdir / "windows.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["faiss_id", "transcript_id", "start", "end"])
        writer.writerows(mapping)
    if packed is not None:
        embeddings, records = packed
        np.savez_compressed(outdir / "embeddings.npz", **embeddings)
        with (outdir / "records.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["transcript_id", "sequence", "secondary_structure"])
            writer.writerows(records)
        meta["has_residue_embeddings"] = True
        meta["n_packed_records"] = len(records)
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, nargs="+", required=True)
    parser.add_argument("--manifests", type=Path, nargs="+", required=True)
    parser.add_argument("--embeddings", type=Path, nargs="*")
    parser.add_argument("--graph-metadata", type=Path, nargs="*")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--backend", choices=("faiss", "scann", "ngt", "cuvs"), default="faiss")
    parser.add_argument("--index-type", default="flatip")
    parser.add_argument("--ngt-index-type", default="ngt")
    parser.add_argument("--ngt-edge-size-for-creation", type=int, default=10)
    parser.add_argument("--ngt-edge-size-for-search", type=int, default=40)
    parser.add_argument("--ngt-num-threads", type=int, default=8)
    parser.add_argument("--ngt-max-no-of-edges", type=int)
    parser.add_argument("--ngt-num-of-search-objects", type=int, default=20)
    parser.add_argument("--ngt-search-range-coefficient", type=float)
    parser.add_argument("--ngt-blob-search-range-coefficient", type=float)
    parser.add_argument("--ngt-search-radius", type=float)
    parser.add_argument("--ngt-result-expansion", type=float)
    parser.add_argument("--ngt-exploration-size", type=int, default=256)
    parser.add_argument("--ngt-exact-result-expansion", type=float, default=0.0)
    parser.add_argument("--ngt-num-of-probes", type=int, default=1)
    parser.add_argument("--ngt-qg-subvector-dimensions", type=int)
    parser.add_argument("--ngt-qbg-subvectors", type=int)
    parser.add_argument("--ngt-qbg-cluster-data-type", default="PQ4")
    parser.add_argument("--cuvs-index-type", default="cagra")
    parser.add_argument("--cuvs-n-lists", type=int)
    parser.add_argument("--cuvs-n-probes", type=int)
    parser.add_argument("--cuvs-pq-bits", type=int, default=8)
    parser.add_argument("--cuvs-pq-dim", type=int, default=0)
    parser.add_argument("--cuvs-intermediate-graph-degree", type=int, default=128)
    parser.add_argument("--cuvs-graph-degree", type=int, default=64)
    parser.add_argument("--cuvs-build-algo", default="nn_descent")
    parser.add_argument("--cuvs-itopk-size", type=int, default=64)
    parser.add_argument("--nlist", type=int)
    parser.add_argument("--nprobe", type=int)
    parser.add_argument("--pq-m", type=int, default=16)
    parser.add_argument("--pq-nbits", type=int, default=8)
    parser.add_argument("--pq-m-refine", type=int, default=4)
    parser.add_argument("--hnsw-m", type=int, default=32)
    parser.add_argument("--hnsw-ef-construction", type=int, default=40)
    parser.add_argument("--hnsw-ef-search", type=int, default=16)
    parser.add_argument("--lsh-nbits", type=int)
    parser.add_argument("--sq-type", default="8bit")
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-device", type=int, default=0)
    parser.add_argument("--scann-reorder", type=int, default=100)
    parser.add_argument("--scann-ah-dim", type=int, default=2)
    parser.add_argument("--scann-anisotropic", type=float, default=0.2)
    parser.add_argument("--scann-soar", action="store_true")
    parser.add_argument("--scann-leaves", type=int)
    parser.add_argument("--scann-leaves-to-search", type=int)
    parser.add_argument("--num-neighbors", type=int, default=100)
    return parser.parse_args(argv)


def options_from_args(args: argparse.Namespace) -> IndexOptions:
    return IndexOptions(
        index_type=args.index_type,
        nlist=args.nlist,
        nprobe=args.nprobe,
        pq_m=args.pq_m,
        pq_nbits=args.pq_nbits,
        pq_m_refine=args.pq_m_refine,
        hnsw_m=args.hnsw_m,
        hnsw_ef_construction=args.hnsw_ef_construction,
        hnsw_ef_search=args.hnsw_ef_search,
        lsh_nbits=args.lsh_nbits,
        sq_type=args.sq_type,
        gpu=args.gpu,
        gpu_device=args.gpu_device,
        scann_reorder=args.scann_reorder,
        scann_ah_dim=args.scann_ah_dim,
        scann_anisotropic=args.scann_anisotropic,
        scann_soar=args.scann_soar,
        scann_num_neighbors=args.num_neighbors,
        scann_leaves=args.scann_leaves,
        scann_leaves_to_search=args.scann_leaves_to_search,
    )


def cuvs_options_from_args(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "n_lists": args.cuvs_n_lists,
        "n_probes": args.cuvs_n_probes,
        "pq_bits": args.cuvs_pq_bits,
        "pq_dim": args.cuvs_pq_dim,
        "intermediate_graph_degree": args.cuvs_intermediate_graph_degree,
        "graph_degree": args.cuvs_graph_degree,
        "build_algo": args.cuvs_build_algo,
        "itopk_size": args.cuvs_itopk_size,
    }


def ngt_options_from_args(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "edge_size_for_creation": args.ngt_edge_size_for_creation,
        "edge_size_for_search": args.ngt_edge_size_for_search,
        "num_threads": args.ngt_num_threads,
        "max_no_of_edges": args.ngt_max_no_of_edges,
        "num_of_search_objects": args.ngt_num_of_search_objects,
        "search_range_coefficient": args.ngt_search_range_coefficient,
        "blob_search_range_coefficient": args.ngt_blob_search_range_coefficient,
        "search_radius": args.ngt_search_radius,
        "result_expansion": args.ngt_result_expansion,
        "exploration_size": args.ngt_exploration_size,
        "exact_result_expansion": args.ngt_exact_result_expansion,
        "num_of_probes": args.ngt_num_of_probes,
        "qg_subvector_dimensions": args.ngt_qg_subvector_dimensions,
        "qbg_subvectors": args.ngt_qbg_subvectors,
        "qbg_cluster_data_type": args.ngt_qbg_cluster_data_type,
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        shards = pair_shards(args.windows, args.manifests)
        index, mapping, meta = build_index(
            shards,
            options_from_args(args),
            backend=args.backend,
            ngt_index_type=args.ngt_index_type,
            ngt_options=ngt_options_from_args(args),
            cuvs_index_type=args.cuvs_index_type,
            cuvs_options=cuvs_options_from_args(args),
        )
        packed = None
        if args.embeddings or args.graph_metadata:
            if not args.embeddings or not args.graph_metadata:
                raise ValueError("--embeddings and --graph-metadata must be passed together")
            packed = pack_records(args.embeddings, args.graph_metadata)
    except (OSError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    write_database(args.outdir, index, mapping, meta, packed)
    keys = ("n_windows", "n_records", "window_dim", "index_type", "metric")
    print(json.dumps({"outdir": str(args.outdir), **{k: meta[k] for k in keys if k in meta}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
