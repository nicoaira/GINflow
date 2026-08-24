#!/usr/bin/env python3
"""Pack original residue embeddings and graph records for searchable databases."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np

PACKED_VECTORS_KEY = "__ginflow_vectors__"
PACKED_OFFSETS_KEY = "__ginflow_offsets__"
PACKED_IDS_KEY = "__ginflow_ids__"
VECTORS_NPY = "embeddings.vectors.npy"
OFFSETS_NPY = "embeddings.offsets.npy"
IDS_TXT = "embeddings.ids.txt"
NPZ_NAME = "embeddings.npz"
PACKED_KEYS = frozenset({PACKED_VECTORS_KEY, PACKED_OFFSETS_KEY, PACKED_IDS_KEY})


SLICE_ID_RE = re.compile(r"^(?P<base>.+):(?P<start>\d+)-(?P<end>\d+)$")


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
        elif char == ")" and stack:
            opening = stack.pop()
            partners[opening] = index
            partners[index] = opening
    return partners


def close_window_structure(structure: str, start: int, end: int) -> str:
    partners = pair_table(structure)
    output: list[str] = []
    for index in range(start, end):
        partner = partners[index]
        if partner < start or partner >= end:
            output.append(".")
        else:
            output.append("(" if partner > index else ")")
    return "".join(output)


def subject_sequence(identifier: str, sequence: str, structure: str) -> tuple[str, str]:
    if len(sequence) != len(structure):
        raise ValueError(f"{identifier} sequence/structure length mismatch")
    span = parse_slice_id(identifier)
    if span is None:
        return sequence, structure
    start, end = span
    if end - start == len(sequence):
        return sequence, close_window_structure(structure, 0, len(structure))
    if not (0 <= start < end <= len(sequence)):
        raise ValueError(f"{identifier} slice [{start}, {end}) is outside a {len(sequence)} nt sequence")
    return sequence[start:end], close_window_structure(structure, start, end)


def pack_records(
    embedding_paths: list[Path], metadata_paths: list[Path]
) -> tuple[dict[str, np.ndarray], list[tuple[str, str, str]]]:
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


def _stack_embeddings(
    embeddings: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    identifiers = list(embeddings)
    if not identifiers:
        raise ValueError("no embeddings to pack")
    parts = [np.ascontiguousarray(embeddings[identifier]) for identifier in identifiers]
    dimension = parts[0].shape[1] if parts[0].ndim == 2 else -1
    for identifier, part in zip(identifiers, parts):
        if part.ndim != 2 or part.shape[1] != dimension:
            raise ValueError(f"{identifier} is not a 2-D embedding of dimension {dimension}")
    vectors = np.concatenate(parts, axis=0)
    offsets = np.zeros(len(identifiers) + 1, dtype=np.int64)
    for index, part in enumerate(parts):
        offsets[index + 1] = offsets[index] + part.shape[0]
    return vectors, offsets, identifiers


def _dict_from_packed(
    vectors: np.ndarray,
    offsets: np.ndarray,
    identifiers: list[str],
    ids: set[str] | None,
    require_all: bool,
) -> dict[str, np.ndarray]:
    offsets = np.asarray(offsets, dtype=np.int64)
    if offsets.ndim != 1 or offsets.shape[0] != len(identifiers) + 1:
        raise ValueError("packed embedding offsets do not match identifiers")
    index = {identifier: position for position, identifier in enumerate(identifiers)}
    if ids is None:
        wanted = identifiers
    else:
        missing = [identifier for identifier in ids if identifier not in index]
        if missing and require_all:
            raise ValueError(
                "missing embedding ids: " + ", ".join(missing[:8])
            )
        wanted = [identifier for identifier in ids if identifier in index]
    arrays: dict[str, np.ndarray] = {}
    for identifier in wanted:
        position = index[identifier]
        arrays[identifier] = vectors[int(offsets[position]):int(offsets[position + 1])]
    return arrays


def _load_npy_pack(
    directory: Path, ids: set[str] | None, require_all: bool
) -> dict[str, np.ndarray]:
    vectors = np.load(directory / VECTORS_NPY, mmap_mode="r")
    offsets = np.load(directory / OFFSETS_NPY)
    identifiers = (directory / IDS_TXT).read_text().splitlines()
    return _dict_from_packed(vectors, offsets, identifiers, ids, require_all)


def _is_packed_npz(archive: np.lib.npyio.NpzFile) -> bool:
    names = set(archive.files)
    return PACKED_VECTORS_KEY in names and PACKED_OFFSETS_KEY in names


def _load_npz(
    path: Path, ids: set[str] | None, require_all: bool
) -> dict[str, np.ndarray]:
    with np.load(path) as archive:
        if _is_packed_npz(archive):
            vectors = np.asarray(archive[PACKED_VECTORS_KEY])
            offsets = np.asarray(archive[PACKED_OFFSETS_KEY])
            if PACKED_IDS_KEY in archive.files:
                identifiers = [str(value) for value in np.asarray(archive[PACKED_IDS_KEY])]
            else:
                raise ValueError(f"{path} packed embeddings are missing identifiers")
            return _dict_from_packed(vectors, offsets, identifiers, ids, require_all)
        names = [name for name in archive.files if name not in PACKED_KEYS]
        if ids is None:
            wanted = names
        else:
            available = set(names)
            missing = [identifier for identifier in ids if identifier not in available]
            if missing and require_all:
                raise ValueError(
                    "missing embedding ids: " + ", ".join(missing[:8])
                )
            wanted = [identifier for identifier in ids if identifier in available]
        arrays: dict[str, np.ndarray] = {}
        for name in wanted:
            arrays[name] = np.asarray(archive[name])
        return arrays


def load_residue_embeddings(
    source: Path,
    ids: set[str] | None = None,
    require_all: bool = True,
) -> dict[str, np.ndarray]:
    """Load residue embeddings from a database directory or an NPZ file.

    New databases store a concatenated memmap plus a compact NPZ. Legacy
    per-id NPZ archives are still read, optionally restricted to ``ids``.
    """
    source = Path(source)
    if source.is_dir():
        npy = source / VECTORS_NPY
        ids_path = source / IDS_TXT
        offsets_path = source / OFFSETS_NPY
        npz = source / NPZ_NAME
        if npy.is_file() and ids_path.is_file() and offsets_path.is_file():
            return _load_npy_pack(source, ids, require_all)
        if npz.is_file():
            return _load_npz(npz, ids, require_all)
        raise ValueError(f"{source} has no residue embeddings")
    if not source.is_file():
        raise ValueError(f"{source} is not an embeddings file")
    return _load_npz(source, ids, require_all)


def load_embedding_files(
    paths: list[Path],
    ids: set[str] | None = None,
) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for path in paths:
        chunk = load_residue_embeddings(path, ids=ids, require_all=False)
        for identifier, value in chunk.items():
            if identifier in arrays:
                raise ValueError(f"duplicate embedding id {identifier!r} in {path}")
            arrays[identifier] = value
    if ids is not None:
        missing = ids - set(arrays)
        if missing:
            raise ValueError("missing embedding ids: " + ", ".join(sorted(missing)[:8]))
    return arrays


def write_residue_embeddings(outdir: Path, embeddings: dict[str, np.ndarray]) -> None:
    vectors, offsets, identifiers = _stack_embeddings(embeddings)
    np.save(outdir / VECTORS_NPY, vectors)
    np.save(outdir / OFFSETS_NPY, offsets)
    (outdir / IDS_TXT).write_text("\n".join(identifiers) + ("\n" if identifiers else ""))
    np.savez_compressed(
        outdir / NPZ_NAME,
        **{
            PACKED_VECTORS_KEY: vectors,
            PACKED_OFFSETS_KEY: offsets,
            PACKED_IDS_KEY: np.asarray(identifiers),
        },
    )


def write_records(
    outdir: Path,
    embeddings: dict[str, np.ndarray],
    records: list[tuple[str, str, str]],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    write_residue_embeddings(outdir, embeddings)
    with (outdir / "records.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "sequence", "secondary_structure"])
        writer.writerows(records)
