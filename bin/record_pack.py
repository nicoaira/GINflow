#!/usr/bin/env python3
"""Pack original residue embeddings and graph records for searchable databases."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np


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


def write_records(
    outdir: Path,
    embeddings: dict[str, np.ndarray],
    records: list[tuple[str, str, str]],
) -> None:
    np.savez_compressed(outdir / "embeddings.npz", **embeddings)
    with (outdir / "records.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "sequence", "secondary_structure"])
        writer.writerows(records)
