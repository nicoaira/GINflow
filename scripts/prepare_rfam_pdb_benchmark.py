#!/usr/bin/env python3
"""Prepare the Rfam/PDB family-retrieval benchmark used by GINflow.

The script selects feasible families from the bundled Rfam/PDB molecule dumps,
builds and calibrates covariance models from the current Rfam seed for every
selected family and its clanmates, removes significant hits from the two
Rouskin backgrounds, and appends twelve labelled family mates per family.

The query TSVs contain the fields GINFINITY needs; cleaned background TSVs add
an ``rfam`` provenance column. The selection and cleaning audit trails live
beside them in the output directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


REQUIRED_COLUMNS = ("transcript_id", "sequence", "secondary_structure")
RFAM_TRANSCRIPT_RE = re.compile(r"^(RF\d{5})")
BIN_ORDER = ("gt70", "55_70", "lt55")
BIN_LABELS = {
    "gt70": ">70%",
    "55_70": "55-70%",
    "lt55": "<55%",
}
RFAM_SEED_URL = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz"
RFAM_CLANIN_URL = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"


@dataclass(frozen=True)
class SeedBlock:
    accession: str
    name: str
    text: str


@dataclass
class Molecule:
    source_id: str
    sequence: str
    structure: str
    core_start: int
    core_end: int
    columns: list[int]


@dataclass
class FamilySelection:
    family: str
    name: str
    clan: str
    clan_members: list[str]
    query: Molecule
    query_id: str
    query_index: int
    pids: np.ndarray
    available: dict[str, int]
    targets: list[tuple[str, int, Molecule, float]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=Path("."),
        help="Directory containing one molecules.jsonl directory per Rfam family.",
    )
    parser.add_argument(
        "--background-6k",
        type=Path,
        default=Path("tests/data/rouskin_sample_6k.tsv"),
    )
    parser.add_argument(
        "--background-30k",
        type=Path,
        default=Path("tests/data/rouskin_sample_30k.tsv"),
    )
    parser.add_argument(
        "--rfam-seed",
        type=Path,
        required=True,
        help="Current Rfam.seed or Rfam.seed.gz used to build the CMs.",
    )
    parser.add_argument(
        "--clanin",
        type=Path,
        required=True,
        help="Current Rfam.clanin used to expand selected families to clanmates.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("tests/data/rfam_pdb_benchmark"),
    )
    parser.add_argument("--n-families", type=int, default=20)
    parser.add_argument(
        "--cores",
        type=int,
        default=14,
        help="Maximum total CPU workers used for CM calibration/search (1-14).",
    )
    parser.add_argument(
        "--calibration-length-mb",
        type=float,
        default=0.01,
        help="Random sequence length passed to cmcalibrate (-L).",
    )
    parser.add_argument("--calibration-seed", type=int, default=42)
    parser.add_argument("--evalue", type=float, default=0.01)
    parser.add_argument(
        "--cmscan-max",
        action="store_true",
        help="Disable cmscan heuristic filters. This is much slower than default cmscan.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        help="Keep intermediate seeds, CMs, calibration logs, and cmscan output here.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow replacing files in an existing output directory.",
    )
    return parser.parse_args()


def numeric_family_id(accession: str) -> int:
    try:
        return int(accession.removeprefix("RF"))
    except ValueError:
        return math.inf


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_stockholm_blocks(path: Path) -> dict[str, SeedBlock]:
    opener = gzip.open if path.name.endswith(".gz") else open
    blocks: dict[str, SeedBlock] = {}
    current: list[str] = []
    with opener(path, "rt", encoding="latin-1") as handle:
        for line in handle:
            if line.startswith("# STOCKHOLM"):
                if current:
                    raise ValueError(f"unterminated Stockholm block before {path}")
                current = [line]
                continue
            if not current:
                continue
            current.append(line)
            if line.startswith("#=GF AC"):
                accession = line.split()[2].split(".", 1)[0]
            if line.startswith("#=GF ID"):
                name = line.split()[2]
            if line.strip() == "//":
                if "accession" not in locals() or "name" not in locals():
                    raise ValueError(f"Stockholm block missing AC or ID in {path}")
                if accession in blocks:
                    raise ValueError(f"duplicate Rfam accession {accession} in {path}")
                blocks[accession] = SeedBlock(accession, name, "".join(current))
                current = []
                del accession, name
    if current:
        raise ValueError(f"unterminated Stockholm block at end of {path}")
    return blocks


def read_clanin(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    name_to_clan: dict[str, str] = {}
    clan_to_names: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or not fields[0].startswith("CL"):
            continue
        clan, names = fields[0], fields[1:]
        clan_to_names[clan] = names
        for name in names:
            if name in name_to_clan and name_to_clan[name] != clan:
                raise ValueError(f"family {name} occurs in multiple clans")
            name_to_clan[name] = clan
    return name_to_clan, clan_to_names


def validate_structure(sequence: str, structure: str, identifier: str) -> None:
    if len(sequence) != len(structure):
        raise ValueError(f"{identifier}: sequence/structure length mismatch")
    if any(base.upper() not in "ACGUT" for base in sequence):
        raise ValueError(f"{identifier}: unsupported base in sequence")
    depth = 0
    for char in structure:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise ValueError(f"{identifier}: unbalanced structure")
        elif char != ".":
            raise ValueError(f"{identifier}: structure is not dot-bracket .()")
    if depth:
        raise ValueError(f"{identifier}: unbalanced structure")


def load_family(path: Path) -> list[Molecule]:
    molecules: list[Molecule] = []
    seen: set[str] = set()
    with (path / "molecules.jsonl").open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            record = json.loads(line)
            source_id = str(record["id"])
            if source_id in seen:
                raise ValueError(f"{path}: duplicate molecule id {source_id}")
            seen.add(source_id)
            sequence = str(record["sequence"])
            structure = str(record["structure"])
            # GINFINITY's public input contract is strict A/C/G/U. Keep the
            # source record untouched in the audit's family count, but do not
            # make an invalid molecule part of a runnable benchmark.
            if any(base.upper() not in "ACGUT" for base in sequence):
                continue
            validate_structure(sequence, structure, source_id)
            core_start = int(record["core_start"]) - 1
            core_length = int(record["core_length"])
            core_end = core_start + core_length
            columns = [int(value) for value in record["columns"]]
            if not (0 <= core_start < core_end <= len(sequence)):
                raise ValueError(f"{path}:{line_number}: invalid core coordinates")
            if len(columns) != core_length or len(set(columns)) != len(columns):
                raise ValueError(f"{source_id}: columns do not describe the core")
            core = sequence[core_start:core_end]
            if len(core) != len(columns):
                raise ValueError(f"{source_id}: core/columns length mismatch")
            molecules.append(
                Molecule(source_id, sequence, structure, core_start, core_end, columns)
            )
    if not molecules:
        raise ValueError(f"no molecules found in {path}")
    return molecules


def identity_matrix(molecules: Sequence[Molecule]) -> np.ndarray:
    width = max(max(molecule.columns) for molecule in molecules) + 1
    aligned = np.full((len(molecules), width), "", dtype="U1")
    for row, molecule in enumerate(molecules):
        core = molecule.sequence[molecule.core_start : molecule.core_end].upper().replace("T", "U")
        aligned[row, np.asarray(molecule.columns, dtype=int)] = list(core)

    present = aligned != ""
    result = np.zeros((len(molecules), len(molecules)), dtype=np.float64)
    for row in range(len(molecules)):
        shared = present & present[row]
        denominator = shared.sum(axis=1)
        matches = ((aligned == aligned[row]) & shared).sum(axis=1)
        result[row] = np.divide(
            matches,
            denominator,
            out=np.zeros(len(molecules), dtype=np.float64),
            where=denominator > 0,
        )
        result[row, row] = -1.0
    return result


def identity_bin(value: float) -> str:
    if value > 0.70:
        return "gt70"
    if value >= 0.55:
        return "55_70"
    return "lt55"


def evenly_spread(indices: Sequence[int], pids: np.ndarray, molecules: Sequence[Molecule], count: int) -> list[int]:
    ordered = sorted(indices, key=lambda index: (float(pids[index]), molecules[index].source_id))
    if len(ordered) <= count:
        return ordered
    positions = [int(round(position * (len(ordered) - 1) / (count - 1))) for position in range(count)]
    return [ordered[position] for position in positions]


def select_family(
    family: str,
    name: str,
    clan: str,
    clan_members: list[str],
    molecules: list[Molecule],
) -> FamilySelection:
    pids = identity_matrix(molecules)
    candidates: list[tuple[tuple[int, int, int, str], int, dict[str, int]]] = []
    for query_index, molecule in enumerate(molecules):
        counts = Counter(identity_bin(float(pids[query_index, index])) for index in range(len(molecules)) if index != query_index)
        available = {bin_name: counts[bin_name] for bin_name in BIN_ORDER}
        slots = sum(min(available[bin_name], 4) for bin_name in BIN_ORDER)
        minimum = min(available.values())
        total = sum(available.values())
        # Maximize feasible slots, then robustness of the weakest bin, then total
        # available mates; source id is the stable final tie-breaker.
        key = (-slots, -minimum, -total, molecule.source_id)
        candidates.append((key, query_index, available))
    _, query_index, available = min(candidates)
    query = molecules[query_index]
    query_id = f"{family}__query"
    targets: list[tuple[str, int, Molecule, float]] = []
    for bin_name in BIN_ORDER:
        indices = [
            index
            for index in range(len(molecules))
            if index != query_index and identity_bin(float(pids[query_index, index])) == bin_name
        ]
        for rank, target_index in enumerate(evenly_spread(indices, pids[query_index], molecules, min(4, len(indices))), 1):
            target_id = f"{family}__target_{bin_name}_{rank:02d}"
            targets.append((target_id, rank, molecules[target_index], float(pids[query_index, target_index])))
    return FamilySelection(
        family,
        name,
        clan,
        clan_members,
        query,
        query_id,
        query_index,
        pids[query_index],
        available,
        targets,
    )


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or any(column not in reader.fieldnames for column in REQUIRED_COLUMNS):
            raise ValueError(f"{path}: required columns are {REQUIRED_COLUMNS}")
        rows = []
        for row in reader:
            if not row["transcript_id"]:
                raise ValueError(f"{path}: empty transcript_id")
            validate_structure(row["sequence"], row["secondary_structure"], row["transcript_id"])
            rows.append(row)
    if len({row["transcript_id"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: transcript_id values are not unique")
    return rows


def write_structure_table(
    path: Path,
    rows: Iterable[dict[str, str]],
    extra_columns: Sequence[str] = (),
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [*REQUIRED_COLUMNS, *extra_columns]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in fieldnames})


def write_rows(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def run_command(command: Sequence[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): {' '.join(command)}; see {log_path}")


def build_calibrated_models(
    seed_blocks: dict[str, SeedBlock],
    model_names: Sequence[str],
    workdir: Path,
    cores: int,
    calibration_length_mb: float,
    calibration_seed: int,
) -> tuple[Path, dict[str, SeedBlock]]:
    model_dir = workdir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    cm_paths: list[Path] = []
    selected_blocks: dict[str, SeedBlock] = {}
    for name in model_names:
        block = next(block for block in seed_blocks.values() if block.name == name)
        selected_blocks[name] = block
        sto_path = model_dir / f"{block.accession}.sto"
        cm_path = model_dir / f"{block.accession}.cm"
        sto_path.write_text(block.text, encoding="latin-1")
        run_command(
            ["cmbuild", "--hand", "-F", "-n", block.accession, str(cm_path), str(sto_path)],
            workdir / "logs" / f"cmbuild_{block.accession}.log",
        )
        cm_paths.append(cm_path)

    def calibrate(cm_path: Path) -> None:
        run_command(
            [
                "cmcalibrate",
                "--cpu",
                "1",
                "-L",
                str(calibration_length_mb),
                "--seed",
                str(calibration_seed),
                str(cm_path),
            ],
            workdir / "logs" / f"cmcalibrate_{cm_path.stem}.log",
        )

    with ThreadPoolExecutor(max_workers=cores) as executor:
        futures = [executor.submit(calibrate, cm_path) for cm_path in cm_paths]
        for future in as_completed(futures):
            future.result()

    combined = workdir / "rfam_selected_clans.cm"
    header = cm_paths[0].read_text(encoding="latin-1").splitlines(keepends=True)[0]
    with combined.open("w", encoding="latin-1") as output:
        for index, cm_path in enumerate(cm_paths):
            lines = cm_path.read_text(encoding="latin-1").splitlines(keepends=True)
            if index == 0:
                output.write("".join(lines))
            else:
                output.write(header)
                output.write("".join(lines[1:]))
    run_command(["cmpress", "-F", str(combined)], workdir / "logs" / "cmpress.log")
    return combined, selected_blocks


def write_union_fasta(backgrounds: Sequence[tuple[str, Path, list[dict[str, str]]]], path: Path) -> dict[str, dict[str, str]]:
    union: dict[str, dict[str, str]] = {}
    for _, _, rows in backgrounds:
        for row in rows:
            identifier = row["transcript_id"]
            existing = union.get(identifier)
            current = {"sequence": row["sequence"], "secondary_structure": row["secondary_structure"]}
            if existing is not None and existing != current:
                raise ValueError(f"background id {identifier} has conflicting records")
            union[identifier] = current
    with path.open("w", encoding="utf-8") as handle:
        for identifier, row in union.items():
            if any(char.isspace() for char in identifier) or ">" in identifier:
                raise ValueError(f"background id cannot be represented in FASTA: {identifier!r}")
            handle.write(f">{identifier}\n{row['sequence']}\n")
    return union


def parse_cmscan_hits(path: Path, evalue_threshold: float) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 17:
                raise ValueError(f"unexpected cmscan tblout row: {line.rstrip()}")
            evalue = float(fields[15])
            if evalue < evalue_threshold:
                hits.append(
                    {
                        "model_name": fields[0],
                        "model_accession": fields[1],
                        "transcript_id": fields[2],
                        "evalue": evalue,
                    }
                )
    return hits


def target_rows(selections: Sequence[FamilySelection]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for selection in selections:
        for target_id, _, molecule, _ in selection.targets:
            rows.append(
                {
                    "transcript_id": target_id,
                    "sequence": molecule.sequence,
                    "secondary_structure": molecule.structure,
                    "rfam": selection.family,
                }
            )
    return rows


def background_rfam(transcript_id: str) -> str:
    match = RFAM_TRANSCRIPT_RE.match(transcript_id)
    return match.group(1) if match else "NA"


def write_selection_files(outdir: Path, selections: Sequence[FamilySelection]) -> None:
    query_rows = [
        {
            "transcript_id": selection.query_id,
            "sequence": selection.query.sequence,
            "secondary_structure": selection.query.structure,
            "rfam": selection.family,
        }
        for selection in selections
    ]
    write_structure_table(outdir / "queries.tsv", query_rows, extra_columns=("rfam",))
    with (outdir / "queries_with_windows.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[*REQUIRED_COLUMNS, "rfam", "start", "end"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row, selection in zip(query_rows, selections):
            writer.writerow(
                {
                    **row,
                    "start": selection.query.core_start,
                    "end": selection.query.core_end,
                }
            )

    selection_rows: list[dict[str, object]] = []
    family_rows: list[dict[str, object]] = []
    for selection in selections:
        selected_counts = Counter(identity_bin(pid) for _, _, _, pid in selection.targets)
        family_rows.append(
            {
                "rfam_family": selection.family,
                "family_name": selection.name,
                "clan": selection.clan,
                "clan_members": ",".join(selection.clan_members),
                "query_id": selection.query_id,
                "query_source_id": selection.query.source_id,
                "query_core_start": selection.query.core_start,
                "query_core_end": selection.query.core_end,
                "family_molecules": len(selection.pids),
                "available_gt70": selection.available["gt70"],
                "available_55_70": selection.available["55_70"],
                "available_lt55": selection.available["lt55"],
                "selected_gt70": selected_counts["gt70"],
                "selected_55_70": selected_counts["55_70"],
                "selected_lt55": selected_counts["lt55"],
                "selected_targets": len(selection.targets),
            }
        )
        for target_id, rank, molecule, pid in selection.targets:
            selection_rows.append(
                {
                    "rfam_family": selection.family,
                    "family_name": selection.name,
                    "clan": selection.clan,
                    "query_id": selection.query_id,
                    "query_source_id": selection.query.source_id,
                    "target_id": target_id,
                    "target_source_id": molecule.source_id,
                    "identity": f"{pid:.8f}",
                    "identity_bin": BIN_LABELS[identity_bin(pid)],
                    "bin_rank": rank,
                    "target_core_start": molecule.core_start,
                    "target_core_end": molecule.core_end,
                }
            )
    write_rows(
        outdir / "family_summary.tsv",
        [
            "rfam_family", "family_name", "clan", "clan_members", "query_id", "query_source_id",
            "query_core_start", "query_core_end", "family_molecules", "available_gt70", "available_55_70",
            "available_lt55", "selected_gt70", "selected_55_70", "selected_lt55", "selected_targets",
        ],
        family_rows,
    )
    write_rows(
        outdir / "selections.tsv",
        [
            "rfam_family", "family_name", "clan", "query_id", "query_source_id", "target_id",
            "target_source_id", "identity", "identity_bin", "bin_rank", "target_core_start", "target_core_end",
        ],
        selection_rows,
    )


def write_model_manifest(
    outdir: Path,
    selected: Sequence[FamilySelection],
    seed_blocks: dict[str, SeedBlock],
) -> list[str]:
    selected_ids = {selection.family for selection in selected}
    names: set[str] = set()
    for selection in selected:
        names.update(selection.clan_members)
    rows = []
    for name in sorted(names, key=lambda value: numeric_family_id(seed_blocks_by_name(seed_blocks)[value].accession)):
        block = seed_blocks_by_name(seed_blocks)[name]
        rows.append(
            {
                "rfam_accession": block.accession,
                "family_name": name,
                "clan": next((selection.clan for selection in selected if name in selection.clan_members), f"SINGLETON:{name}"),
                "selected_family": "yes" if block.accession in selected_ids else "no",
            }
        )
    write_rows(outdir / "selected_models.tsv", ["rfam_accession", "family_name", "clan", "selected_family"], rows)
    return [row["family_name"] for row in rows]


def seed_blocks_by_name(seed_blocks: dict[str, SeedBlock]) -> dict[str, SeedBlock]:
    return {block.name: block for block in seed_blocks.values()}


def clean_backgrounds(
    outdir: Path,
    backgrounds: Sequence[tuple[str, Path, list[dict[str, str]]]],
    hits: Sequence[dict[str, object]],
    selections: Sequence[FamilySelection],
) -> None:
    by_identifier: dict[str, list[dict[str, object]]] = defaultdict(list)
    for hit in hits:
        by_identifier[str(hit["transcript_id"])].append(hit)
    targets = target_rows(selections)
    hit_rows: list[dict[str, object]] = []
    for identifier, identifier_hits in sorted(by_identifier.items()):
        for hit in sorted(identifier_hits, key=lambda value: (float(value["evalue"]), str(value["model_accession"]))):
            hit_rows.append(
                {
                    "transcript_id": identifier,
                    "model_name": hit["model_name"],
                    "model_accession": hit["model_accession"],
                    "evalue": f"{float(hit['evalue']):.6g}",
                }
            )
    write_rows(outdir / "cmscan_hits.tsv", ["transcript_id", "model_name", "model_accession", "evalue"], hit_rows)

    removed_rows = []
    for identifier, identifier_hits in sorted(by_identifier.items()):
        best = min(identifier_hits, key=lambda value: float(value["evalue"]))
        removed_rows.append(
            {
                "transcript_id": identifier,
                "hit_count": len(identifier_hits),
                "best_model": best["model_accession"],
                "best_evalue": f"{float(best['evalue']):.6g}",
            }
        )
    write_rows(outdir / "removed_background_ids.tsv", ["transcript_id", "hit_count", "best_model", "best_evalue"], removed_rows)

    summaries = []
    for label, path, rows in backgrounds:
        retained = [row for row in rows if row["transcript_id"] not in by_identifier]
        output_rows = [
            {
                "transcript_id": row["transcript_id"],
                "sequence": row["sequence"],
                "secondary_structure": row["secondary_structure"],
                "rfam": background_rfam(row["transcript_id"]),
            }
            for row in retained
        ] + targets
        output_path = outdir / f"{path.stem}.cleaned.tsv"
        write_structure_table(output_path, output_rows, extra_columns=("rfam",))
        summaries.append(
            {
                "background": label,
                "source": str(path),
                "input_rows": len(rows),
                "cmscan_hit_rows": sum(1 for row in rows if row["transcript_id"] in by_identifier),
                "removed_rows": len(rows) - len(retained),
                "retained_background_rows": len(retained),
                "appended_target_rows": len(targets),
                "output_rows": len(output_rows),
                "output": str(output_path),
            }
        )
    write_rows(
        outdir / "cleaning_summary.tsv",
        [
            "background", "source", "input_rows", "cmscan_hit_rows", "removed_rows", "retained_background_rows",
            "appended_target_rows", "output_rows", "output",
        ],
        summaries,
    )


def write_readme(
    outdir: Path,
    selected: Sequence[FamilySelection],
    model_names: Sequence[str],
    backgrounds: Sequence[tuple[str, Path, list[dict[str, str]]]],
    union_count: int,
    args: argparse.Namespace,
) -> None:
    total_targets = sum(len(selection.targets) for selection in selected)
    family_lines = []
    for selection in selected:
        counts = Counter(identity_bin(pid) for _, _, _, pid in selection.targets)
        family_lines.append(
            f"| {selection.family} | {selection.name} | {selection.clan} | {selection.query_id} | "
            f"{counts['gt70']}/{counts['55_70']}/{counts['lt55']} |"
        )
    background_lines = []
    for label, path, rows in backgrounds:
        background_lines.append(f"| {label} | `{path}` | {len(rows)} | `{path.stem}.cleaned.tsv` |")
    text = f"""# GINflow Rfam/PDB retrieval benchmark

This directory contains a 20-family retrieval benchmark built from
`research-gine-rna-encoder/dataset/rfam_pdb_families` and the two Rouskin
background tables. It is generated by
[`scripts/prepare_rfam_pdb_benchmark.py`](../../../scripts/prepare_rfam_pdb_benchmark.py).

## Files

- `queries.tsv` — one query per family with only GINflow's mandatory columns.
- `queries_with_windows.tsv` — the same queries plus `start`/`end` core
  windows. Coordinates are 0-based, half-open, as required by GINflow.
- `rouskin_sample_6k.cleaned.4k.tsv` and `rouskin_sample_30k.cleaned.17k.tsv` —
  clan-cleaned backgrounds followed by the twelve selected mates per family.
- `family_summary.tsv` and `selections.tsv` — query, target, identity-bin, and
  source-record audit trails.
- `selected_models.tsv`, `cmscan_hits.tsv`, `removed_background_ids.tsv`, and
  `cleaning_summary.tsv` — the Rfam model and background-cleaning audit.
- `run_metadata.json` — source checksums and exact preparation parameters.

## Selection protocol

Families were sorted by numeric Rfam accession and the first {len(selected)}
families with a query supporting all three requested bins were selected. For
each family, the query maximizes the number of available four-member slots,
then the weakest-bin count, then total available mates; the source molecule ID
breaks ties. Four targets per bin were selected at evenly spaced identity ranks
within that bin. Identity is computed over shared Rfam seed match columns:
`>70%`, `55-70%` (including both boundaries), and `<55%`.

All selected families have 12 targets in this build ({total_targets} targets
total), with four targets in each bin. The source Rfam/PDB molecule chosen as
the query is not appended to either background.

| family | name | clan | query id | selected targets: >70 / 55-70 / <55 |
|---|---|---|---|---|
{os.linesep.join(family_lines)}

## Clan-cleaning protocol

The current Rfam seed was split into one Stockholm alignment per selected
family and per member of its Rfam clan. Each alignment was built with
`cmbuild --hand`, calibrated with `cmcalibrate` before search, combined into a
pressed CM database, and searched with `cmscan` at E-value `< {args.evalue}`.
The calibration used `-L {args.calibration_length_mb}` Mb and seed
`{args.calibration_seed}`; this is recorded because full default calibration is
substantially slower. The scan used the default Infernal heuristic filters
unless `--cmscan-max` was supplied. The union of both backgrounds contained
{union_count} unique molecules and was scanned once; the 6k table is an exact
subset of the 30k table, so this is equivalent to scanning each table.

| background | source | input rows | cleaned output |
|---|---|---:|---|
{os.linesep.join(background_lines)}

Any input molecule with at least one reported hit below the threshold to any
selected-family or clanmate CM was removed. The selected 12 family mates are
then appended intentionally as positives. No sequences or structures were
refolded or otherwise modified.

## GINflow usage

Build a database from either cleaned background:

```bash
nextflow run . -profile docker \\
  --input tests/data/rfam_pdb_benchmark/rouskin_sample_30k.cleaned.17k.tsv \\
  --outdir results/rfam-pdb-30k
```

Search it with the mandatory query sheet:

```bash
nextflow run . -profile docker \\
  --query tests/data/rfam_pdb_benchmark/queries.tsv \\
  --database results/rfam-pdb-30k/index \\
  --outdir results/rfam-pdb-30k-search
```

Use `queries_with_windows.tsv` for a core-only query run. Its `start` and
`end` values select the Rfam-aligned part without the molecule's flanks.

## Rfam sources

- {RFAM_SEED_URL}
- {RFAM_CLANIN_URL}

See `run_metadata.json` for SHA-256 checksums and the exact command settings.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    if not 1 <= args.cores <= 14:
        raise ValueError("--cores must be between 1 and 14")
    if args.calibration_length_mb <= 0 or args.evalue <= 0:
        raise ValueError("calibration length and E-value must be positive")
    if args.n_families <= 0:
        raise ValueError("--n-families must be positive")
    for path in [args.dataset_root, args.background_6k, args.background_30k, args.rfam_seed, args.clanin]:
        if not path.exists():
            raise FileNotFoundError(path)

    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.force:
        raise FileExistsError(f"{args.outdir} is not empty; use --force to replace generated files")
    args.outdir.mkdir(parents=True, exist_ok=True)
    for child in args.outdir.iterdir():
        if child.is_file():
            child.unlink()
        elif child.is_dir():
            shutil.rmtree(child)

    seed_blocks = read_stockholm_blocks(args.rfam_seed)
    name_to_clan, clan_to_names = read_clanin(args.clanin)
    family_dirs = sorted((path for path in args.dataset_root.iterdir() if path.is_dir()), key=lambda path: numeric_family_id(path.name))
    family_info = {}
    all_molecules = {}
    for family_dir in family_dirs:
        metadata = json.loads((family_dir / "family.json").read_text(encoding="utf-8"))
        family = metadata["family"]
        name = metadata["name"]
        molecules = load_family(family_dir)
        family_info[family] = (name, molecules)
        all_molecules[family] = molecules

    feasible: list[FamilySelection] = []
    for family in sorted(family_info, key=numeric_family_id):
        name, molecules = family_info[family]
        clan = name_to_clan.get(name, f"SINGLETON:{name}")
        clan_members = sorted(clan_to_names.get(clan, [name]), key=lambda value: numeric_family_id(seed_blocks_by_name(seed_blocks)[value].accession))
        selection = select_family(family, name, clan, clan_members, molecules)
        if all(selection.available[bin_name] >= 4 for bin_name in BIN_ORDER):
            feasible.append(selection)
    selected = feasible[: args.n_families]
    if not selected:
        raise RuntimeError("no feasible families found")

    backgrounds = [
        ("rouskin_6k", args.background_6k, read_table(args.background_6k)),
        ("rouskin_30k", args.background_30k, read_table(args.background_30k)),
    ]
    all_background_ids = [row["transcript_id"] for _, _, rows in backgrounds for row in rows]
    if len(all_background_ids) != len(set(all_background_ids)):
        # The shipped 6k table is a subset of 30k. Duplicate IDs are allowed
        # across the two inputs, but must represent the same sequence/structure.
        pass

    work_context = tempfile.TemporaryDirectory(prefix="rfam-pdb-benchmark-")
    if args.workdir:
        workdir = args.workdir
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        workdir = Path(work_context.name)
    try:
        model_names = sorted(
            {name for selection in selected for name in selection.clan_members},
            key=lambda value: numeric_family_id(seed_blocks_by_name(seed_blocks)[value].accession),
        )
        combined_cm, selected_blocks = build_calibrated_models(
            seed_blocks,
            model_names,
            workdir,
            args.cores,
            args.calibration_length_mb,
            args.calibration_seed,
        )
        union_fasta = workdir / "background_union.fa"
        union = write_union_fasta(backgrounds, union_fasta)
        tblout = workdir / "cmscan.tblout"
        scan_command = [
            "cmscan",
            "--cpu",
            str(args.cores),
            "--noali",
            "--notextw",
            "--tblout",
            str(tblout),
            "-E",
            str(args.evalue),
            "--incE",
            str(args.evalue),
        ]
        if args.cmscan_max:
            scan_command.append("--max")
        scan_command.extend([str(combined_cm), str(union_fasta)])
        run_command(scan_command, workdir / "logs" / "cmscan.log")
        hits = parse_cmscan_hits(tblout, args.evalue)

        write_selection_files(args.outdir, selected)
        write_model_manifest(args.outdir, selected, seed_blocks)
        clean_backgrounds(args.outdir, backgrounds, hits, selected)

        metadata = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "script": "scripts/prepare_rfam_pdb_benchmark.py",
            "rfam_seed_url": RFAM_SEED_URL,
            "rfam_clanin_url": RFAM_CLANIN_URL,
            "rfam_seed": str(args.rfam_seed),
            "rfam_seed_sha256": sha256(args.rfam_seed),
            "rfam_clanin": str(args.clanin),
            "rfam_clanin_sha256": sha256(args.clanin),
            "dataset_root": str(args.dataset_root),
            "backgrounds": {label: {"path": str(path), "sha256": sha256(path), "rows": len(rows)} for label, path, rows in backgrounds},
            "selected_families": [selection.family for selection in selected],
            "selected_family_count": len(selected),
            "selected_model_count": len(model_names),
            "selected_model_names": model_names,
            "background_union_unique_molecules": len(union),
            "cmscan_command": scan_command,
            "cmscan_evalue_threshold": args.evalue,
            "cmscan_mode": "max" if args.cmscan_max else "default",
            "cores": args.cores,
            "cmcalibrate_length_mb": args.calibration_length_mb,
            "cmcalibrate_seed": args.calibration_seed,
            "cmscan_reported_hit_rows_below_threshold": len(hits),
            "cmscan_removed_unique_molecules": len({str(hit["transcript_id"]) for hit in hits}),
            "target_count": sum(len(selection.targets) for selection in selected),
        }
        (args.outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        write_readme(args.outdir, selected, model_names, backgrounds, len(union), args)
        print(json.dumps({
            "selected_families": [selection.family for selection in selected],
            "selected_models": len(model_names),
            "union_background_molecules": len(union),
            "cmscan_hit_rows": len(hits),
            "cmscan_removed_molecules": len({str(hit['transcript_id']) for hit in hits}),
            "targets_added_per_background": sum(len(selection.targets) for selection in selected),
            "outdir": str(args.outdir),
        }, indent=2))
    finally:
        if args.workdir is None:
            work_context.cleanup()


if __name__ == "__main__":
    main()
