#!/usr/bin/env python3
"""Convert benchmark window selections into a GINflow query structures TSV."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


REQUIRED_STRUCTURE_COLUMNS = {"transcript_id", "sequence", "secondary_structure"}
REQUIRED_QUERY_COLUMNS = {"transcript_id", "window_offset"}


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = [
            {key: (value or "") for key, value in row.items()}
            for row in reader
            if any(value for value in row.values())
        ]
    return fieldnames, rows


def parse_nonnegative_int(value: str, field: str, row_number: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"row {row_number}: {field} must be a non-negative integer") from exc
    if parsed < 0:
        raise ValueError(f"row {row_number}: {field} must be a non-negative integer")
    return parsed


def convert(
    structures_path: Path,
    selections_path: Path,
    output_path: Path,
    window_size: int,
    full_molecules: bool = False,
) -> dict[str, object]:
    if window_size < 1:
        raise ValueError("--window-size must be positive")

    structure_columns, structure_rows = read_table(structures_path)
    missing = REQUIRED_STRUCTURE_COLUMNS - set(structure_columns)
    if missing:
        raise ValueError(f"{structures_path} is missing columns: {', '.join(sorted(missing))}")

    records: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(structure_rows, start=2):
        identifier = row["transcript_id"].strip()
        sequence = row["sequence"].strip().upper().replace("T", "U")
        structure = row["secondary_structure"].strip()
        if not identifier:
            raise ValueError(f"{structures_path} row {row_number}: transcript_id is empty")
        if identifier in records:
            raise ValueError(f"{structures_path} row {row_number}: duplicate transcript_id {identifier!r}")
        if len(sequence) != len(structure):
            raise ValueError(
                f"{structures_path} row {row_number}: sequence and secondary_structure lengths differ"
            )
        records[identifier] = {**row, "sequence": sequence, "secondary_structure": structure}

    query_columns, query_rows = read_table(selections_path)
    missing = REQUIRED_QUERY_COLUMNS - set(query_columns)
    if missing:
        raise ValueError(f"{selections_path} is missing columns: {', '.join(sorted(missing))}")

    provenance_columns = [
        column
        for column in query_columns
        if column not in structure_columns and column not in {"start", "end"}
    ]
    output_columns = [*structure_columns]
    if not full_molecules:
        for column in ("start", "end", *provenance_columns):
            if column not in output_columns:
                output_columns.append(column)

    output_rows: list[dict[str, str]] = []
    seen_identifiers: set[str] = set()
    for row_number, selection in enumerate(query_rows, start=2):
        identifier = selection["transcript_id"].strip()
        if identifier not in records:
            raise ValueError(f"{selections_path} row {row_number}: unknown transcript_id {identifier!r}")
        start = parse_nonnegative_int(selection["window_offset"], "window_offset", row_number)
        end = start + window_size
        record = records[identifier]
        if end > len(record["sequence"]):
            raise ValueError(
                f"{selections_path} row {row_number}: window [{start}, {end}) is outside "
                f"{identifier!r} length {len(record['sequence'])}"
            )

        if full_molecules:
            if identifier not in seen_identifiers:
                output_rows.append({column: record.get(column, "") for column in structure_columns})
                seen_identifiers.add(identifier)
        else:
            output_row = {column: record.get(column, "") for column in structure_columns}
            output_row["start"] = str(start)
            output_row["end"] = str(end)
            output_row.update({column: selection.get(column, "") for column in provenance_columns})
            output_rows.append(output_row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    return {
        "source_records": len(records),
        "query_rows": len(query_rows),
        "output_records": len(output_rows),
        "window_size": window_size,
        "full_molecules": full_molecules,
        "output": str(output_path),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures", type=Path, required=True, help="Full-molecule structures TSV")
    parser.add_argument("--selections", type=Path, required=True, help="Query selection TSV")
    parser.add_argument("--output", type=Path, required=True, help="Pipeline-compatible query structures TSV")
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument(
        "--full-molecules",
        action="store_true",
        help="Write one deduplicated full-molecule query per selected transcript instead of sliced subjects",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = convert(args.structures, args.selections, args.output, args.window_size, args.full_molecules)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
