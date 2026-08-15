#!/usr/bin/env python3
"""Merge per-query alignment TSVs/texts and rank by database E-value."""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

CLUSTER_HEAD = re.compile(
    r"^# cluster\s+(\S+)\s+(\S+)\s+vs\s+(\S+)",
)


def parse_float(value: str, default: float = math.inf) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def load_texts(paths: list[Path]) -> dict[tuple[str, str, str], str]:
    blocks: dict[tuple[str, str, str], list[str]] = {}
    key: tuple[str, str, str] | None = None
    current: list[str] = []
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        for line in path.read_text().splitlines():
            match = CLUSTER_HEAD.match(line)
            if match:
                if key is not None:
                    blocks[key] = current
                key = (match.group(1), match.group(2), match.group(3))
                current = [line]
                continue
            if key is not None:
                current.append(line)
        if key is not None:
            blocks[key] = current
            key = None
            current = []
    return {item: "\n".join(lines).rstrip() for item, lines in blocks.items()}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, nargs="*", default=[])
    parser.add_argument("--alignment-text", type=Path, nargs="*", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alignment-text-output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None
    for path in args.alignments:
        if not path.is_file():
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            if fieldnames is None:
                fieldnames = list(reader.fieldnames)
            rows.extend(row for row in reader if any(row.values()))
    if fieldnames is None:
        fieldnames = [
            "cluster_id", "query_id", "target_id", "score", "bit_score",
            "evalue", "evalue_pair", "query_start", "query_end",
            "target_start", "target_end", "query_length", "target_length",
            "match_count", "aligned_columns", "seed_count", "max_seed_score",
            "query_sequence", "query_structure", "target_sequence", "target_structure",
            "query_aligned", "target_aligned",
        ]

    rows.sort(key=lambda row: (
        parse_float(row.get("evalue")),
        -parse_float(row.get("score"), 0.0),
        row.get("query_id", ""),
        row.get("target_id", ""),
        str(row.get("cluster_id", "")),
    ))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

    if args.alignment_text_output:
        texts = load_texts(args.alignment_text)
        blocks = []
        for row in rows:
            key = (row.get("cluster_id", ""), row.get("query_id", ""), row.get("target_id", ""))
            block = texts.get(key)
            if block:
                blocks.append(block)
        args.alignment_text_output.write_text("\n\n".join(blocks) + ("\n" if blocks else ""))

    print(f"merged {len(rows)} alignments from {len(args.alignments)} tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
