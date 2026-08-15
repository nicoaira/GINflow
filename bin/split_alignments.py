#!/usr/bin/env python3
"""Split a TSV with a query_id column into one file per query, preserving order."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


def safe_name(value: str) -> str:
    cleaned = SAFE_NAME.sub("_", value).strip("_")
    return cleaned or "rna"


def unique_stem(query_id: str, used: dict[str, str]) -> str:
    stem = safe_name(query_id)
    if stem not in used or used[stem] == query_id:
        used[stem] = query_id
        return stem
    index = 2
    while f"{stem}_{index}" in used:
        index += 1
    alias = f"{stem}_{index}"
    used[alias] = query_id
    return alias


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "--alignments", dest="input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        with args.input.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or "query_id" not in reader.fieldnames:
                raise ValueError(f"{args.input} must have a query_id column")
            fieldnames = list(reader.fieldnames)
            groups: dict[str, list[dict[str, str]]] = {}
            order: list[str] = []
            for row in reader:
                query_id = row.get("query_id", "")
                if query_id not in groups:
                    groups[query_id] = []
                    order.append(query_id)
                groups[query_id].append(row)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    used: dict[str, str] = {}
    written = 0
    for query_id in order:
        path = args.outdir / f"{unique_stem(query_id, used)}.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(groups[query_id])
        written += 1
    print(f"wrote {written} query tables ({sum(len(rows) for rows in groups.values())} pairs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
