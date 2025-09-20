#!/usr/bin/env python3
"""Merge node embedding chunk TSV files into a single table.

The manifest is a tab-delimited file with two columns:
    - chunk identifier (used for deterministic ordering)
    - path to the chunk TSV file
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge node embedding TSV chunks")
    parser.add_argument("--manifest", required=True, help="Tab-delimited manifest describing the chunk files")
    parser.add_argument("--output", required=True, help="Path to write the merged TSV")
    return parser.parse_args()


def load_manifest(path: str) -> List[Tuple[str, Path]]:
    manifest_path = Path(path)
    if not manifest_path.exists():
        raise SystemExit(f"Manifest not found: {manifest_path}")
    entries: List[Tuple[str, Path]] = []
    for line in manifest_path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split('\t', 1)
        if len(parts) != 2:
            raise SystemExit(f"Malformed manifest line: {line}")
        chunk_id, chunk_path = parts
        entries.append((chunk_id, Path(chunk_path)))
    return entries


def merge_chunks(manifest: List[Tuple[str, Path]], output_path: Path) -> None:
    if not manifest:
        raise SystemExit("Manifest contained no entries")

    entries = []
    for chunk_id, chunk_path in manifest:
        chunk_id = str(chunk_id)
        if not chunk_path.exists():
            raise SystemExit(f"Embedding chunk not found for {chunk_id}: {chunk_path}")
        entries.append({"chunk_id": chunk_id, "path": chunk_path})

    entries.sort(key=lambda item: item["chunk_id"])

    header_written = False
    header_reference = None
    with output_path.open("w", encoding="utf-8") as out_handle:
        for item in entries:
            chunk_path = item["path"]
            with chunk_path.open("r", encoding="utf-8") as src:
                header = src.readline()
                if not header:
                    continue
                if not header_written:
                    out_handle.write(header)
                    header_written = True
                    header_reference = header.strip()
                else:
                    # ensure headers are consistent when present
                    if header_reference and header.strip() and header.strip() != header_reference:
                        raise SystemExit(
                            "Unexpected header mismatch between embedding chunks: "
                            f"chunk {item['chunk_id']} header '{header.strip()}' != '{header_reference}'"
                        )
                for line in src:
                    out_handle.write(line)

    if not header_written:
        # file was empty, create zero-length file with header? just touch
        output_path.touch()


def main() -> None:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    merge_chunks(manifest, Path(args.output))


if __name__ == "__main__":
    main()
