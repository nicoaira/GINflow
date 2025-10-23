#!/usr/bin/env python3
"""Split merged query vectors and metadata into separate files per query.

This enables parallelization of downstream processes (FAISS querying, clustering,
alignment) across multiple query structures.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split query vectors and metadata by query ID"
    )
    parser.add_argument(
        "--query-vectors", required=True, help="Input .npy file with all query vectors"
    )
    parser.add_argument(
        "--query-metadata",
        required=True,
        help="Input .tsv file with query window metadata",
    )
    parser.add_argument(
        "--id-column",
        default="transcript_id",
        help="Column identifying query transcripts",
    )
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for split files"
    )
    return parser.parse_args()


def safe_filename(query_id: str) -> str:
    """Convert query ID to filesystem-safe string."""
    return "".join(c if c.isalnum() or c in ("_", "-") else "_" for c in query_id)


def main() -> None:
    args = parse_args()

    # Load query metadata
    query_meta = pd.read_csv(args.query_metadata, sep="\t", dtype={args.id_column: str})

    if query_meta.empty:
        print("No query windows found. Exiting.")
        return

    # Load query vectors
    query_vectors = np.load(args.query_vectors)

    if len(query_meta) != len(query_vectors):
        raise SystemExit(
            f"Mismatch: {len(query_meta)} metadata rows vs {len(query_vectors)} vectors"
        )

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Group by query ID
    id_col = args.id_column
    grouped = query_meta.groupby(id_col, sort=False)

    query_info: List[Dict[str, str]] = []

    for query_id, group in grouped:
        safe_id = safe_filename(str(query_id))
        indices = group.index.tolist()

        # Extract vectors for this query
        query_vecs = query_vectors[indices]

        # Write vectors
        vec_path = output_dir / f"query_{safe_id}_vectors.npy"
        np.save(vec_path, query_vecs)

        # Write metadata
        meta_path = output_dir / f"query_{safe_id}_metadata.tsv"
        # Reset window indices to be 0-based for each query
        group_copy = group.copy()
        group_copy["window_index"] = range(len(group_copy))
        group_copy.to_csv(meta_path, sep="\t", index=False)

        query_info.append(
            {
                "query_id": str(query_id),
                "safe_id": safe_id,
                "num_windows": len(group),
                "vector_file": str(vec_path.name),
                "metadata_file": str(meta_path.name),
            }
        )

        print(
            f"Query '{query_id}' -> {len(group)} windows "
            f"({vec_path.name}, {meta_path.name})"
        )

    # Write summary JSON
    summary_path = output_dir / "query_split_manifest.json"
    summary_path.write_text(
        json.dumps(
            {"num_queries": len(query_info), "queries": query_info}, indent=2
        )
        + "\n"
    )

    print(f"\nSplit {len(query_info)} queries into {output_dir}")
    print(f"Manifest: {summary_path}")


if __name__ == "__main__":
    main()
