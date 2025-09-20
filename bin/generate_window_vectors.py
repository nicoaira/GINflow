#!/usr/bin/env python3
"""Build sliding window vectors from per-node embeddings.

The script expects a TSV with columns:
    <id_column>  node_index  embedding_vector
where ``embedding_vector`` is a comma-separated list of floats. Rows should be
sorted by node index; if not, they will be sorted prior to windowing.

Two outputs are produced under ``--output-dir``:
    - database_windows.npy / database_windows.tsv
    - query_windows.npy / query_windows.tsv

The numpy arrays contain float32 window vectors (each row normalized to unit
length). The TSV files provide metadata for each window so downstream modules
can map FAISS hits back to transcripts and positions.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterator, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate normalized sliding window vectors from node embeddings")
    parser.add_argument("--embeddings", required=True, help="TSV with node-level embeddings")
    parser.add_argument("--queries", required=True, help="CSV/TSV listing query IDs (expects a column named 'id' unless overridden)")
    parser.add_argument("--id-column", default="transcript_id", help="Identifier column shared between embeddings and query metadata")
    parser.add_argument("--node-index-column", default="node_index", help="Column giving 0-based node index per transcript")
    parser.add_argument("--embedding-column", default="embedding_vector", help="Column containing comma-separated embedding values")
    parser.add_argument("--window-size", type=int, required=True, help="Sliding window length in nodes")
    parser.add_argument("--stride", type=int, default=1, help="Stride between window starts (default: 1)")
    parser.add_argument("--output-dir", required=True, help="Directory for generated vectors + metadata")
    parser.add_argument("--normalize", action="store_true", help="Re-normalize each concatenated window (default true)")
    parser.add_argument("--no-normalize", dest="normalize", action="store_false")
    parser.set_defaults(normalize=True)
    parser.add_argument("--query-column", default="id", help="Column in --queries identifying query RNAs")
    parser.add_argument("--metadata-json", default=None, help="Optional path to write window statistics as JSON")
    parser.add_argument("--chunk-id", default=None, help="Optional chunk identifier for chunked execution")
    return parser.parse_args()


def load_query_ids(path: str, column: str) -> set[str]:
    ext = Path(path).suffix.lower()
    sep = "\t" if ext in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    if column not in df.columns:
        raise SystemExit(f"Query file {path} missing required column '{column}'")
    return set(df[column].dropna().astype(str))


def parse_embedding_string(series: pd.Series) -> np.ndarray:
    """Convert a Series of comma-delimited embedding strings into a 2D float array."""
    return np.stack(series.astype(str).str.split(',').map(lambda xs: np.asarray(xs, dtype=np.float32)))


def iter_windows(vecs: np.ndarray, window: int, stride: int) -> Iterator[Tuple[int, np.ndarray]]:
    limit = vecs.shape[0] - window + 1
    for start in range(0, max(0, limit), stride):
        end = start + window
        yield start, vecs[start:end]


def flatten_and_normalize(block: np.ndarray, normalize: bool) -> np.ndarray:
    flat = block.reshape(-1)
    if not normalize:
        return flat
    norm = np.linalg.norm(flat)
    if norm == 0:
        return flat
    return flat / norm


def count_windows(length: int, window_size: int, stride: int) -> int:
    limit = length - window_size + 1
    if limit <= 0:
        return 0
    return ((limit - 1) // stride) + 1


def main() -> None:
    args = parse_args()
    embeddings = pd.read_csv(args.embeddings, sep='\t', dtype={args.id_column: str})

    expected_cols = {args.id_column, args.node_index_column, args.embedding_column}
    missing = expected_cols - set(embeddings.columns)
    if missing:
        raise SystemExit(f"Embeddings file {args.embeddings} missing columns: {', '.join(sorted(missing))}")

    embeddings = embeddings.rename(columns={args.node_index_column: 'node_index', args.embedding_column: 'embedding_vector'})
    if embeddings.empty:
        raise SystemExit("No rows found in embeddings file; cannot build windows")

    embeddings[args.id_column] = embeddings[args.id_column].astype(str)
    try:
        embeddings['node_index'] = embeddings['node_index'].astype(int)
    except Exception as exc:  # pragma: no cover - defensive
        raise SystemExit(f"Failed to coerce node indices to integers: {exc}")
    embeddings = embeddings.sort_values([args.id_column, 'node_index'], kind='mergesort').reset_index(drop=True)

    try:
        vector_dim = len(str(embeddings['embedding_vector'].iloc[0]).split(','))
    except Exception as exc:  # pragma: no cover - defensive guard
        raise SystemExit(f"Failed to infer embedding dimension: {exc}")

    query_ids = load_query_ids(args.queries, args.query_column)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    window_dim = args.window_size * vector_dim

    meta_columns = [args.id_column, 'window_index', 'window_start', 'window_end', 'window_size', 'transcript_length']
    group_lengths = embeddings.groupby(args.id_column, sort=True)['node_index'].count()

    db_total = 0
    query_total = 0
    for transcript_id, length in group_lengths.items():
        window_count = count_windows(int(length), args.window_size, args.stride)
        if transcript_id in query_ids:
            query_total += window_count
        else:
            db_total += window_count
    db_path = outdir / 'database_windows.npy'
    query_path = outdir / 'query_windows.npy'

    if db_total:
        db_vectors = np.lib.format.open_memmap(db_path, mode='w+', dtype=np.float32, shape=(db_total, window_dim))
    else:
        np.save(db_path, np.zeros((0, window_dim), dtype=np.float32))
        db_vectors = None

    if query_total:
        query_vectors = np.lib.format.open_memmap(query_path, mode='w+', dtype=np.float32, shape=(query_total, window_dim))
    else:
        np.save(query_path, np.zeros((0, window_dim), dtype=np.float32))
        query_vectors = None

    db_meta_path = outdir / 'database_windows.tsv'
    query_meta_path = outdir / 'query_windows.tsv'
    db_meta_path.write_text('\t'.join(meta_columns) + '\n')
    query_meta_path.write_text('\t'.join(meta_columns) + '\n')

    db_index = 0
    query_index = 0

    with db_meta_path.open('a', encoding='utf-8') as db_meta_handle, query_meta_path.open('a', encoding='utf-8') as query_meta_handle:
        for transcript_id, group in embeddings.groupby(args.id_column, sort=True):
            group = group.sort_values('node_index')
            length = len(group)
            if length < args.window_size:
                continue
            vectors = parse_embedding_string(group['embedding_vector'])
            is_query = transcript_id in query_ids
            metadata_handle = query_meta_handle if is_query else db_meta_handle
            for window_idx, (start, block) in enumerate(iter_windows(vectors, window=args.window_size, stride=args.stride)):
                flat = flatten_and_normalize(block, args.normalize).astype(np.float32, copy=False)
                meta = {
                    args.id_column: transcript_id,
                    'window_index': int(window_idx),
                    'window_start': int(start),
                    'window_end': int(start + args.window_size),
                    'window_size': args.window_size,
                    'transcript_length': int(length),
                }
                line = '\t'.join(str(meta[col]) for col in meta_columns) + '\n'
                metadata_handle.write(line)
                if is_query:
                    if query_vectors is not None:
                        query_vectors[query_index] = flat
                    query_index += 1
                else:
                    if db_vectors is not None:
                        db_vectors[db_index] = flat
                    db_index += 1

    if db_vectors is not None:
        db_vectors.flush()
        del db_vectors
    if query_vectors is not None:
        query_vectors.flush()
        del query_vectors

    if db_total != db_index:
        raise SystemExit(f"Database window count mismatch (expected {db_total}, wrote {db_index})")
    if query_total != query_index:
        raise SystemExit(f"Query window count mismatch (expected {query_total}, wrote {query_index})")

    stats = {
        'database_windows': db_index,
        'query_windows': query_index,
        'window_size': args.window_size,
        'stride': args.stride,
        'queries': len(query_ids),
        'transcripts': len(group_lengths),
    }
    if args.chunk_id is not None:
        stats['chunk_id'] = args.chunk_id
    if args.metadata_json:
        Path(args.metadata_json).write_text(json.dumps(stats, indent=2))

    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
