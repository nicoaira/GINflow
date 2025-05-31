#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd
import faiss

# Maximum index memory usage (bytes)
MAX_INDEX_BYTES = 16 * 1024**3

def parse_args():
    p = argparse.ArgumentParser(
        description="Build a FAISS index from embeddings.tsv, using memmap + chunked adds"
    )
    p.add_argument('--input',        required=True,
                   help="Path to embeddings.tsv (must contain 'embedding_vector', 'window_start','window_end','seq_len')")
    p.add_argument('--id-column',    required=True,
                   help="Name of the transcript/window ID column")
    p.add_argument('--query',        required=True,
                   help="ID value to treat as “query” (will be excluded from DB)")
    p.add_argument('--index-path',   required=True,
                   help="Where to write the FAISS index file")
    p.add_argument('--mapping-path', required=True,
                   help="Where to write the DB metadata (with seq_len)")
    return p.parse_args()

def infer_size_and_dim(tsv_path, vec_col='embedding_vector'):
    """Read one line to discover D, then count total rows for N."""
    with open(tsv_path, 'r') as f:
        header = f.readline().rstrip('\n').split('\t')
        sample = f.readline().rstrip('\n').split('\t')
        vec_str = sample[ header.index(vec_col) ]
        D = len(vec_str.split(','))
    # count total data lines
    with open(tsv_path, 'r') as f:
        N = sum(1 for _ in f) -1
    return N, D

def main():
    args = parse_args()

    # 1) figure out N,D
    N, D = infer_size_and_dim(args.input)
    est_index_bytes = N * D * 4
    print(f"Found {N:,} vectors of dimension {D} → index will be ~{est_index_bytes/1024**3:.1f} GiB")
    if est_index_bytes > MAX_INDEX_BYTES:
        print("Warning: estimated index size exceeds 16 GiB; you may run out of RAM.")

    # 2) build a disk‐backed memmap
    memmap_path = args.index_path + '.memmap.npy'
    print(f"Creating memmap at {memmap_path}")
    mm = np.memmap(memmap_path, dtype='float32', mode='w+',
                   shape=(N, D))

    # 3) read in chunks, fill memmap & collect metadata
    meta_cols = [ args.id_column, 'window_start','window_end','seq_len','embedding_vector' ]
    reader = pd.read_csv(
        args.input,
        sep='\t',
        usecols=meta_cols,
        dtype={args.id_column:str},
        chunksize=10_000
    )

    idx = 0
    id_list, ws_list, we_list, sl_list = [], [], [], []
    for chunk in reader:
        M = len(chunk)
        # parse embedding strings → float32 array
        emb_mat = np.stack(
            chunk['embedding_vector']
                 .str.split(',')
                 .map(lambda xs: np.array(xs, dtype=np.float32))
        )
        mm[idx:idx+M] = emb_mat
        # store metadata
        id_list .extend(chunk[args.id_column].tolist())
        ws_list .extend(chunk['window_start'].tolist())
        we_list .extend(chunk['window_end'].tolist())
        sl_list .extend(chunk['seq_len'].tolist())
        idx += M
        print(f"  loaded {idx}/{N:,}")

    # flush to disk
    mm.flush()

    # 4) write mapping table (so index row i → mapping.iloc[i])
    map_df = pd.DataFrame({
        args.id_column: id_list,
        'window_start': ws_list,
        'window_end':   we_list,
        'seq_len':      sl_list
    })
    map_df.to_csv(args.mapping_path, sep='\t', index=False)
    print(f"Wrote mapping to {args.mapping_path}")

    # 5) build the FAISS index *incrementally*, chunk by chunk
    print("Building FAISS index in chunks …")
    index = faiss.IndexFlatL2(D)
    reader2 = np.memmap(memmap_path, dtype='float32', mode='r', shape=(N, D))
    for i in range(0, N, 10_000):
        j = min(N, i+10_000)
        index.add(reader2[i:j])
        print(f"  added rows {i}-{j-1}")

    # 6) save
    faiss.write_index(index, args.index_path)
    print(f"FAISS index saved to {args.index_path}")

if __name__ == '__main__':
    main()
