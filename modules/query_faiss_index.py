#!/usr/bin/env python3
import argparse
import pandas as pd
import faiss
import numpy as np
import sys

def main():
    p = argparse.ArgumentParser(
        description="Query a FAISS index to get top-K neighbors, streaming only the query rows"
    )
    p.add_argument('--input',        required=True,
                   help="Same embeddings.tsv used to build the index")
    p.add_argument('--id-column',    required=True,
                   help="Name of the transcript/window ID column")
    p.add_argument('--query',        required=True,
                   help="ID value to treat as “query”")
    p.add_argument('--index-path',   required=True,
                   help="Path to the faiss index file")
    p.add_argument('--mapping-path', required=True,
                   help="DB metadata file (from build_faiss_index.py)")
    p.add_argument('--top-k',        type=int, default=100,
                   help="How many nearest neighbors per query window")
    p.add_argument('--output',       required=True,
                   help="Path to write distances.tsv")
    args = p.parse_args()

    # ─── 1) Read mapping (DB side) ───
    db_meta = pd.read_csv(args.mapping_path, sep='\t')
    db_records = db_meta.to_dict('records')

    # ─── 2) Gather only the query’s windows ───
    q_meta = []
    q_vecs_list = []
    with open(args.input, 'r') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            id_i    = header.index(args.id_column)
            start_i = header.index('window_start')
            end_i   = header.index('window_end')
            len_i   = header.index('seq_len')
            vec_i   = header.index('embedding_vector')
        except ValueError as e:
            sys.exit(f"ERROR: missing column in {args.input}: {e}")

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if parts[id_i] != args.query:
                continue
            q_meta.append({
                args.id_column: parts[id_i],
                'window_start': int(parts[start_i]),
                'window_end':   int(parts[end_i]),
                'seq_len':      int(parts[len_i])
            })
            q_vecs_list.append(np.fromstring(parts[vec_i], sep=',', dtype='float32'))

    if not q_vecs_list:
        sys.exit(f"ERROR: no rows found where {args.id_column} == {args.query}")

    q_vecs = np.stack(q_vecs_list)

    # ─── 3) Load FAISS index & over-fetch ───
    index = faiss.read_index(args.index_path)
    extra = len(q_meta)  # allow dropping up to one “self” per window
    D, I = index.search(q_vecs, args.top_k + extra)      # ← ask for a few extra

    # ─── 4) Build output rows, exactly top_k per query window ───
    rows = []
    for qi, (dists, neighs) in enumerate(zip(D, I)):
        qrec = q_meta[qi]
        kept = 0
        for dist, dbi in zip(dists, neighs):
            nrec = db_records[dbi]
            if nrec[args.id_column] == args.query:
                continue                                   # ← skip self
            rows.append({
                f"{args.id_column}_1": qrec[args.id_column],
                "window_start_1":      qrec["window_start"],
                "window_end_1":        qrec["window_end"],
                "seq_len_1":           qrec["seq_len"],

                f"{args.id_column}_2": nrec[args.id_column],
                "window_start_2":      nrec["window_start"],
                "window_end_2":        nrec["window_end"],
                "seq_len_2":           nrec["seq_len"],

                "distance":            float(dist)
            })
            kept += 1
            if kept >= args.top_k:
                break                                      # ← stop at exactly K

    # ─── 5) Write distances.tsv ───
    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Written {len(rows)} distances to {args.output}")

if __name__ == '__main__':
    main()
