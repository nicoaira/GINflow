#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
from dask import dataframe as dd


def compute_pair_score(args):
    group, alpha1, beta1, alpha2, beta2, gamma = args
    # now read the lengths we carried through
    L = int(group['len_1'].iat[0])
    M = int(group['len_2'].iat[0])
    M_mat = np.zeros((L, M))
    denom = np.zeros((L, M))
    for _, row in group.iterrows():
        s1, e1 = int(row['window_start_1']), int(row['window_end_1'])
        s2, e2 = int(row['window_start_2']), int(row['window_end_2'])
        r = float(row['rnk'])
        f_num = r**(-alpha1)
        f_den = r**(-alpha2)
        for i in range(s1, e1+1):
            for j in range(s2, e2+1):
                delta = abs((i-s1) - (j-s2))
                M_mat[i, j]   += f_num * np.exp(-beta1 * delta)
                denom[i, j]   += f_den * np.exp(-beta2 * delta)
    mask = denom > 0
    M_mat[mask] /= denom[mask]**gamma
    return int(round(M_mat.sum() / 1e5))


def find_contigs(df_grp):
    """
    Collapse a set of window pairs into connected components.
    Two windows belong to the same contig iff
      [s1,e1] overlaps [s1',e1']   AND
      [s2,e2] overlaps [s2',e2'].
    The algorithm is O(n²) but n is tiny (top-percentile windows),
    so it is fast and – unlike the old sweep – always finds the
    full transitive closure.
    """
    n = len(df_grp)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # brute-force pairwise test (n is small, so this is OK)
    s1  = df_grp['window_start_1'].values
    e1  = df_grp['window_end_1'].values
    s2  = df_grp['window_start_2'].values
    e2  = df_grp['window_end_2'].values

    for i in range(n):
        for j in range(i + 1, n):
            if (s1[i] <= e1[j] and s1[j] <= e1[i]) and \
               (s2[i] <= e2[j] and s2[j] <= e2[i]):
                union(i, j)

    # gather indices for each connected component
    comps = {}
    for i in range(n):
        r = find(i)
        comps.setdefault(r, []).append(i)

    # return a DataFrame per contig
    return [df_grp.iloc[idx_list].reset_index(drop=True)
            for idx_list in comps.values()]

def main():
    parser = argparse.ArgumentParser(
        description='Compute aggregated pair scores for windows or contigs')
    parser.add_argument('--input',             required=True,
                        help='Path to the sorted distances TSV')
    parser.add_argument('--id-column',         default='exon_id',
                        help='Name of the original ID column (without _1/_2).')
    parser.add_argument('--percentile',        type=float, default=0.01,
                        help='Percentile of top windows to use')
    parser.add_argument('--mode',              choices=['global','contigs'], default='global',
                        help='Aggregation mode')
    parser.add_argument('--num-workers',       type=int, default=1,
                        help='Workers for parallel operations')
    parser.add_argument('--alpha1',            type=float, default=1.0)
    parser.add_argument('--beta1',             type=float, default=0.1)
    parser.add_argument('--alpha2',            type=float, default=1.0)
    parser.add_argument('--beta2',             type=float, default=0.1)
    parser.add_argument('--gamma',             type=float, default=0.5)
    parser.add_argument('--output',            default='exon_pair_scores.tsv',
                        help='Output TSV path for aggregated scores')
    parser.add_argument('--output-unaggregated', action='store_true',
                        help='Also output unaggregated windows')
    args = parser.parse_args()

    # Dynamically compute ID column names
    id1_col = f"{args.id_column}_1"
    id2_col = f"{args.id_column}_2"

    # 1) load & filter only the columns we need
    usecols = [
        id1_col, 'window_start_1', 'window_end_1', 'seq_len_1',
        id2_col, 'window_start_2', 'window_end_2', 'seq_len_2',
        'distance'
    ]
    ddf = dd.read_csv(
        args.input, sep='\t',
        usecols=usecols,
        dtype={id1_col: str, id2_col: str}
    )
    total = ddf.shape[0].compute()
    top_n = max(1, int(total * args.percentile / 100.0))

    df = ddf.head(top_n, compute=True).reset_index(drop=True)

    # 2) annotate windows & lengths
    df['window_distance'] = df['distance']
    df['window_rank']     = np.arange(1, len(df) + 1)
    df['len_1']           = df['seq_len_1'].astype(int)
    df['len_2']           = df['seq_len_2'].astype(int)

    if args.output_unaggregated:
        df_unagg = df.copy()
        # ensure the 'rnk' column exists on df_unagg too
        df_unagg['rnk'] = df_unagg['window_rank']

    df = df[[
        id1_col, 'window_start_1', 'window_end_1',
        id2_col, 'window_start_2', 'window_end_2',
        'len_1', 'len_2', 'window_distance', 'window_rank'
    ]]
    df['rnk'] = df['window_rank']

    # 3) unify ID order on each row
    def unify(r):
        a, b = r[id1_col], r[id2_col]
        if a <= b:
            return r
        return pd.Series({
            id1_col:          b,
            'window_start_1': r['window_start_2'],
            'window_end_1':   r['window_end_2'],
            id2_col:          a,
            'window_start_2': r['window_start_1'],
            'window_end_2':   r['window_end_1'],
            'len_1':          r['len_2'],
            'len_2':          r['len_1'],
            'window_distance': r['window_distance'],
            'window_rank':    r['window_rank'],
            'rnk':            r['rnk']
        })

    df = df.apply(unify, axis=1)
    if args.output_unaggregated:
        df_unagg = df_unagg.apply(unify, axis=1)

    # 4) group by ID-pairs and collapse into contigs
    df['pair'] = list(zip(df[id1_col], df[id2_col]))
    aggregated = []
    unagg = []
    old_id = 0

    for (e1, e2), grp in df.groupby('pair'):
        grp = grp.reset_index(drop=True)
        contigs = [grp] if args.mode == 'global' else find_contigs(grp)
        for sub in contigs:
            old_id += 1
            nwin = len(sub)
            G = compute_pair_score((
                sub, args.alpha1, args.beta1,
                args.alpha2, args.beta2, args.gamma
            ))
            if args.mode == 'global':
                aggregated.append((e1, e2, old_id, nwin, G))
            else:
                cs1, ce1 = sub.window_start_1.min(), sub.window_end_1.max()
                cs2, ce2 = sub.window_start_2.min(), sub.window_end_2.max()
                aggregated.append((e1, cs1, ce1, e2, cs2, ce2, old_id, nwin, G))

            if args.output_unaggregated:
                for _, r in sub.iterrows():
                    rec = r.to_dict()
                    rec['old_contig_id'] = old_id
                    unagg.append(rec)

    # 5) build & sort the contig-level table
    if args.mode == 'global':
        cols_out = [id1_col, id2_col, 'old_contig_id', 'n_collapsed_windows', 'score']
    else:
        cols_out = [
            id1_col, 'contig_start_1', 'contig_end_1',
            id2_col, 'contig_start_2', 'contig_end_2',
            'old_contig_id', 'n_collapsed_windows', 'score'
        ]

    df_agg = pd.DataFrame(aggregated, columns=cols_out)
    df_agg.sort_values('score', ascending=False, inplace=True)
    df_agg['contig_id']   = np.arange(1, len(df_agg) + 1)
    df_agg['contig_rank'] = df_agg['contig_id']

    if args.mode == 'global':
        final_cols = [id1_col, id2_col, 'contig_id', 'contig_rank', 'n_collapsed_windows', 'score']
    else:
        final_cols = [
            id1_col, 'contig_start_1', 'contig_end_1',
            id2_col, 'contig_start_2', 'contig_end_2',
            'contig_id', 'contig_rank', 'n_collapsed_windows', 'score'
        ]

    df_agg[final_cols].to_csv(args.output, sep='\t', index=False)
    print(f"Written aggregated to {args.output}")

    # 6) write the unaggregated windows if requested
    if args.output_unaggregated:
        base, ext = os.path.splitext(args.output)
        fn = f"{base}.unaggregated{ext}"
        mapping = df_agg.set_index('old_contig_id')[[
            'contig_id', 'contig_rank', 'n_collapsed_windows', 'score'
        ]].to_dict('index')

        rows = []
        for rec in unagg:
            m = mapping[rec['old_contig_id']]
            rec.update({
                'contig_id':            m['contig_id'],
                'contig_rank':          m['contig_rank'],
                'n_collapsed_windows':  m['n_collapsed_windows'],
                'contig_score':        m['score']
            })
            rec.pop('old_contig_id')
            rows.append(rec)

        df_un = pd.DataFrame(rows)
        df_un.sort_values('contig_score', ascending=False, inplace=True)
        df_un.to_csv(fn, sep='\t', index=False)
        print(f"Written unaggregated to {fn}")

if __name__ == '__main__':
    main()