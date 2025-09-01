#!/usr/bin/env python3
import argparse
import pandas as pd


def parse_pairs(structure: str) -> dict:
    stack = []
    pairs = {}
    for i, ch in enumerate(structure, start=1):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
    return pairs


def has_cross_pair(pairs: dict, s1: int, e1: int, s2: int, e2: int) -> bool:
    for i in range(s1, e1 + 1):
        j = pairs.get(i)
        if j and s2 <= j <= e2:
            return True
    for i in range(s2, e2 + 1):
        j = pairs.get(i)
        if j and s1 <= j <= e1:
            return True
    return False


def overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def _detect_structure_col(cols, prefix: str) -> str:
    """Return the structure column name for a given prefix ('query'|'subject').
    Supports both '*_secondary_structure' and '*_structure', and falls back to
    any column containing 'structure' with the given prefix.
    """
    candidates = [f"{prefix}_secondary_structure", f"{prefix}_structure"]
    for c in candidates:
        if c in cols:
            return c
    # Fallback: any prefixed column that contains 'structure'
    for c in cols:
        if c.startswith(f"{prefix}_") and 'structure' in c:
            return c
    raise KeyError(f"No structure column found for prefix '{prefix}' (expected one of: {', '.join(candidates)})")


def merge_group(df_grp: pd.DataFrame, max_ov: int, q_colname: str = None, s_colname: str = None) -> list:
    # Resolve structure columns deterministically if provided; else detect
    if q_colname and s_colname:
        q_col = q_colname
        s_col = s_colname
    else:
        # Detect structure columns (query/subject) robustly across schema variants
        q_col = _detect_structure_col(df_grp.columns, 'query')
        s_col = _detect_structure_col(df_grp.columns, 'subject')
    q_struct = df_grp[q_col].iloc[0]
    s_struct = df_grp[s_col].iloc[0]
    q_pairs = parse_pairs(q_struct)
    s_pairs = parse_pairs(s_struct)

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

    qs = df_grp['query_contig_start'].astype(int).values
    qe = df_grp['query_contig_end'].astype(int).values
    ss = df_grp['subject_contig_start'].astype(int).values
    se = df_grp['subject_contig_end'].astype(int).values

    for i in range(n):
        for j in range(i + 1, n):
            if overlap(qs[i], qe[i], qs[j], qe[j]) > max_ov:
                continue
            if overlap(ss[i], se[i], ss[j], se[j]) > max_ov:
                continue
            if not has_cross_pair(q_pairs, qs[i], qe[i], qs[j], qe[j]):
                continue
            if not has_cross_pair(s_pairs, ss[i], se[i], ss[j], se[j]):
                continue
            union(i, j)

    comps = {}
    for i in range(n):
        r = find(i)
        comps.setdefault(r, []).append(i)
    return [df_grp.iloc[idxs].reset_index(drop=True) for idxs in comps.values()]


def main():
    ap = argparse.ArgumentParser(description='Merge contigs by secondary structure connectivity')
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--max-contig-overlap', type=int, default=10)
    ap.add_argument('--structure-column-name', dest='scol', default=None,
                    help="Base name of the structure column in meta map (e.g. 'secondary_structure' or 'rnafold_dotbracket'). "
                         "When provided, expects columns 'query_<name>' and 'subject_<name>'.")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    meta_cols = [c for c in df.columns if c not in {
        'query_id', 'query_contig_start', 'query_contig_end',
        'subject_id', 'subject_contig_start', 'subject_contig_end',
        'contig_id', 'contig_rank', 'n_collapsed_windows', 'score'
    }]

    # pre-resolve structure column names if provided
    q_scol = s_scol = None
    if args.scol:
        q_scol = f"query_{args.scol}"
        s_scol = f"subject_{args.scol}"
        # Validate presence early for clearer errors
        missing = [c for c in (q_scol, s_scol) if c not in df.columns]
        if missing:
            raise KeyError(f"Expected structure columns not found: {', '.join(missing)}")

    rows = []
    for (q, s), grp in df.groupby(['query_id', 'subject_id']):
        merged = merge_group(grp, args.max_contig_overlap, q_scol, s_scol)
        for sub in merged:
            row = {
                'query_id': q,
                'subject_id': s,
                'query_contig_start': ','.join(map(str, sub['query_contig_start'].tolist())),
                'query_contig_end': ','.join(map(str, sub['query_contig_end'].tolist())),
                'subject_contig_start': ','.join(map(str, sub['subject_contig_start'].tolist())),
                'subject_contig_end': ','.join(map(str, sub['subject_contig_end'].tolist())),
                'contig_ids': ','.join(map(str, sub['contig_id'].tolist())),
                'contig_rank': ','.join(map(str, sub['contig_rank'].tolist())),
                # Preserve per-contig window counts as a comma-separated list
                'n_collapsed_windows': ','.join(map(str, sub['n_collapsed_windows'].tolist())),
                'score': int(sub['score'].sum()),
            }
            for mc in meta_cols:
                row[mc] = sub[mc].iloc[0]
            rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df.sort_values('score', ascending=False, inplace=True)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Written aggregated to {args.output}")

if __name__ == '__main__':
    main()
