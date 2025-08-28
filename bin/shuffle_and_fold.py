#!/usr/bin/env python3
import argparse, random, subprocess, pandas as pd
from collections import defaultdict

def dinuc_shuffle(seq: str) -> str:
    if len(seq) < 2:
        return seq
    edges = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        edges[a].append(b)
    for k in edges:
        random.shuffle(edges[k])
    start = seq[0]
    path = [start]
    out = []
    while path:
        v = path[-1]
        if edges[v]:
            path.append(edges[v].pop())
        else:
            out.append(path.pop())
    shuffled = ''.join(out[::-1])
    if len(shuffled) != len(seq):
        return dinuc_shuffle(seq)
    return shuffled

def main():
    ap = argparse.ArgumentParser(description="Shuffle sequence dinucleotides and fold with RNAfold")
    ap.add_argument('--meta', required=True, help='id_meta TSV with sequences')
    ap.add_argument('--id-column', required=True, help='ID column name')
    ap.add_argument('--query', required=True, help='Query ID to shuffle')
    ap.add_argument('--new-id', required=True, help='ID to assign to shuffled sequence')
    ap.add_argument('--output', required=True, help='Output TSV path')
    ap.add_argument('--structure-column-name', default='secondary_structure',
                    help='Name for the output structure column')
    args = ap.parse_args()

    df = pd.read_csv(args.meta, sep='\t', dtype=str)
    row = df[df[args.id_column] == args.query]
    if row.empty:
        raise SystemExit(f"ID {args.query} not found in meta file")
    seq_col = [c for c in row.columns if 'sequence' in c.lower()]
    if not seq_col:
        raise SystemExit('No sequence column found in meta file')
    seq = row.iloc[0][seq_col[0]]
    shuffled = dinuc_shuffle(seq)

    proc = subprocess.run(['RNAfold', '--noPS'], input=(shuffled+'\n').encode(), stdout=subprocess.PIPE, check=True)
    lines = proc.stdout.decode().strip().splitlines()
    if len(lines) < 2:
        raise SystemExit('RNAfold produced unexpected output')
    structure = lines[1].split()[0]
    L = len(shuffled)

    with open(args.output, 'w') as fh:
        fh.write(
            f"{args.id_column}\tsequence\t{args.structure_column_name}\twindow_start\twindow_end\tseq_len\n"
        )
        fh.write(
            f"{args.new_id}\t{shuffled}\t{structure}\t0\t{L}\t{L}\n"
        )

if __name__ == '__main__':
    main()
