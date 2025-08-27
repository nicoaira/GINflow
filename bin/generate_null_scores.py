#!/usr/bin/env python3
import argparse
import pandas as pd
import random
import subprocess


def dinuc_shuffle(seq: str) -> str:
    pairs = [seq[i:i+2] for i in range(0, len(seq)-1, 2)]
    random.shuffle(pairs)
    tail = seq[-1] if len(seq) % 2 else ''
    return ''.join(pairs) + tail


def fold_energy(seq: str) -> float:
    proc = subprocess.run(
        ['RNAfold', '--noPS'],
        input=f"{seq}\n", text=True, capture_output=True
    )
    out = proc.stdout.strip().splitlines()
    if len(out) < 2:
        return None
    line = out[1]
    try:
        energy = line.split('(')[-1].split(')')[0]
        return float(energy)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description='Generate null score distribution by dinucleotide shuffling and folding')
    ap.add_argument('--queries', required=True, help='CSV with query ids')
    ap.add_argument('--meta-map', required=True, help='TSV with id and sequence columns')
    ap.add_argument('--iterations', type=int, default=100, help='Shuffles per query')
    ap.add_argument('--output', required=True, help='Output TSV with score column')
    args = ap.parse_args()

    meta = pd.read_csv(args.meta_map, sep='\t', dtype=str)
    id_col = meta.columns[0]
    seq_col = next((c for c in meta.columns if 'sequence' in c.lower()), None)
    if seq_col is None:
        raise SystemExit('No sequence column found in meta map')

    queries = pd.read_csv(args.queries)['id']
    scores = []
    for qid in queries:
        row = meta[meta[id_col] == qid]
        if row.empty:
            continue
        seq = row.iloc[0][seq_col]
        for _ in range(args.iterations):
            shuf = dinuc_shuffle(seq)
            energy = fold_energy(shuf)
            if energy is not None:
                scores.append({'score': energy})

    pd.DataFrame(scores).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
