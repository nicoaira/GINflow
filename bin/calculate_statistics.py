import argparse
import numpy as np
import pandas as pd
from scipy.stats import gumbel_r

def calculate_e_value(score, loc, scale, db_size, query_length):
    p_value = 1 - gumbel_r.cdf(score, loc=loc, scale=scale)
    e_value = p_value * db_size * query_length
    return e_value

def main():
    parser = argparse.ArgumentParser(description="Calculate E-values for alignment scores.")
    parser.add_argument("--alignments_tsv", required=True, help="Path to the alignments TSV file.")
    parser.add_argument("--output_tsv", required=True, help="Path to the output TSV file with E-values.")
    parser.add_argument("--db_size", required=True, type=int, help="Total size of the database.")
    parser.add_argument("--loc", required=True, type=float, help="Location parameter of the Gumbel distribution.")
    parser.add_argument("--scale", required=True, type=float, help="Scale parameter of the Gumbel distribution.")
    args = parser.parse_args()

    alignments_df = pd.read_csv(args.alignments_tsv, sep='\t')
    
    alignments_df['e_value'] = alignments_df.apply(
        lambda row: calculate_e_value(row['score'], args.loc, args.scale, args.db_size, row['q_len']),
        axis=1
    )

    alignments_df.to_csv(args.output_tsv, sep='\t', index=False)

if __name__ == "__main__":
    main()
