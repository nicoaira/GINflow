import argparse
import json
import numpy as np
from scipy.stats import gumbel_r

def main():
    parser = argparse.ArgumentParser(description="Fit a Gumbel distribution to background scores.")
    parser.add_argument("--alignment_stats", required=True, help="Path to the alignment_stats.json file.")
    parser.add_argument("--output_json", required=True, help="Path to the output JSON file for Gumbel parameters.")
    args = parser.parse_args()

    with open(args.alignment_stats, 'r') as f:
        stats = json.load(f)
        background_scores = stats.get('background_scores')

    if not background_scores or len(background_scores) < 2:
        # Handle case where there are no/few background scores, output defaults
        loc, scale = 0.0, 1.0
    else:
        loc, scale = gumbel_r.fit(background_scores)

    gumbel_params = {'loc': loc, 'scale': scale}

    with open(args.output_json, 'w') as f:
        json.dump(gumbel_params, f)

if __name__ == "__main__":
    main()
