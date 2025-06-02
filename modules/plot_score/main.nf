#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PLOT_SCORE {
    when   { params.plot_score_distribution }
    tag    "plot_score"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path enriched_all

    output:
    path "score_distribution.png"

    script:
    """
    python3 - << 'PY'
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv('pairs_scores_all_contigs.tsv', sep='\\t')
plt.hist(df['score'], bins=${params.score_bins})
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.title('Contig Score Distribution')
plt.tight_layout()
plt.savefig('score_distribution.png')
PY
    """
}