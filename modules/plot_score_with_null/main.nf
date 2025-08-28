#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PLOT_SCORE_WITH_NULL {
    tag { "plot_score_with_null_${query_id}" }

    label 'lightweight'

    when   { params.plot_score_distribution }
    publishDir "${params.outdir}/queries_results/${query_id}/plots", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-plot-score:latest' :
        'nicoaira/ginflow-plot-score:latest' }"

    input:
    tuple val(query_id), path(enriched_all), path(null_scores)

    output:
    path "score_distribution.png"

    script:
    """
    python3 - << 'PY'
import pandas as pd
import matplotlib.pyplot as plt

real_df = pd.read_csv('pairs_scores_all_contigs.tsv', sep='\t')
null_df = pd.read_csv('null_scores.tsv', sep='\t')

plt.figure(figsize=(6,4))
bins = ${params.score_bins}

# Plot real distribution
plt.hist(real_df['score'], bins=bins, alpha=0.6, label='Real', color='#1f77b4', density=True)

# Plot null distribution overlay
if 'score' in null_df.columns and len(null_df) > 0:
    plt.hist(null_df['score'], bins=bins, alpha=0.6, label='Null', color='#ff7f0e', density=True)

plt.xlabel('Score')
plt.ylabel('Density')
plt.title('Contig Score vs Null Distribution')
plt.legend()
plt.tight_layout()
plt.savefig('score_distribution.png', dpi=150)
PY
    """
}

