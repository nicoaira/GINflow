#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PLOT_SCORE {
    tag "plot_score_${query_id}"

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-plot-score:latest' :
        'nicoaira/ginflow-plot-score:latest' }"

    when   { params.plot_score_distribution }
    publishDir "${params.outdir}/queries_results/${query_id}/plots", mode: 'copy'

    input:
    tuple val(query_id), path(enriched_agg)

    output:
    path "score_distribution.png"

    script:
    """
    python3 - << 'PY'
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv('pairs_scores_all_contigs.aggregated.tsv', sep='\t')
plt.hist(df['score'], bins=${params.score_bins})
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.title('Contig Score Distribution')
plt.tight_layout()
plt.savefig('score_distribution.png')
PY
    """
}
