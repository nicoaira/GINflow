#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SORT_DISTANCES {
    tag "sort_distances"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path distances

    output:
    path "distances.sorted.tsv", emit: sorted_distances

    script:
    """
    python3 - << 'PY'
import pandas as pd
df = pd.read_csv('distances.tsv', sep='\\t')
df.sort_values('distance', inplace=True)
df.to_csv('distances.sorted.tsv', sep='\\t', index=False)
PY
    """
}