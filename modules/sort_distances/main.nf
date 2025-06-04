#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SORT_DISTANCES {
    tag "sort_distances"
    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-slim-2.2.2:latest' :
        'amancevice/pandas:slim-2.2.2' }"        

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