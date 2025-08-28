#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SORT_DISTANCES_NULL {
    tag { "sort_distances_null_${query_id}" }

    label 'high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-2.2.2:latest' :
        'amancevice/pandas:2.2.2' }"

    input:
    tuple val(query_id), path(distances)

    output:
    tuple val(query_id), path("distances.sorted.tsv"), emit: sorted_distances

    script:
    """
    python3 - << 'PY'
import pandas as pd
df = pd.read_csv('distances.tsv', sep='\t')
df.sort_values('distance', inplace=True)
df.to_csv('distances.sorted.tsv', sep='\t', index=False)
PY
    """
}

