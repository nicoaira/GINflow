#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SPLIT_QUERIES {
    tag "split_queries"

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-split-queries:latest' :
        'docker.io/nicoaira/ginflow-split-queries:latest' }"

    input:
    path query_vectors
    path query_metadata

    output:
    path 'queries/*_vectors.npy', emit: query_vectors
    path 'queries/*_metadata.tsv', emit: query_metadata
    path 'queries/query_split_manifest.json', emit: manifest

    script:
    """
    python3 ${baseDir}/bin/split_queries.py \
        --query-vectors ${query_vectors} \
        --query-metadata ${query_metadata} \
        --id-column ${params.id_column} \
        --output-dir queries
    """
}
