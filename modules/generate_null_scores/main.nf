#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_NULL_SCORES {
    tag "generate_null_scores"

    label 'lightweight'

    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'oras://quay.io/nicoaira/viennarna:latest' : 'viennarna/viennarna:2.6.4' }"

    input:
    path queries
    path meta_map

    output:
    path "null_scores.tsv", emit: null_scores

    script:
    """
    python3 ${baseDir}/bin/generate_null_scores.py \
        --queries $queries \
        --meta-map $meta_map \
        --iterations ${params.null_iterations} \
        --output null_scores.tsv
    """
}
