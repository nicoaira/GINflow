#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_EMBEDDINGS {
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ginflow:0.1.0' :
        'ginflow:0.1.0' }"

    tag "merge_embeddings"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path batch_embeddings

    output:
    path "embeddings.tsv", emit: embeddings

    script:
    """
    head -n1 ${batch_embeddings[0]} > embeddings.tsv
    for f in ${batch_embeddings.join(' ')}; do
      tail -n +2 \$f >> embeddings.tsv
    done
    """
}