#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_EMBEDDINGS {
  
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/official-debian-bullseye-slim:latest' :
        'ubuntu:22.04' }"

    tag "merge_embeddings"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path batch_embeddings, stageAs: 'batch_*.tsv'

    output:
    path "embeddings.tsv", emit: embeddings

    script:
    """
    # Get the first file to extract header
    first_file=\$(ls batch_*.tsv | head -n 1)
    head -n 1 "\$first_file" > embeddings.tsv
    
    # Append all data rows (skip header from each file)
    for f in batch_*.tsv; do
      tail -n +2 "\$f" >> embeddings.tsv
    done
    """
}