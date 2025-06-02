#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_EMBEDDINGS {
    tag "merge_embeddings"
    container ''
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