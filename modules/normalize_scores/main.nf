#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process NORMALIZE_SCORES {
    tag "normalize_scores"

    label 'lightweight'

    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-2.2.2:latest' :
        'amancevice/pandas:2.2.2' }"

    input:
    path scores
    path null_dist

    output:
    path "pairs_scores_all_contigs.normalized.tsv"

    script:
    """
    python3 ${baseDir}/bin/normalize_scores.py \
      --scores $scores \
      --null-distribution $null_dist \
      --output pairs_scores_all_contigs.normalized.tsv
    """
}
