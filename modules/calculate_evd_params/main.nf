#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CALCULATE_EVD_PARAMS {
    tag "calculate_evd_params"

    label 'high_cpu'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'evd_params.json'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-calculate-evd-params:latest' :
        'docker.io/nicoaira/ginflow-calculate-evd-params:latest' }"

    input:
    path embeddings

    output:
    path 'evd_params.json', emit: evd_params

    script:
    def evalueArgs = (params.calculate_evalue == false) ? "" : ""
    """
    python3 ${baseDir}/bin/calculate_evd_params.py \
        --embeddings ${embeddings} \
        --id-column ${params.id_column} \
        --node-index-column node_index \
        --embedding-column embedding_vector \
        --background-samples ${params.background_samples ?: 10000} \
        --evd-samples ${params.evd_samples ?: 1000} \
        --gamma ${params.alignment_gamma ?: 1.5} \
        --band-width ${params.band_width ?: 96} \
        --gap-open ${params.gap_open ?: 12} \
        --gap-extend ${params.gap_extend ?: 2} \
        --xdrop ${params.xdrop ?: 50} \
        --score-min ${params.score_min ?: -4} \
        --score-max ${params.score_max ?: 8} \
        --random-seed ${params.random_seed ?: 42} \
        --output evd_params.json
    """
}
