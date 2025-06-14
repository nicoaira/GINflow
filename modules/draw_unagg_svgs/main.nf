#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DRAW_UNAGG_SVGS {
    tag "draw_window_svgs"

    label 'high_cpu'
    
    publishDir "${params.outdir}/drawings/unagg_windows", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-draw-structures:latest' :
        'docker.io/nicoaira/ginflow-draw-structures:latest' }"

    input:
    path top_unagg_tsv

    output:
    path 'individual_svgs', emit: window_individual

    script:
    """
    mkdir -p individual_svgs
    python3 ${baseDir}/bin/draw_structures.py \
      --tsv ${top_unagg_tsv} --outdir individual_svgs \
      --pair-type "window" \
      --id-column ${params.id_column} \
      --width 500 --height 500 --highlight-colour "#00FF99" \
      --num-workers ${task.cpus}
    """
}