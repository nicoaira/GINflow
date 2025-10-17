#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DRAW_UNAGG_SVGS_R4RNA {
    tag "draw_window_svgs_r4rna_${query_id}"

    label 'high_cpu'

    publishDir "${params.outdir}/queries_results/${query_id}/drawings/unagg_windows", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-draw-structures-r4rna:latest' :
        'docker.io/nicoaira/ginflow-draw-structures-r4rna:latest' }"

    input:
    tuple val(query_id), path(top_windows_tsv)

    output:
    tuple val(query_id), path('individual_svgs'), emit: window_individual

    script:
    """
    mkdir -p individual_svgs
    python3 ${baseDir}/bin/draw_structures_r4rna.py \
      --tsv ${top_windows_tsv} \
      --outdir individual_svgs \
      --pair-type "window" \
      --id-column ${params.id_column} \
      --highlight-colour "#FF0000" \
      --num-workers ${task.cpus}
    """
}
