#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DRAW_CONTIG_SVGS {
    tag "draw_contig_svgs"
    cpus 16
    publishDir "${params.outdir}/drawings/contigs", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    // TODO: update singularity container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-draw-svgs:latest' :
        'docker.io/nicoaira/ginflow-draw-structures:latest' }"

    input:
    path top_contigs_tsv

    output:
    path 'individual_svgs', emit: contig_individual

    script:
    """
    mkdir -p individual_svgs
    python3 ${baseDir}/bin/draw_structures.py \
      --tsv ${top_contigs_tsv} --outdir individual_svgs \
      --pair-type "contig" \
      --id-column ${params.id_column} \
      --width 500 --height 500 --highlight-colour "#00FF99" \
      --num-workers ${task.cpus}
    """
}