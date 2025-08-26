#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DRAW_CONTIG_SVGS {
    tag "draw_contig_svgs_${query_id}"

    label 'high_cpu'

    publishDir "${params.outdir}/queries_results/${query_id}/drawings/contigs", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-draw-structures:latest' :
        'docker.io/nicoaira/ginflow-draw-structures:latest' }"

    input:
    path top_contigs_tsv
    val query_id

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
