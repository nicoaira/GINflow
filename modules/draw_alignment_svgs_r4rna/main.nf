#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DRAW_ALIGNMENT_SVGS_R4RNA {
    tag "draw_alignment_svgs_r4rna"

    label 'high_cpu'

    publishDir "${params.outdir}/drawings/alignments", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-draw-structures-r4rna:latest' :
        'docker.io/nicoaira/ginflow-draw-structures-r4rna:latest' }"

    input:
    path alignments_tsv

    output:
    path 'individual_svgs', emit: alignment_individual

    script:
    """
    mkdir -p individual_svgs
    python3 ${baseDir}/bin/draw_structures_r4rna.py \
      --tsv ${alignments_tsv} \
      --outdir individual_svgs \
      --pair-type "alignment" \
      --id-column ${params.id_column} \
      --highlight-colour "#FF0000" \
      --num-workers ${task.cpus}
    """
}
