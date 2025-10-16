#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_REPORT {
    tag "generate_report"

    label 'lightweight'

    publishDir "${params.outdir}/reports", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-generate-report:latest' :
        'nicoaira/ginflow-generate-report:latest' }"

    input:
    tuple path(alignments_tsv), path(alignment_individual)

    output:
    path "alignments_report.html"

    script:
    """
    python3 ${baseDir}/bin/generate_report.py \
      --pairs ${alignments_tsv} \
      --svg-dir ${alignment_individual} \
      --id-column ${params.id_column} \
      --output alignments_report.html
    """
}
