#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_UNAGGREGATED_REPORT {
    tag "gen_unagg_report"

    label 'lightweight'
    
    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-generate-report:latest' :
        'nicoaira/ginflow-generate-report:latest' }"

    input:
    path top_contigs_unagg_tsv
    path window_individual

    output:
    path "pairs_contigs_report.unaggregated.html"

    script:
    """
    python3 ${baseDir}/bin/generate_report.py \
      --pairs ${top_contigs_unagg_tsv} \
      --svg-dir ${window_individual} \
      --id-column ${params.id_column} \
      --output pairs_contigs_report.unaggregated.html
    """
}