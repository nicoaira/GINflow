#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_UNAGGREGATED_REPORT {
    tag "gen_unagg_report_${query_id}"

    label 'lightweight'

    publishDir "${params.outdir}/queries_results/${query_id}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-generate-report:latest' :
        'nicoaira/ginflow-generate-report:latest' }"

    input:
    tuple val(query_id), path(top_windows_tsv), path(window_individual)

    output:
    path "pairs_contigs_report.unaggregated.html"

    script:
    """
    python3 ${baseDir}/bin/generate_report.py \
      --pairs ${top_windows_tsv} \
      --svg-dir ${window_individual} \
      --id-column ${params.id_column} \
      --output pairs_contigs_report.unaggregated.html
    """
}
