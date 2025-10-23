#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_REPORT {
    tag { "query_" + alignments_tsv.baseName.replaceAll(/^alignments_/, '') }

    label 'lightweight'

    publishDir "${params.outdir}/query_results", mode: 'copy',
        saveAs: { filename ->
            def qid = filename.replaceAll(/^alignments_report_/, '').replaceAll(/\.html$/, '')
            "${qid}/${filename}"
        }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-generate-report:latest' :
        'nicoaira/ginflow-generate-report:latest' }"

    input:
    tuple path(alignments_tsv), path(alignment_individual)

    output:
    path "alignments_report_*.html"

    script:
    def query_id = alignments_tsv.baseName.replaceAll(/^alignments_/, '')
    """
    python3 ${baseDir}/bin/generate_report.py \
      --pairs ${alignments_tsv} \
      --svg-dir ${alignment_individual} \
      --id-column ${params.id_column} \
      --output alignments_report_${query_id}.html
    """
}
