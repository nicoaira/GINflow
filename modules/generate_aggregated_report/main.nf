#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_AGGREGATED_REPORT {
    tag "gen_agg_report"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path top_contigs_tsv
    path contig_individual

    output:
    path "pairs_contigs_report.html"

    script:
    """
    python3 ${baseDir}/bin/generate_report.py \
      --pairs ${top_contigs_tsv} \
      --svg-dir ${contig_individual} \
      --id-column ${params.id_column} \
      --output pairs_contigs_report.html
    """
}