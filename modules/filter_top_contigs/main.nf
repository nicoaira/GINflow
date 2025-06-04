#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FILTER_TOP_CONTIGS {
    tag "filter_top_contigs"
    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-slim-2.2.2:latest' :
        'amancevice/pandas:slim-2.2.2' }"        

    input:
    path enriched_all
    path enriched_unagg

    output:
    path "pairs_scores_top_contigs.tsv",              emit: top_contigs
    path "pairs_scores_top_contigs.unaggregated.tsv", emit: top_contigs_unagg

    script:
    """
    python3 - << 'PY'
import pandas as pd
top = ${params.top_n}
all_df = pd.read_csv('pairs_scores_all_contigs.tsv', sep='\\t')
un_df  = pd.read_csv('pairs_scores_all_contigs.unaggregated.tsv', sep='\\t')
sel    = all_df[all_df.contig_rank <= top]
ids    = sel.contig_id.unique()
sel.to_csv('pairs_scores_top_contigs.tsv', sep='\\t', index=False)
un_df[un_df.contig_id.isin(ids)].to_csv('pairs_scores_top_contigs.unaggregated.tsv', sep='\\t', index=False)
PY
    """
}