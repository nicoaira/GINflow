#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FILTER_TOP_CONTIGS {
    tag "filter_top_contigs_${query_id}"

    label 'lightweight'

    publishDir "${params.outdir}/queries_results/${query_id}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-2.2.2:latest' :
        'amancevice/pandas:2.2.2' }"

    input:
    tuple val(query_id), path(enriched_agg), path(windows)

    output:
    tuple val(query_id), path("pairs_scores_top_contigs.tsv"),          emit: top_contigs
    tuple val(query_id), path("pairs_scores_top_contigs.windows.tsv"),  emit: top_contigs_unagg

    script:
    """
    python3 - << 'PY'
import pandas as pd
top = ${params.top_n}
all_df = pd.read_csv('pairs_scores_all_contigs.aggregated.tsv', sep='\t')
un_df  = pd.read_csv('pairs_scores_all_contigs.windows.tsv', sep='\t')
all_df = all_df.sort_values('score', ascending=False)
sel    = all_df.head(top)
ids = set()
for ids_str in sel['contig_ids']:
    ids.update(map(int, str(ids_str).split(',')))
sel.to_csv('pairs_scores_top_contigs.tsv', sep='\t', index=False)
un_df[un_df.contig_id.isin(ids)].to_csv('pairs_scores_top_contigs.windows.tsv', sep='\t', index=False)
PY
    """
}
