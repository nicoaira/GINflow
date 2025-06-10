#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AGGREGATE_SCORE {
    tag "aggregate_score"

    label 'lightweight'

    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-aggregate-score:latest' :
        'nicoaira/ginflow-aggregate-score:latest' }"

    input:
    path sorted_distances
    path meta_map

    output:
    path "pairs_scores_all_contigs.tsv",              emit: enriched_all
    path "pairs_scores_all_contigs.unaggregated.tsv", emit: enriched_unagg

    script:
    """
    python3 ${baseDir}/bin/aggregated_score.py \
      --input distances.sorted.tsv \
      --id-column ${params.id_column} \
      --alpha1 ${params.alpha1} --alpha2 ${params.alpha2} \
      --beta1 ${params.beta1}   --beta2 ${params.beta2} \
      --gamma ${params.gamma}   --percentile ${params.percentile} \
      --mode contigs \
      --output raw_contigs.tsv \
      --output-unaggregated

    python3 - << 'PY'
import pandas as pd
idc      = '${params.id_column}'

# meta_map has IDs + *all* sequence and structure columns
meta     = pd.read_csv('id_meta.tsv', sep='\\t', dtype=str)
meta_cols = [c for c in meta.columns if c != idc]

m1 = meta.rename(columns={idc:f'{idc}_1', **{c:f'{c}_1' for c in meta_cols}})
m2 = meta.rename(columns={idc:f'{idc}_2', **{c:f'{c}_2' for c in meta_cols}})

all_df   = pd.read_csv('raw_contigs.tsv', sep='\\t').merge(m1,on=f'{idc}_1').merge(m2,on=f'{idc}_2')
unagg_df = pd.read_csv('raw_contigs.unaggregated.tsv', sep='\\t').merge(m1,on=f'{idc}_1').merge(m2,on=f'{idc}_2')

all_df.to_csv('pairs_scores_all_contigs.tsv', sep='\\t', index=False)
unagg_df.to_csv('pairs_scores_all_contigs.unaggregated.tsv', sep='\\t', index=False)
PY
    """
}