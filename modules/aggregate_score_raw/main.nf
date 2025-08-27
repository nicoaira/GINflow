#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AGGREGATE_SCORE_RAW {
    tag { "aggregate_score_raw_${query_id}" }

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-aggregate-score:latest' :
        'nicoaira/ginflow-aggregate-score:latest' }"

    input:
    tuple val(query_id), path(sorted_distances)

    output:
    tuple val(query_id), path('scores.tsv'), emit: scores

    script:
    """
    python3 ${baseDir}/bin/aggregated_score.py \
      --input ${sorted_distances} \
      --id-column ${params.id_column} \
      --alpha1 ${params.alpha1} --alpha2 ${params.alpha2} \
      --beta1 ${params.beta1}   --beta2 ${params.beta2} \
      --gamma ${params.gamma}   --percentile ${params.percentile} \
      --mode contigs \
      --output scores.tsv
    """
}
