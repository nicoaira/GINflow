#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SHUFFLE_AND_FOLD {
    tag { "shuffle_and_fold_${orig_id}" }

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-shuffle-and-fold:latest' :
        'nicoaira/ginflow-shuffle-and-fold:latest' }"

    input:
    tuple val(orig_id), val(new_id)
    path meta_map

    output:
    tuple val(new_id), path("${new_id}.tsv"), emit: shuffled

    script:
    """
    python3 ${baseDir}/bin/shuffle_and_fold.py \
      --meta ${meta_map} \
      --id-column ${params.id_column} \
      --structure-column-name ${params.structure_column_name} \
      --query ${orig_id} \
      --new-id ${new_id} \
      --output ${new_id}.tsv
    """
}
