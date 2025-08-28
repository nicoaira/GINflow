#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process QUERY_FAISS_INDEX_BY_TUPLE {
    tag  { "query_faiss_index_${query_id}" }

    label 'medium_memory'

    // No publishDir here; null results are not copied to outdir

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-query-faiss-index:latest' :
        'nicoaira/ginflow-query-faiss-index:latest' }"

    input:
    tuple val(query_id), path(embeddings)
    path faiss_idx
    path faiss_map

    output:
    tuple val(query_id), path("distances.tsv"), emit: distances

    script:
    """
    python3 ${baseDir}/bin/query_faiss_index.py \
      --input ${embeddings} \
      --id-column ${params.id_column} \
      --query ${query_id} \
      --index-path ${faiss_idx} \
      --mapping-path ${faiss_map} \
      --top-k ${params.faiss_k} \
      --output distances.tsv
    """
}
