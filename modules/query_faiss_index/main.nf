#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process QUERY_FAISS_INDEX {
    tag  "query_faiss_index_${query_id}"

    label params.use_gpu ? 'gpu' : 'medium_memory'
    maxForks = 1

    publishDir "${params.outdir}/queries_results/${query_id}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-query-faiss-index:latest' :
        'nicoaira/ginflow-query-faiss-index:latest' }"

    input:
    path embeddings
    path faiss_idx
    path faiss_map
    val query_id

    output:
    tuple val(query_id), path("distances.tsv"), emit: distances

    script:
    def gpuFlag = params.use_gpu ? '--use-gpu' : ''
    """
    python3 ${baseDir}/bin/query_faiss_index.py \
      --input ${embeddings} \
      --id-column ${params.id_column} \
      --query ${query_id} \
      --index-path ${faiss_idx} \
      --mapping-path ${faiss_map} \
      --top-k ${params.faiss_k} \
      --output distances.tsv ${gpuFlag}
    """
}
