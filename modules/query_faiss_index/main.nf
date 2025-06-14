#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process QUERY_FAISS_INDEX {
    tag  "query_faiss_index"

    label 'medium_memory'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-query-faiss-index:latest' :
        'nicoaira/ginflow-query-faiss-index:latest' }"

    input:
    path embeddings
    path faiss_idx
    path faiss_map

    output:
    path "distances.tsv", emit: distances

    script:
    """
    python3 ${baseDir}/bin/query_faiss_index.py \
      --input embeddings.tsv \
      --id-column ${params.id_column} \
      --query ${params.query} \
      --index-path faiss.index \
      --mapping-path faiss_mapping.tsv \
      --top-k ${params.faiss_k} \
      --output distances.tsv
    """
}