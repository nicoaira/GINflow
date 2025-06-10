#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BUILD_FAISS_INDEX {
    tag "build_faiss_index"

    label 'mediumweight'

    publishDir "${params.outdir}/faiss_index", mode: 'copy'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-build-faiss-index:latest' :
        'nicoaira/ginflow-build-faiss-index:latest' }"

    input:
    path embeddings

    output:
    path "faiss.index",       emit: faiss_idx
    path "faiss_mapping.tsv", emit: faiss_map

    script:
    """
    python3 ${baseDir}/bin/build_faiss_index.py \
      --input embeddings.tsv \
      --id-column ${params.id_column} \
      --query ${params.query} \
      --index-path faiss.index \
      --mapping-path faiss_mapping.tsv
    """
}