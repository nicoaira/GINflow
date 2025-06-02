#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BUILD_FAISS_INDEX {
    tag    "build_faiss_index"
    cpus   1
    memory '16 GB'
    publishDir "${params.outdir}/faiss_index", mode: 'copy'

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