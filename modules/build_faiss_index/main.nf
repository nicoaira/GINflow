#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BUILD_FAISS_INDEX {
    tag "build_faiss_index"

    label 'high_cpu'

    publishDir "${params.outdir}/faiss_index", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    tuple path(vectors), path(metadata)

    output:
    path "faiss.index",       emit: faiss_idx
    path "faiss_mapping.tsv", emit: faiss_map
    path "faiss_index_stats.json", optional: true, emit: faiss_stats

    script:
    def indexType = params.index_type ?: 'flat_ip'
    def metric    = params.faiss_metric ?: 'ip'
    def nlist     = params.faiss_nlist ?: 1024
    def pqM       = params.faiss_pq_m   ?: 32
    def pqBits    = params.faiss_pq_bits ?: 8
    def opqM      = params.faiss_opq_m ?: 64
    def hnswM     = params.faiss_hnsw_m ?: 32
    def hnswEfc   = params.faiss_hnsw_efc ?: 200
    def gpuFlag   = params.faiss_use_gpu ? '--use-gpu' : ''
    def nprobe    = params.faiss_nprobe ? "--nprobe ${params.faiss_nprobe}" : ''
    def addBatch  = params.faiss_add_batch_size ? "--add-batch-size ${params.faiss_add_batch_size}" : ''
    """
    python3 ${baseDir}/bin/build_faiss_index.py \
      --vectors ${vectors} \
      --metadata ${metadata} \
      --index-type ${indexType} \
      --metric ${metric} \
      --nlist ${nlist} \
      --pq-m ${pqM} \
      --pq-bits ${pqBits} \
      --opq-m ${opqM} \
      --hnsw-m ${hnswM} \
      --hnsw-efc ${hnswEfc} \
      --output-index faiss.index \
      --output-mapping faiss_mapping.tsv \
      --stats-json faiss_index_stats.json \
      ${gpuFlag} \
      ${nprobe} \
      ${addBatch}
    """
}
