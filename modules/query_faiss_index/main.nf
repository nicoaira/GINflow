#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process QUERY_FAISS_INDEX {
    tag  "query_faiss_index"

    label params.faiss_use_gpu ? 'gpu' : 'medium_memory'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'seeds.tsv'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-query-faiss-index:latest' :
        'nicoaira/ginflow-query-faiss-index:latest' }"

    input:
    path faiss_idx
    path database_vectors
    path database_metadata
    path query_vectors
    path query_metadata

    output:
    path "seeds.tsv", emit: seeds
    path "query_stats.json", optional: true, emit: query_stats

    script:
    def metric = params.faiss_metric ?: 'ip'
    def gpuFlag = params.faiss_use_gpu ? '--use-gpu' : ''
    def nprobe = params.faiss_nprobe ? "--nprobe ${params.faiss_nprobe}" : ''
    def hnswEfs = params.faiss_hnsw_efs ? "--hnsw-ef-search ${params.faiss_hnsw_efs}" : ''
    def rescoreDefault = params.index_type in ['hnsw', 'hnswsq8']
    def rescoreEffective = params.faiss_exact_rescore != null ? params.faiss_exact_rescore : rescoreDefault
    def exactRescore = rescoreEffective ? '--exact-rescore' : ''
    """
    python3 ${baseDir}/bin/query_faiss_index.py \
      --index ${faiss_idx} \
      --database-vectors ${database_vectors} \
      --database-metadata ${database_metadata} \
      --query-vectors ${query_vectors} \
      --query-metadata ${query_metadata} \
      --id-column ${params.id_column} \
      --top-k ${params.faiss_k ?: 50} \
      --similarity-threshold ${params.seed_similarity_threshold ?: 0.7} \
      --metric ${metric} \
      ${gpuFlag} \
      ${nprobe} \
      ${hnswEfs} \
      ${exactRescore} \
      --output seeds.tsv \
      --stats-json query_stats.json
    """
}
