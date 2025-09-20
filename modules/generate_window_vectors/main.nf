#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_WINDOW_VECTORS {
    tag { "window_vec_${chunk_id}" }

    label 'high_cpu'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(chunk_id), path(node_embeddings), path(queries)

    output:
    tuple val(chunk_id), path('database_windows.npy'), path('database_windows.tsv'), path('query_windows.npy'), path('query_windows.tsv'), path('window_stats.json'), emit: window_chunk

    script:
    def stride = params.window_stride ?: 1
    def output_dir = 'windows'
    def queryColArg = params.query_column ? "--query-column ${params.query_column}" : null
    def chunkArg = chunk_id ? "--chunk-id ${chunk_id}" : null
    def cliArgs = [
        "--embeddings ${node_embeddings}",
        "--queries ${queries}",
        "--id-column ${params.id_column}",
        "--window-size ${params.window_size ?: 11}",
        "--stride ${stride}",
        "--output-dir ${output_dir}",
        "--metadata-json window_stats.json"
    ] + [queryColArg, chunkArg].findAll { it }
    def cliBlock = cliArgs.collect { "        ${it}" }.join(' \\\n')
    """
    mkdir -p ${output_dir}
    python3 ${baseDir}/bin/generate_window_vectors.py \
${cliBlock}

    mv ${output_dir}/database_windows.npy ./database_windows.npy
    mv ${output_dir}/database_windows.tsv ./database_windows.tsv
    mv ${output_dir}/query_windows.npy ./query_windows.npy
    mv ${output_dir}/query_windows.tsv ./query_windows.tsv
    """
}
