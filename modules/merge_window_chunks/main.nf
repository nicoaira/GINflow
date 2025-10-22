#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_WINDOW_CHUNKS {
    tag 'merge_window_chunks'

    label 'high_cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-merge-window-chunks:latest' :
        'docker.io/nicoaira/ginflow-merge-window-chunks:latest' }"

    input:
    val chunk_records

    output:
    path 'database_windows.npy', emit: database_vectors
    path 'database_windows.tsv', emit: database_metadata
    path 'query_windows.npy',    emit: query_vectors
    path 'query_windows.tsv',    emit: query_metadata
    path 'window_stats.json',    emit: window_stats
    path 'chunk_manifest.tsv',   emit: chunk_manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!chunk_records) {
        error 'No window vector chunks produced; merging cannot continue'
    }
    def manifestLines = chunk_records.collect { item ->
        [item[0], item[1], item[2], item[3], item[4], item[5]].collect { it.toString() }.join('\t')
    }
    def manifestContent = (
        ['chunk_id\tdatabase_windows\tdatabase_metadata\tquery_windows\tquery_metadata\twindow_stats'] +
        manifestLines
    ).join('\n')
    """
    cat <<'EOF' > window_chunks.tsv
${manifestContent}
EOF

    python3 ${baseDir}/bin/merge_window_chunks.py \
        --manifest window_chunks.tsv \
        --database-vectors database_windows.npy \
        --database-metadata database_windows.tsv \
        --query-vectors query_windows.npy \
        --query-metadata query_windows.tsv \
        --stats window_stats.json \
        --chunk-manifest chunk_manifest.tsv

    rm -f window_chunks.tsv
    """
}
