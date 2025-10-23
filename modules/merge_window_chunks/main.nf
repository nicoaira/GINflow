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
    tuple val(chunk_ids),
          path(db_windows_list, stageAs: 'db_windows_*.npy'),
          path(db_metadata_list, stageAs: 'db_metadata_*.tsv'),
          path(query_windows_list, stageAs: 'query_windows_*.npy'),
          path(query_metadata_list, stageAs: 'query_metadata_*.tsv'),
          path(stats_list, stageAs: 'stats_*.json')

    output:
    path 'database_windows.npy', emit: database_vectors
    path 'database_windows.tsv', emit: database_metadata
    path 'query_windows.npy',    emit: query_vectors
    path 'query_windows.tsv',    emit: query_metadata
    path 'window_stats.json',    emit: window_stats
    path 'chunk_manifest.tsv',   emit: chunk_manifest

    script:
    if (!chunk_ids) {
        error 'No window vector chunks produced; merging cannot continue'
    }
    // All files are now staged in the work directory by Nextflow
    // Build manifest using staged file names
    def manifestLines = [chunk_ids, db_windows_list, db_metadata_list, query_windows_list, query_metadata_list, stats_list]
        .transpose()
        .collect { row -> row.join('\t') }
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
