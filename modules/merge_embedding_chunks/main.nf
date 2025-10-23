#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_EMBEDDING_CHUNKS {
    tag 'merge_embedding_chunks'

    label 'high_cpu'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-merge-embedding-chunks:latest' :
        'docker.io/nicoaira/ginflow-merge-embedding-chunks:latest' }"

    input:
    tuple val(batch_ids), path(embedding_files, stageAs: 'embedding_chunk_*.tsv')

    output:
    path 'node_embeddings.tsv', emit: node_embeddings
    path 'embedding_chunks.tsv', emit: embedding_manifest

    script:
    if (!batch_ids || !embedding_files) {
        error 'No embedding chunks provided; cannot merge embeddings'
    }
    // embedding_files are now staged in the work directory by Nextflow
    // Build manifest using staged file names
    def manifestLines = [batch_ids, embedding_files].transpose().collect { id, file ->
        "${id}\t${file}"
    }
    def manifestContent = manifestLines.join('\n')
    """
    cat <<'EOF' > embedding_chunks.tsv
${manifestContent}
EOF

    python3 ${baseDir}/bin/merge_embedding_chunks.py \
        --manifest embedding_chunks.tsv \
        --output node_embeddings.tsv
    """
}
