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
    tuple val(chunk_id), path(embedding_file) from chunk_records.collect()

    output:
    path 'node_embeddings.tsv', emit: node_embeddings
    path 'embedding_chunks.tsv', emit: embedding_manifest

    script:
    def manifestContent = chunk_id.zip(embedding_file).collect{ it.join('\t') }.join('\n')
    """
    cat <<'EOF' > embedding_chunks.tsv
${manifestContent}
EOF

    python3 ${baseDir}/bin/merge_embedding_chunks.py \
        --manifest embedding_chunks.tsv \
        --output node_embeddings.tsv
    """
}
