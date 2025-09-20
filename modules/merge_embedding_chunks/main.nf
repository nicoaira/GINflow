#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_EMBEDDING_CHUNKS {
    tag 'merge_embedding_chunks'

    label 'high_cpu'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*'

    conda "${moduleDir}/environment.yml"

    input:
    val chunk_records

    output:
    path 'node_embeddings.tsv', emit: node_embeddings
    path 'embedding_chunks.tsv', emit: embedding_manifest

    script:
    if (!chunk_records) {
        error 'No embedding chunks provided; cannot merge embeddings'
    }
    def manifestLines = chunk_records.collect { item -> "${item[0]}\t${item[1]}" }
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
