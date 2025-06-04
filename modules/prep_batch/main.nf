#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PREP_BATCH {
    tag "prep_batch"

    conda "${moduleDir}/environment.yml"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/official-debian-bullseye-slim:latest' :
        'debian:bullseye-slim' }"
        
    input:
    val rows

    output:
    path "batch_${task.index}.tsv", emit: batch_file

    script:
    def header = rows[0].keySet().join("\t")
    def lines = rows.collect { r -> r.values().join("\t") }.join("\n")
    """
    cat > batch_${task.index}.tsv << 'EOF'
${header}
${lines}
EOF
    """
}