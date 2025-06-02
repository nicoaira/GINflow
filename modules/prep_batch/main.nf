#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PREP_BATCH {
    tag "prep_batch"
    
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