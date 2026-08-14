process BUILD_RNA_GRAPHS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52aef85397c34ad4b14f3068d1881c5b09fd9df24ae5697900bb598e2dc892d0/data' :
        'community.wave.seqera.io/library/ginfinity:1.0.1--ec3ad584ab4c66db' }"

    input:
    tuple val(meta), path(structures)

    output:
    tuple val(meta), path("*.safetensors"), path("*.json"), emit: graphs
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ginfinity \\
        build-graphs \\
        --input ${structures} \\
        --output ${prefix}.safetensors \\
        --metadata ${prefix}.json \\
        --id-column ${params.id_column} \\
        --sequence-column ${params.sequence_column} \\
        --structure-column ${params.structure_column} \\
        --delimiter '${params.delimiter}' \\
        --checksum \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity: \$(ginfinity --version 2>&1 | sed 's/^ginfinity //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.safetensors
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity: \$(ginfinity --version 2>&1 | sed 's/^ginfinity //')
    END_VERSIONS
    """
}
