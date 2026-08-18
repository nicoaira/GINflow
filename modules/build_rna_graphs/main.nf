process BUILD_RNA_GRAPHS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9af1d5c2bf9cb29fac6bef4b49ac1c9d3bf0c0885dc4f8f3f44109e8832737bb/data' :
        'community.wave.seqera.io/library/ginfinity:1.2.1--18305a26eeaa0fcb' }"

    input:
    tuple val(meta), path(structures)

    output:
    tuple val(meta), path("*.safetensors"), path("*.json"), emit: graphs
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args   ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def slice_args     = params.no_slices ? '--no-slices' : "--start-column ${params.start_column} --end-column ${params.end_column}"
    def neighbour_args = (!params.no_slices && params.keep_paired_neighbours) ? "--keep-paired-neighbours --context-hops ${params.context_hops}" : ''
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
        ${slice_args} \\
        ${neighbour_args} \\
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
