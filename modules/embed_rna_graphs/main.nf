process EMBED_RNA_GRAPHS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30bd42d8856c1e86b0d195aa80b1100979719942fc21b56759e3b7735e28b6dd/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9af1d5c2bf9cb29fac6bef4b49ac1c9d3bf0c0885dc4f8f3f44109e8832737bb/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/ginfinity_pytorch-gpu_cuda-version:2bfaf62f521c8200' : 'community.wave.seqera.io/library/ginfinity:1.2.1--18305a26eeaa0fcb') }"

    input:
    tuple val(meta), path(graphs), path(metadata)

    output:
    tuple val(meta), path("*.npz"), path("*.manifest.json"), emit: embeddings
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def use_gpu   = task.accelerator
    def device    = use_gpu ? 'cuda' : 'cpu'
    def cuda_flag = use_gpu ? '--allow-nondeterministic-cuda' : ''
    def precision_flag = params.ginfinity_full_precision ? '--full-precision' : ''
    """
    ginfinity \\
        embed-graphs \\
        --input ${graphs} \\
        --metadata ${metadata} \\
        --output ${prefix}.npz \\
        --manifest ${prefix}.manifest.json \\
        --device ${device} \\
        --max-batch-nodes ${params.max_batch_nodes} \\
        --max-batch-edges ${params.max_batch_edges} \\
        --verify-checksum \\
        --checksum \\
        ${cuda_flag} \\
        ${precision_flag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity: \$(ginfinity --version 2>&1 | sed 's/^ginfinity //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.npz
    touch ${prefix}.manifest.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity: \$(ginfinity --version 2>&1 | sed 's/^ginfinity //')
    END_VERSIONS
    """
}
