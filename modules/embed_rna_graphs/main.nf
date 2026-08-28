process EMBED_RNA_GRAPHS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/be8e6a92287e437f6ad12a1ec6c00ee175dad4d3aa20d8aff902c740e3ae1ed4/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94aa66bb7d6f358fd34bdcc02f0cd69a44ca39a2bcee482e047a7ece445174a/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/ginfinity_pytorch-gpu_cuda-version:9442b2ea6d70a13b' : 'community.wave.seqera.io/library/ginfinity:1.2.2--4f957fe3f833bed2') }"

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
