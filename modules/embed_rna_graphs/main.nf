process EMBED_RNA_GRAPHS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5889c2a8c201d0ceb8323fe5e2409eb86bfd97eb5ab96f066058a3c01521524/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52aef85397c34ad4b14f3068d1881c5b09fd9df24ae5697900bb598e2dc892d0/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/ginfinity_pytorch-gpu_cuda-version:d36bdbe90847a4ec' : 'community.wave.seqera.io/library/ginfinity:1.0.1--ec3ad584ab4c66db') }"

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
