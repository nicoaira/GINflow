process EMBED_RNA_GRAPHS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b4bdbc4205729d72e9cd7bba9ec992386928246f0fd0546f74d771adcd2b7b9/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/6585e5731a92a6725c71af35c8462cfd871fbd4e9798e8a8144eec43343062ee/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/ginfinity_pytorch-gpu_cuda-version:055fb9718f25c342' : 'community.wave.seqera.io/library/ginfinity:1.1.0--55be3b8c50d410c9') }"

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
