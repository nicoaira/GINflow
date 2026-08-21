process FIT_NODE_QUANTIZER {
    tag "node-quantization"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path embeddings, stageAs: 'embeddings/*'
    path manifests, stageAs: 'manifests/*'

    output:
    path "quantization", emit: quantization
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    quantize_nodes.py \\
        --embeddings embeddings/*.npz \\
        --manifests manifests/*.manifest.json \\
        --outdir quantization \\
        --mode ${params.quantize} \\
        --pq-m ${params.pq_m} \\
        --pq-nbits ${params.pq_nbits} \\
        --sample-size ${params.pq_sample_size} \\
        --niter ${params.pq_niter} \\
        --opq-iters ${params.opq_iters} \\
        --seed ${params.pq_seed} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p quantization/nodes
    echo '{"mode":"${params.quantize}"}' > quantization/quantizer.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: stub
    END_VERSIONS
    """
}

process APPLY_NODE_QUANTIZER {
    tag "node-quantization-apply"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path embeddings, stageAs: 'embeddings/*'
    path manifests, stageAs: 'manifests/*'
    path quantizer, stageAs: 'quantizer'

    output:
    path "quantization", emit: quantization
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    quantize_nodes.py \\
        --embeddings embeddings/*.npz \\
        --manifests manifests/*.manifest.json \\
        --outdir quantization \\
        --mode ${params.quantize} \\
        --quantizer-dir quantizer \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p quantization/nodes
    echo '{"mode":"${params.quantize}","applied":true}' > quantization/quantizer.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: stub
    END_VERSIONS
    """
}
