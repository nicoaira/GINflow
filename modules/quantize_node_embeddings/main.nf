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
    quantize_node_embeddings.py \\
        --embeddings embeddings/*.npz \\
        --manifests manifests/*.manifest.json \\
        --outdir quantization \\
        --k ${params.node_quantization_k} \\
        --sample-size ${params.node_quantization_sample_size} \\
        --niter ${params.node_quantization_niter} \\
        --seed ${params.node_quantization_seed} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)" 2>/dev/null || echo unavailable)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p quantization/nodes
    touch quantization/centroids.npy
    touch quantization/similarity.npy
    echo '{}' > quantization/quantization.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
