process BUILD_FAISS_INDEX {
    tag "faiss"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path windows, stageAs: 'windows/*'
    path manifests, stageAs: 'manifests/*'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'

    output:
    path "faiss", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    build_faiss.py \\
        --windows windows/*.windows.npz \\
        --manifests manifests/*.windows.manifest.json \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir faiss \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p faiss
    touch faiss/index.faiss
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    echo '{}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
