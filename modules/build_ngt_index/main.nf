process BUILD_NGT_INDEX {
    tag "ngt"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_ngt:ff987f0bb9f59553' :
        'community.wave.seqera.io/library/python_numpy_ngt:9a0ca7a46e9c18b2' }"

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
        --backend ngt \\
        --ngt-index-type ${params.ngt_index} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngt: \$(python3 -c "from importlib.metadata import version; print(version('ngt'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p faiss/ngt
    touch faiss/ngt/.stub
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    echo '{}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngt: \$(python3 -c "from importlib.metadata import version; print(version('ngt'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
