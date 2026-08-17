process BUILD_SCANN_INDEX {
    tag "scann"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data' :
        'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91' }"

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
    def args         = task.ext.args ?: ''
    def scann_soar   = params.scann_soar ? '--scann-soar' : ''
    def scann_leaves = params.scann_leaves != null ? "--scann-leaves ${params.scann_leaves}" : ''
    def scann_lts    = params.scann_leaves_to_search != null ? "--scann-leaves-to-search ${params.scann_leaves_to_search}" : ''
    """
    build_faiss.py \\
        --windows windows/*.windows.npz \\
        --manifests manifests/*.windows.manifest.json \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir faiss \\
        --backend scann \\
        --index-type scann \\
        --scann-reorder ${params.scann_reorder} \\
        --scann-ah-dim ${params.scann_ah_dim} \\
        --scann-anisotropic ${params.scann_anisotropic} \\
        --num-neighbors ${params.seed_k} \\
        ${scann_soar} \\
        ${scann_leaves} \\
        ${scann_lts} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scann: \$(python3 -c "from importlib.metadata import version; print(version('scann'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p faiss/scann
    touch faiss/scann/.stub
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    echo '{}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scann: \$(python3 -c "from importlib.metadata import version; print(version('scann'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
