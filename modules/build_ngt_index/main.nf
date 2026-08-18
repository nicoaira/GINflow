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
    def max_edges = params.ngt_max_no_of_edges != null ? "--ngt-max-no-of-edges ${params.ngt_max_no_of_edges}" : ''
    def search_epsilon = params.ngt_search_range_coefficient != null ? "--ngt-search-range-coefficient ${params.ngt_search_range_coefficient}" : ''
    def blob_epsilon = params.ngt_blob_search_range_coefficient != null ? "--ngt-blob-search-range-coefficient ${params.ngt_blob_search_range_coefficient}" : ''
    def search_radius = params.ngt_search_radius != null ? "--ngt-search-radius ${params.ngt_search_radius}" : ''
    def result_expansion = params.ngt_result_expansion != null ? "--ngt-result-expansion ${params.ngt_result_expansion}" : ''
    def qg_subvectors = params.ngt_qg_subvector_dimensions != null ? "--ngt-qg-subvector-dimensions ${params.ngt_qg_subvector_dimensions}" : ''
    def qbg_subvectors = params.ngt_qbg_subvectors != null ? "--ngt-qbg-subvectors ${params.ngt_qbg_subvectors}" : ''
    """
    build_faiss.py \\
        --windows windows/*.windows.npz \\
        --manifests manifests/*.windows.manifest.json \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir faiss \\
        --backend ngt \\
        --ngt-index-type ${params.ngt_index} \\
        --ngt-edge-size-for-creation ${params.ngt_edge_size_for_creation} \\
        --ngt-edge-size-for-search ${params.ngt_edge_size_for_search} \\
        --ngt-num-threads ${params.ngt_num_threads} \\
        --ngt-num-of-search-objects ${params.ngt_num_of_search_objects} \\
        --ngt-exploration-size ${params.ngt_exploration_size} \\
        --ngt-exact-result-expansion ${params.ngt_exact_result_expansion} \\
        --ngt-num-of-probes ${params.ngt_num_of_probes} \\
        --ngt-qbg-cluster-data-type ${params.ngt_qbg_cluster_data_type} \\
        ${max_edges} \\
        ${search_epsilon} \\
        ${blob_epsilon} \\
        ${search_radius} \\
        ${result_expansion} \\
        ${qg_subvectors} \\
        ${qbg_subvectors} \\
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
