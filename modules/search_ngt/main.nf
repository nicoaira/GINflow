process SEARCH_NGT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_ngt:ff987f0bb9f59553' :
        'community.wave.seqera.io/library/python_numpy_ngt:9a0ca7a46e9c18b2' }"

    input:
    tuple val(meta), path(windows), path(manifest)
    path database

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_edges = params.ngt_max_no_of_edges != null ? "--ngt-max-no-of-edges ${params.ngt_max_no_of_edges}" : ''
    def search_epsilon = params.ngt_search_range_coefficient != null ? "--ngt-search-range-coefficient ${params.ngt_search_range_coefficient}" : ''
    def blob_epsilon = params.ngt_blob_search_range_coefficient != null ? "--ngt-blob-search-range-coefficient ${params.ngt_blob_search_range_coefficient}" : ''
    def search_radius = params.ngt_search_radius != null ? "--ngt-search-radius ${params.ngt_search_radius}" : ''
    def result_expansion = params.ngt_result_expansion != null ? "--ngt-result-expansion ${params.ngt_result_expansion}" : ''
    def edge_size = params.ngt_edge_size_for_search != null ? "--ngt-edge-size-for-search ${params.ngt_edge_size_for_search}" : ''
    """
    search_faiss.py \\
        --windows ${windows} \\
        --manifest ${manifest} \\
        --database ${database} \\
        --output ${prefix}.seeds.tsv \\
        --k ${params.seed_k} \\
        --min-similarity ${params.seed_min_similarity} \\
        --ngt-num-of-search-objects ${params.ngt_num_of_search_objects} \\
        --ngt-exploration-size ${params.ngt_exploration_size} \\
        --ngt-exact-result-expansion ${params.ngt_exact_result_expansion} \\
        --ngt-num-of-probes ${params.ngt_num_of_probes} \\
        ${max_edges} \\
        ${edge_size} \\
        ${search_epsilon} \\
        ${blob_epsilon} \\
        ${search_radius} \\
        ${result_expansion} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngt: \$(python3 -c "from importlib.metadata import version; print(version('ngt'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngt: \$(python3 -c "from importlib.metadata import version; print(version('ngt'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
