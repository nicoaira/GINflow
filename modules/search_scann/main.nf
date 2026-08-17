process SEARCH_SCANN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data' :
        'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91' }"

    input:
    tuple val(meta), path(windows), path(manifest)
    path database

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def scann_lts = params.scann_leaves_to_search != null ? "--scann-leaves-to-search ${params.scann_leaves_to_search}" : ''
    """
    search_faiss.py \\
        --windows ${windows} \\
        --manifest ${manifest} \\
        --database ${database} \\
        --output ${prefix}.seeds.tsv \\
        --k ${params.seed_k} \\
        --min-similarity ${params.seed_min_similarity} \\
        --scann-reorder ${params.scann_reorder} \\
        ${scann_lts} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scann: \$(python3 -c "from importlib.metadata import version; print(version('scann'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scann: \$(python3 -c "from importlib.metadata import version; print(version('scann'))")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
