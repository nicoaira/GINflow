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
    """
    search_faiss.py \\
        --windows ${windows} \\
        --manifest ${manifest} \\
        --database ${database} \\
        --output ${prefix}.seeds.tsv \\
        --k ${params.seed_k} \\
        --min-similarity ${params.seed_min_similarity} \\
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
