process ESTIMATE_EVD {
    tag "evd"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e74cfb57864deb17f17a966196bf9cde064bb56aa8df1fb6fe2d906008ba4148/data' :
        'community.wave.seqera.io/library/python_ginfinity-sw:b934c840797fac86' }"

    input:
    path database
    path parameters

    output:
    path "evd.json",     emit: evd
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    estimate_evd.py \\
        --database ${database} \\
        --parameters ${parameters} \\
        --output evd.json \\
        --samples ${params.evd_samples} \\
        --max-length ${params.evd_max_length} \\
        --max-cells ${params.align_max_cells} \\
        --max-alignments ${params.align_max_alignments} \\
        --min-score ${params.align_min_score} \\
        --min-match-count ${params.align_min_match_count} \\
        --seed ${params.evd_seed} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo '{"lambda": 0.1, "K": 0.1, "database_residues": 1, "n_samples": 0}' > evd.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
