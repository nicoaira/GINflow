process ESTIMATE_EVD {
    tag "evd"
    cpus { Math.max(1, params.align_cpus.toString().toInteger()) }
    time { 8.h * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7db9e2e04ab9a21b4e3bf856e5c3a9e3eb1909cf3e251b366dca6e10ae8684a/data' :
        'community.wave.seqera.io/library/python_ginfinity-sw:151e772020737622' }"

    input:
    path database
    val align_mu
    val align_sigma
    val align_gamma
    val align_score_min
    val align_score_max
    val align_gap_open
    val align_gap_extend
    val align_score_offset

    output:
    path "evd.json",     emit: evd
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export NUMBA_NUM_THREADS=${task.cpus}
    estimate_evd.py \\
        --database ${database} \\
        --output evd.json \\
        --mu ${align_mu} \\
        --sigma ${align_sigma} \\
        --gamma ${align_gamma} \\
        --score-min ${align_score_min} \\
        --score-max ${align_score_max} \\
        --gap-open ${align_gap_open} \\
        --gap-extend ${align_gap_extend} \\
        --score-offset ${align_score_offset} \\
        --samples ${params.evd_samples} \\
        --max-length ${params.evd_max_length} \\
        --max-cells ${params.align_max_cells} \\
        --max-alignments ${params.align_max_alignments} \\
        --min-score ${params.align_min_score} \\
        --min-match-count ${params.align_min_match_count} \\
        --seed ${params.evd_seed} \\
        --workers ${task.cpus} \\
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
