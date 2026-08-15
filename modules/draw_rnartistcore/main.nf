process DRAW_RNARTISTCORE {
    tag "${alignments.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0dddce1e576a7568d4b8f15af94dada49065a224ea16df56564820607d7bba22/data' :
        'community.wave.seqera.io/library/python_openjdk_rnartistcore:938f7ce45857eb48' }"

    input:
    path alignments

    output:
    path "plots_rnartistcore_*", emit: plots
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def plot_dir = "plots_rnartistcore_${alignments.baseName}"
    """
    mkdir -p ${plot_dir}
    plot_rnartistcore.py \\
        --alignments ${alignments} \\
        --outdir ${plot_dir} \\
        --highlight-colour '${params.plot_highlight_colour}' \\
        --max-pairs ${params.plot_max_pairs} \\
        --cpus ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnartistcore: 0.4.6
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """

    stub:
    def plot_dir = "plots_rnartistcore_${alignments.baseName}"
    """
    mkdir -p ${plot_dir}
    touch ${plot_dir}/.stub

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnartistcore: 0.4.6
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """
}
