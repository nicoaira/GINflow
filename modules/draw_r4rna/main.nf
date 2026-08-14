process DRAW_R4RNA {
    tag "r4rna"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b78862b3f32b7102418639ef10cc42a0f24088e9e3f7a42e5a7a78ca653e58c/data' :
        'community.wave.seqera.io/library/python_r-r4rna:bee0f0a84f183a3c' }"

    input:
    path alignments

    output:
    path "plots_r4rna",  emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p plots_r4rna
    plot_r4rna.py \\
        --alignments ${alignments} \\
        --outdir plots_r4rna \\
        --highlight-colour '${params.plot_highlight_colour}' \\
        --max-plots ${params.plot_max} \\
        --cpus ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-r4rna: 2.0.9
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p plots_r4rna
    touch plots_r4rna/.stub

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-r4rna: 2.0.9
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """
}
