process WRITE_REPORT {
    tag "report"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path alignments
    path alignment_text
    path evd
    path clusters
    path seeds
    path plots_rnartist, stageAs: 'plots_rnartist'
    path plots_r4rna,    stageAs: 'plots_r4rna'
    path plots_sw,       stageAs: 'plots_sw'

    output:
    path "report.html", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def rn_arg = plots_rnartist instanceof Collection && plots_rnartist.isEmpty() ? '' : "--plots-rnartist ${plots_rnartist}"
    def r4_arg = plots_r4rna instanceof Collection && plots_r4rna.isEmpty() ? '' : "--plots-r4rna ${plots_r4rna}"
    def sw_arg = plots_sw instanceof Collection && plots_sw.isEmpty() ? '' : "--plots-sw ${plots_sw}"
    """
    write_report.py \\
        --alignments ${alignments} \\
        --alignment-text ${alignment_text} \\
        --evd ${evd} \\
        --clusters ${clusters} \\
        --seeds ${seeds} \\
        ${rn_arg} \\
        ${r4_arg} \\
        ${sw_arg} \\
        --highlight-colour '${params.plot_highlight_colour}' \\
        --logo ${projectDir}/docs/images/ginfinity.png \\
        --output report.html \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """

    stub:
    """
    echo '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/><title>ginflow search report</title></head><body><h1>Search report</h1></body></html>' > report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """
}
