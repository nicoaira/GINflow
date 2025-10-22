process FIT_GUMBEL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas conda-forge::scipy conda-forge::numpy" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-15846967b28854b833d28a818a3dd6f5d5832276:a636322415420b977f56964b1593f3333c9935d3-0' :
        'nicoaira/ginflow-fit-gumbel:latest' }"

    input:
    tuple val(meta), path(alignment_stats)

    output:
    tuple val(meta), path("*.json"), emit: gumbel_params

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_gumbel_params"
    """
    fit_gumbel.py \\
        --alignment_stats $alignment_stats \\
        --output_json ${prefix}.json \\
        $args
    """
}
