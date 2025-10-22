process CALCULATE_STATISTICS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas conda-forge::scipy" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-92153046697e926078267aa93de333f67407a99c:2737c20743a118374b0f352744b373b4e3103e59-0' :
        'nicoaira/ginflow-calculate-statistics:latest' }"

    input:
    tuple val(meta), path(alignments)
    val db_size
    val loc
    val scale

    output:
    tuple val(meta), path("*.tsv"), emit: alignments_with_evalue

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    calculate_statistics.py \\
        --alignments_tsv $alignments \\
        --output_tsv ${prefix}.tsv \\
        --db_size $db_size \\
        --loc $loc \\
        --scale $scale \\
        $args
    """
}
