process DRAW_SW {
    tag "${alignments.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7db9e2e04ab9a21b4e3bf856e5c3a9e3eb1909cf3e251b366dca6e10ae8684a/data' :
        'community.wave.seqera.io/library/python_ginfinity-sw:151e772020737622' }"

    input:
    path alignments
    path clusters
    path query_embeddings, stageAs: 'query_emb/*'
    path database
    path parameters

    output:
    path "plots_sw_*",    emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def plot_dir = "plots_sw_${alignments.baseName}"
    """
    mkdir -p ${plot_dir}
    plot_sw.py \\
        --alignments ${alignments} \\
        --clusters ${clusters} \\
        --parameters ${parameters} \\
        --query-embeddings query_emb/*.npz \\
        --database ${database} \\
        --outdir ${plot_dir} \\
        --pad ${params.align_pad} \\
        --max-cells ${params.align_max_cells} \\
        --max-pairs ${params.plot_max_pairs} \\
        --highlight-colour '${params.plot_highlight_colour}' \\
        --cpus ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def plot_dir = "plots_sw_${alignments.baseName}"
    """
    mkdir -p ${plot_dir}
    touch ${plot_dir}/.stub

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
