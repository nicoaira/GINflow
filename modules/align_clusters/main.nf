process ALIGN_CLUSTERS {
    tag "${clusters.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e74cfb57864deb17f17a966196bf9cde064bb56aa8df1fb6fe2d906008ba4148/data' :
        'community.wave.seqera.io/library/python_ginfinity-sw:b934c840797fac86' }"

    input:
    path clusters
    path query_embeddings, stageAs: 'query_emb/*'
    path query_metadata, stageAs: 'query_meta/*'
    path database
    path parameters
    path evd

    output:
    path "*.alignments.tsv",          emit: alignments
    path "*.alignments.txt",          emit: text
    path "*.alignment_stats.json",    emit: stats
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: clusters.baseName
    """
    align_clusters.py \\
        --clusters ${clusters} \\
        --parameters ${parameters} \\
        --query-embeddings query_emb/*.npz \\
        --query-metadata query_meta/*.json \\
        --database ${database} \\
        --evd ${evd} \\
        --pad ${params.align_pad} \\
        --max-cells ${params.align_max_cells} \\
        --output ${prefix}.alignments.tsv \\
        --alignment-text ${prefix}.alignments.txt \\
        --stats-json ${prefix}.alignment_stats.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: clusters.baseName
    """
    echo -e "cluster_id\\tquery_id\\ttarget_id\\tscore\\ttotal_score\\tmax_score\\tbit_score\\tevalue\\tevalue_pair\\tlog_evalue\\talignment_count\\tquery_start\\tquery_end\\ttarget_start\\ttarget_end\\tquery_length\\ttarget_length\\tmatch_count\\taligned_columns\\tungapped_columns\\tbase_matches\\tstructure_matches\\tseed_count\\tmax_seed_score\\tquery_sequence\\tquery_structure\\ttarget_sequence\\ttarget_structure\\tquery_aligned\\ttarget_aligned" > ${prefix}.alignments.tsv
    touch ${prefix}.alignments.txt
    echo '{}' > ${prefix}.alignment_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
