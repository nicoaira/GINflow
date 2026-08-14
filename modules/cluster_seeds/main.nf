process CLUSTER_SEEDS {
    tag "clusters"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path seeds

    output:
    path "clusters.tsv",        emit: clusters
    path "cluster_members.tsv", emit: members
    path "cluster_stats.json",  emit: stats
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cluster_seeds.py \\
        --seeds ${seeds} \\
        --output-clusters clusters.tsv \\
        --output-members cluster_members.tsv \\
        --cluster-span ${params.cluster_span} \\
        --min-seeds ${params.cluster_min_seeds} \\
        --diagonal-tolerance ${params.cluster_diagonal_tolerance} \\
        --max-diagonal-span ${params.cluster_max_diagonal_span} \\
        --max-seed-rank ${params.cluster_max_seed_rank} \\
        --stats-json cluster_stats.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo -e "cluster_id\\tquery_id\\ttarget_id\\tseed_count\\tquery_start\\tquery_end\\ttarget_start\\ttarget_end\\tmax_score\\tdiagonal\\tdiagonal_min\\tdiagonal_max\\tdiagonal_span" > clusters.tsv
    echo -e "cluster_id\\tquery_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank\\tdiagonal" > cluster_members.tsv
    echo '{}' > cluster_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
