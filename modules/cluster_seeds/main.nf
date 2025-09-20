#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CLUSTER_SEEDS {
    tag "cluster_seeds"

    label 'lightweight'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'seeds*'

    conda "${moduleDir}/environment.yml"

    input:
    path seeds

    output:
    path 'clusters.tsv', emit: clusters
    path 'cluster_members.tsv', emit: cluster_members
    path 'cluster_stats.json', optional: true, emit: cluster_stats

    script:
    def span = params.cluster_span ?: 80
    def minSeeds = params.cluster_min_seeds ?: 2
    """
    python3 ${baseDir}/bin/cluster_seeds.py \
        --seeds ${seeds} \
        --output-clusters clusters.tsv \
        --output-members cluster_members.tsv \
        --cluster-span ${span} \
        --min-seeds ${minSeeds} \
        --stats-json cluster_stats.json
    """
}
