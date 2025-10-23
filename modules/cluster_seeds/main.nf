#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CLUSTER_SEEDS {
    tag { "query_" + seeds.baseName.replaceAll(/^seeds_/, '') }

    label 'lightweight'

    publishDir "${params.outdir}/query_results", mode: 'copy', pattern: 'cluster*.{tsv,json}',
        saveAs: { filename ->
            def qid = filename.replaceAll(/^cluster(s|_members|_stats)_/, '').replaceAll(/\.(tsv|json)$/, '')
            "${qid}/${filename}"
        }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-cluster-seeds:latest' :
        'docker.io/nicoaira/ginflow-cluster-seeds:latest' }"

    input:
    path seeds

    output:
    path 'clusters_*.tsv', emit: clusters
    path 'cluster_members_*.tsv', emit: cluster_members
    path 'cluster_stats_*.json', optional: true, emit: cluster_stats

    script:
    def query_id = seeds.baseName.replaceAll(/^seeds_/, '')
    def span = params.cluster_span ?: 80
    def minSeeds = params.cluster_min_seeds ?: 2
    def diagTol = params.cluster_diagonal_tolerance ?: 12
    def maxDiagSpan = params.cluster_max_diagonal_span ?: 96
    """
    python3 ${baseDir}/bin/cluster_seeds.py \
        --seeds ${seeds} \
        --output-clusters clusters_${query_id}.tsv \
        --output-members cluster_members_${query_id}.tsv \
        --cluster-span ${span} \
        --min-seeds ${minSeeds} \
        --diagonal-tolerance ${diagTol} \
        --max-diagonal-span ${maxDiagSpan} \
        --stats-json cluster_stats_${query_id}.json
    """
}
