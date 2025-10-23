#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ALIGN_CANDIDATES {
    tag { "query_" + cluster_members.baseName.replaceAll(/^cluster_members_/, '') }

    label 'high_cpu'

    publishDir "${params.outdir}/query_results", mode: 'copy', pattern: 'alignment*',
        saveAs: { filename ->
            def qid = filename.replaceAll(/^alignment(s|_stats|_dp|_pairs|_plots)_/, '').replaceAll(/\.(tsv|json|jsonl|txt)$/, '')
            "${qid}/${filename}"
        }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-align-candidates:latest' :
        'docker.io/nicoaira/ginflow-align-candidates:latest' }"

    input:
    path cluster_members
    path cluster_summaries
    path embeddings
    path meta
    path evd_params, stageAs: 'evd_params.json'

    output:
    path 'alignments_*.tsv', emit: alignments
    path 'alignment_stats_*.json', optional: true, emit: alignment_stats
    path 'alignment_dp_*.jsonl', optional: true, emit: alignment_dp
    path 'alignment_pairs_*.txt', optional: true, emit: alignment_text
    path 'alignment_plots_*', optional: true, emit: alignment_plots, type: 'dir'

    script:
    def query_id = cluster_members.baseName.replaceAll(/^cluster_members_/, '')
    def argsSeq = params.sequence_column ? "--sequence-column ${params.sequence_column}" : ''
    def argsStruct = params.structure_column_name ? "--structure-column ${params.structure_column_name}" : ''
    def plotArgs = params.plot_scoring_matrices ? "--plot-scoring-matrices --plot-dir alignment_plots_${query_id}" : ''
    def evalueArgs = (params.calculate_evalue == false) ? "--no-calculate-evalue" : ""
    def evdParamsArg = evd_params.name != 'NO_FILE' ? "--evd-params-json evd_params.json" : ""
    """
    python3 ${baseDir}/bin/align_candidates.py \
        --cluster-members ${cluster_members} \
        --cluster-summaries ${cluster_summaries} \
        --embeddings ${embeddings} \
        --meta ${meta} \
        --id-column ${params.id_column} \
        --node-index-column node_index \
        --embedding-column embedding_vector \
        --background-samples ${params.background_samples ?: 10000} \
        --random-seed ${params.random_seed ?: 42} \
        --gamma ${params.alignment_gamma ?: 1.5} \
        --band-width ${params.band_width ?: 96} \
        --band-buffer ${params.band_buffer ?: 32} \
        --band-max-width ${params.band_max_width ?: 0} \
        --xdrop ${params.xdrop ?: 50} \
        --gap-open ${params.gap_open ?: 12} \
        --gap-extend ${params.gap_extend ?: 2} \
        --padding ${params.alignment_padding ?: 32} \
        --score-min ${params.score_min ?: -4} \
        --score-max ${params.score_max ?: 8} \
        --top-n ${params.top_n ?: 50} \
        --evd-samples ${params.evd_samples ?: 1000} \
        --workers ${task.cpus} \
        --output alignments_${query_id}.tsv \
        --stats-json alignment_stats_${query_id}.json \
        --dp-output alignment_dp_${query_id}.jsonl \
        --alignment-text alignment_pairs_${query_id}.txt \
        ${argsSeq} \
        ${argsStruct} \
        ${evalueArgs} \
        ${evdParamsArg} \
        ${plotArgs}
    """
}
