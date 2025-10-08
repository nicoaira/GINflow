#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ALIGN_CANDIDATES {
    tag "align_candidates"

    label 'high_cpu'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'alignment*'

    conda "${moduleDir}/environment.yml"

    input:
    path cluster_members
    path cluster_summaries
    path embeddings
    path meta

    output:
    path 'alignments.tsv', emit: alignments
    path 'alignment_stats.json', optional: true, emit: alignment_stats
    path 'alignment_dp.jsonl', optional: true, emit: alignment_dp
    path 'alignment_pairs.txt', optional: true, emit: alignment_text
    path 'alignment_plots', optional: true, emit: alignment_plots

    script:
    def argsSeq = params.sequence_column ? "--sequence-column ${params.sequence_column}" : ''
    def argsStruct = params.structure_column_name ? "--structure-column ${params.structure_column_name}" : ''
    def plotArgs = params.plot_scoring_matrices ? "--plot-scoring-matrices --plot-dir alignment_plots" : ''
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
        --output alignments.tsv \
        --stats-json alignment_stats.json \
        --dp-output alignment_dp.jsonl \
        --alignment-text alignment_pairs.txt \
        ${argsSeq} \
        ${argsStruct} \
        ${plotArgs}
    """
}
