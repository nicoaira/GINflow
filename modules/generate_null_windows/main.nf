#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_NULL_WINDOWS {
    tag { "generate_null_windows_batch_${ids.size()}" }

    label 'high_cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity@sha256:a2a43687b5583700d6fb0376de3c24f3dcc81ab10ea86427d51877aef8368ed6' }"

    // initialize BLAS thread limits before running the script
    beforeScript 'export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1'

    input:
    tuple val(ids), path(meta_tsv)

    output:
    tuple val(ids), path("windows_graphs.pt"), path("windows_metadata.tsv"), emit: window_files_null

    script:
    """
    ginfinity-generate-windows \\
      --input ${meta_tsv} \\
      --output-dir . \\
      --id-column ${params.id_column} \\
      --structure-column-name ${params.structure_column_name} \\
      --L ${params.L} \\
      ${params.keep_paired_neighbors ? '--keep-paired-neighbors' : ''} \\
      --mask-threshold ${params.mask_threshold} \\
      --num-workers ${task.cpus}
    """
}
