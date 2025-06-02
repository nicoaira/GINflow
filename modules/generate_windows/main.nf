#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_WINDOWS {
    tag "generate_windows"
    cpus 14

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ginflow-generate-windows:latest' :
        'nicoaira/ginflow-generate-windows:latest' }"

    // initialize BLAS thread limits before running the script
    beforeScript 'export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1'

    input:
    path orig_tsv

    output:
    tuple path("windows_graphs.pt"), path("windows_metadata.tsv"), emit: window_files

    script:
    """
    ginfinity-generate-windows \\
      --input ${orig_tsv} \\
      --output-dir . \\
      --id-column ${params.id_column} \\
      --structure-column-name ${params.structure_column_name} \\
      --L ${params.L} \\
      ${params.keep_paired_neighbors ? '--keep-paired-neighbors' : ''} \\
      --mask-threshold ${params.mask_threshold} \\
      --num-workers ${task.cpus}
    """
}