#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS {
    label params.use_gpu ? 'gpu' : null

    tag { "embeddings batch=${task.index} device=${params.use_gpu ? 'gpu' : 'cpu'}" }
    maxForks = 1
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ginflow-generate-embeddings:latest' :
        'nicoaira/ginflow-generate-embeddings:latest' }"

    input:
    val item   // either a batch_tsv path or a tuple [graphs_pt, metadata_tsv]

    output:
    path "embeddings_batch_${task.index}.tsv", emit: batch_embeddings

    script:
    // Determine device based on use_gpu parameter
    def DEVICE = params.use_gpu ? 'cuda' : 'cpu'

    def common_args = """ \\
      --id-column ${params.id_column} \\
      --output embeddings_batch_${task.index}.tsv \\
      --num-workers ${params.num_workers} \\
      --batch-size ${params.batch_size}
    """

    def cmd
    if (params.subgraphs) {
        // item is a tuple [graphs_pt, metadata_tsv]
        def graphs_pt_file    = item[0]
        def metadata_tsv_file = item[1]
        cmd = """
        echo "Using device: ${DEVICE}"
        
        ginfinity-embed \\
          --graph-pt ${graphs_pt_file} \\
          --meta-tsv ${metadata_tsv_file} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          --device ${DEVICE} \\
          ${common_args}
        """
    } else {
        // item is a batch_tsv file path
        def batch_tsv_file = item
        cmd = """
        echo "Using device: ${DEVICE}"
        
        ginfinity-embed \\
          --input ${batch_tsv_file} \\
          --structure-column-name ${params.structure_column_name} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          --device ${DEVICE} \\
          ${common_args}
        """
    }

    """
    echo "===== GENERATE_EMBEDDINGS: using device=${DEVICE} ====="
    ${cmd}
    """
}