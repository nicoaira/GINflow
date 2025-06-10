#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS {
    tag { item instanceof List ? "embeddings ${item[1].baseName} device=${params.use_gpu ? 'gpu' : 'cpu'}" : "embeddings ${item.baseName} device=${params.use_gpu ? 'gpu' : 'cpu'}" }
    
    label params.use_gpu ? 'gpu' : 'cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.ion/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity@sha256:a2a43687b5583700d6fb0376de3c24f3dcc81ab10ea86427d51877aef8368ed6' }"

    input:
    val item   // either a batch_tsv path or a tuple [graphs_pt, metadata_tsv]

    output:
    path "embeddings.tsv", emit: batch_embeddings

    script:
    def DEVICE = params.use_gpu ? 'cuda' : 'cpu'
    def OUTFILE = "embeddings.tsv"

    def common_args = """ \\
      --id-column ${params.id_column} \\
      --output ${OUTFILE} \\
      --num-workers ${params.num_workers} \\
      --batch-size ${params.inference_batch_size}
    """

    def cmd
    if (item instanceof List) {
        // subgraphs mode: item is a tuple [graphs_pt, metadata_tsv]
        def graphs_pt_file    = item[0]
        def metadata_tsv_file = item[1]
        cmd = """
        ginfinity-embed \\
          --graph-pt ${graphs_pt_file} \\
          --meta-tsv ${metadata_tsv_file} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          --device ${DEVICE} \\
          ${common_args}
        """
    } else {
        // direct mode: item is a batch_tsv file path
        def batch_tsv_file = item
        cmd = """
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
    python - << 'PY'
import torch
print('TORCH_CUDA_IS_AVAILABLE=' + str(torch.cuda.is_available()))
PY
    ${cmd}
    """
}