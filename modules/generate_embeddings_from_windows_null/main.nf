#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS_FROM_WINDOWS_NULL {
    tag { "embeddings_null_windows_batch_${ids.size()} device=${params.use_gpu ? 'gpu' : 'cpu'}" }

    label params.use_gpu ? 'gpu' : 'cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity:latest' }"

    input:
    tuple val(ids), path(graphs_pt), path(metadata_tsv)

    output:
    tuple val(ids), path("embeddings.tsv"), emit: embeddings_for_query

    script:
    def DEVICE = params.use_gpu ? 'cuda' : 'cpu'
    def OUTFILE = "embeddings.tsv"

    def common_args = """ \\
      --id-column ${params.id_column} \\
      --output ${OUTFILE} \\
      --num-workers ${params.num_workers} \\
      --batch-size ${params.inference_batch_size}
    """

    """
    echo "===== GENERATE_EMBEDDINGS_FROM_WINDOWS_NULL: batch size=${ids.size()} device=${DEVICE} ====="
    python - << 'PY'
import torch
print('TORCH_CUDA_IS_AVAILABLE=' + str(torch.cuda.is_available()))
PY
    ginfinity-embed \
      --graph-pt ${graphs_pt} \
      --meta-tsv ${metadata_tsv} \
      --keep-cols ${params.id_column},window_start,window_end \
      --device ${DEVICE} \
      ${common_args}
    """
}
