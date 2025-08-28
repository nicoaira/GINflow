#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS_NULL {
    tag { "embeddings_null ${new_id} device=${params.use_gpu ? 'gpu' : 'cpu'}" }

    label params.use_gpu ? 'gpu' : 'cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity:latest' }"

    input:
    tuple val(new_id), path(meta_tsv)

    output:
    tuple val(new_id), path("embeddings.tsv"), emit: embeddings_for_query

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
    echo "===== GENERATE_EMBEDDINGS_NULL: ${new_id} using device=${DEVICE} ====="
    python - << 'PY'
import torch
print('TORCH_CUDA_IS_AVAILABLE=' + str(torch.cuda.is_available()))
PY
    if ${params.subgraphs} ; then
      echo "[NULL] Generating windows for ${new_id}..."
      ginfinity-generate-windows \
        --input ${meta_tsv} \
        --output-dir . \
        --id-column ${params.id_column} \
        --structure-column-name ${params.structure_column_name} \
        --L ${params.L} \
        ${params.keep_paired_neighbors ? '--keep-paired-neighbors' : ''} \
        --mask-threshold ${params.mask_threshold} \
        --num-workers ${task.cpus}

      echo "[NULL] Embedding windows for ${new_id}..."
      ginfinity-embed \
        --graph-pt windows_graphs.pt \
        --meta-tsv windows_metadata.tsv \
        --keep-cols ${params.id_column},window_start,window_end \
        --device ${DEVICE} \
        ${common_args}
    else
      echo "[NULL] Embedding full sequence for ${new_id}..."
      ginfinity-embed \
        --input ${meta_tsv} \
        --structure-column-name ${params.structure_column_name} \
        --keep-cols ${params.id_column},window_start,window_end \
        --device ${DEVICE} \
        ${common_args}
    fi
    """
}
