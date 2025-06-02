#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS {
    tag { "embeddings batch=${task.index}" }
    maxForks = 1

    input:
    val item   // either a batch_tsv path or a tuple [graphs_pt, metadata_tsv]

    output:
    path "embeddings_batch_${task.index}.tsv", emit: batch_embeddings

    script:
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
        # Detect if CUDA is available and use it
        DEVICE=\$(python3 -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')
        echo "Using detected device: \$DEVICE"
        
        ginfinity-embed \\
          --graph-pt ${graphs_pt_file} \\
          --meta-tsv ${metadata_tsv_file} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          --device \$DEVICE \\
          ${common_args}
        """
    } else {
        // item is a batch_tsv file path
        def batch_tsv_file = item
        cmd = """
        # Detect if CUDA is available and use it
        DEVICE=\$(python3 -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')
        echo "Using detected device: \$DEVICE"
        
        ginfinity-embed \\
          --input ${batch_tsv_file} \\
          --structure-column-name ${params.structure_column_name} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          --device \$DEVICE \\
          ${common_args}
        """
    }

    """
    echo "===== DEBUG: inside GENERATE_EMBEDDINGS ====="
    echo "NVIDIA_VISIBLE_DEVICES=\${NVIDIA_VISIBLE_DEVICES:-not set}"
    echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-not set}"

    # Check PyTorch GPU availability directly (nvidia-smi is not required)
    python3 - << 'PY'
import torch
print("==== DEBUG: torch.cuda.is_available() =>", torch.cuda.is_available())
if torch.cuda.is_available():
    print("==== Number of CUDA devices:", torch.cuda.device_count())
    print("==== CUDA Device Name:", torch.cuda.get_device_name(0))
PY

    ${cmd}
    """
}