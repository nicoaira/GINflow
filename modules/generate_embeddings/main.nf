#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS {
    tag { item instanceof List ? "embeddings ${item[1].baseName} device=${params.use_gpu ? 'gpu' : 'cpu'}" : "embeddings ${item.baseName} device=${params.use_gpu ? 'gpu' : 'cpu'}" }
    
    label params.use_gpu ? 'gpu' : 'cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity:latest' }"

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
        graphs_pt_file    = item[0]
        metadata_tsv_file = item[1]
        cmd = """
        STATUS=\$(python - <<'PY'
from pathlib import Path

meta_path = Path('${metadata_tsv_file}')
out_path = Path('embeddings.tsv')

header_line = ''
has_windows = False
if meta_path.exists():
    with meta_path.open('r') as handle:
        header_line = handle.readline()
        for line in handle:
            if line.strip():
                has_windows = True
                break
else:
    print('RUN_EMBED')
    raise SystemExit(0)

if has_windows:
    print('RUN_EMBED')
    raise SystemExit(0)

header = [item for item in header_line.strip().split('\\t') if item]
if not header:
    header = ['window_id', '${params.id_column}', 'window_start', 'window_end']
if 'embedding_vector' not in header:
    try:
        idx = header.index('window_end') + 1
    except ValueError:
        idx = len(header)
    header = header[:idx] + ['embedding_vector'] + header[idx:]

with out_path.open('w', newline='') as out_handle:
    out_handle.write('\\t'.join(header) + '\\n')

print('SKIP_EMBED')
PY
        )
        if [[ "\$STATUS" == "RUN_EMBED" ]]; then
            python - << 'PY'
import torch
print('TORCH_CUDA_IS_AVAILABLE=' + str(torch.cuda.is_available()))
PY
            ginfinity-embed \\
              --graph-pt ${graphs_pt_file} \\
              --meta-tsv ${metadata_tsv_file} \\
              --keep-cols ${params.id_column},window_start,window_end \\
              --device ${DEVICE} \\
              ${common_args}
        else
            echo "===== GENERATE_EMBEDDINGS: no windows after masking; emitted empty embeddings.tsv ====="
        fi
        """.stripIndent()
    } else {
        // direct mode: item is a batch_tsv file path
        batch_tsv_file = item
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
    """.stripIndent()
}
