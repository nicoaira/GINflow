#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_EMBEDDINGS_FROM_WINDOWS {
    tag { "embeddings ${metadata_tsv_file.baseName} device=${params.use_gpu ? 'gpu' : 'cpu'}" }

    label params.use_gpu ? 'gpu' : 'cpu'
    maxForks = 2
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity:latest' }"

    input:
    tuple path(graphs_pt_file), path(metadata_tsv_file)

    output:
    path "embeddings.tsv", emit: batch_embeddings

    script:
    def DEVICE = params.use_gpu ? 'cuda' : 'cpu'
    def OUTFILE = "embeddings.tsv"

    def common_args_lines = [
        "--id-column ${params.id_column}",
        "--output ${OUTFILE}",
        "--num-workers ${params.num_workers}",
        "--batch-size ${params.inference_batch_size}"
    ]
    if (params.ginfinity_model_path) {
        common_args_lines << "--model-path ${params.ginfinity_model_path}"
    }
    def common_args = common_args_lines.collect { "      ${it}" }.join(" \\\n")

    """
    echo "===== GENERATE_EMBEDDINGS_FROM_WINDOWS: using device=${DEVICE} ====="
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
        echo "===== GENERATE_EMBEDDINGS_FROM_WINDOWS: no windows after masking; emitted empty embeddings.tsv ====="
    fi
    """.stripIndent()
}
