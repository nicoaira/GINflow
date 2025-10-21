#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_NODE_EMBEDDINGS {
    tag { "node_emb_${batch_id}" }

    label params.use_gpu ? 'gpu' : 'high_cpu'

    maxForks = 5

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'docker.io/nicoaira/ginfinity:latest' }"

    input:
    tuple val(batch_id), path(input_tsv)

    output:
    tuple val(batch_id), path('node_embeddings.tsv'), emit: node_embeddings

    script:
    def args_model = params.ginfinity_model_path ? "--model-path ${params.ginfinity_model_path}" : null
    def args_keep  = params.keep_cols ? "--keep-cols ${params.keep_cols}" : null
    def args_struct = params.structure_column_name ? "--structure-column-name ${params.structure_column_name}" : null
    def device = params.use_gpu ? 'cuda' : 'cpu'
    def workers = params.num_workers ?: 4
    def batch = params.split_size ?: (params.inference_batch_size ?: 1024)
    def optionalArgs = [args_struct, args_model, args_keep].findAll { it }
    def cliLines = [
        "--input ${input_tsv}",
        "--output raw_node_embeddings.tsv",
        "--id-column ${params.id_column}",
        "--device ${device}",
        "--num-workers ${workers}"
    ] + optionalArgs + ["--batch-size ${batch}"]
    def cliBlock = cliLines.collect { "        ${it}" }.join(' \\\n')
    """
    ginfinity-generate-node-embeddings \
${cliBlock}

    python3 - <<'PY'
import ast
import pandas as pd
import sys

ID_COL = '${params.id_column}'
POSSIBLE_ORDER_COLS = ['node_index', 'position', 'idx', 'node_idx', 'node_position']

try:
    df = pd.read_csv('raw_node_embeddings.tsv', sep='\t')
except Exception as exc:
    sys.exit(f'Failed to read raw node embeddings: {exc}')

if ID_COL not in df.columns:
    sys.exit(f'Missing expected column {ID_COL} in ginfinity output')

if 'node_embeddings' in df.columns:
    exploded_records = []
    for _, row in df.iterrows():
        raw_embeddings = row['node_embeddings']
        try:
            parsed_embeddings = ast.literal_eval(str(raw_embeddings))
        except Exception as exc:
            sys.exit(f"Failed to parse node_embeddings for {row[ID_COL]}: {exc}")
        if not isinstance(parsed_embeddings, (list, tuple)):
            sys.exit(f"node_embeddings entry for {row[ID_COL]} is not list-like")
        for node_idx, vector in enumerate(parsed_embeddings):
            if not isinstance(vector, (list, tuple)):
                sys.exit(f"Embedding {node_idx} for {row[ID_COL]} is not list-like")
            try:
                values = [float(val) for val in vector]
            except Exception as exc:
                sys.exit(f"Failed to coerce embedding {node_idx} for {row[ID_COL]} to floats: {exc}")
            exploded_records.append({
                ID_COL: str(row[ID_COL]),
                'node_index': int(node_idx),
                'embedding_vector': ','.join(str(val) for val in values),
            })
    if not exploded_records:
        sys.exit('node_embeddings column present but no vectors extracted')
    pd.DataFrame(exploded_records).to_csv('node_embeddings.tsv', sep='\t', index=False)
    sys.exit(0)

order_col = None
for candidate in POSSIBLE_ORDER_COLS:
    if candidate in df.columns:
        order_col = candidate
        break
if order_col is None:
    sys.exit('Unable to locate node index column in ginfinity output')

def infer_embedding_column(frame: pd.DataFrame) -> pd.Series:
    if 'embedding_vector' in frame.columns:
        return frame['embedding_vector'].astype(str)
    embed_cols = [c for c in frame.columns if c.startswith('emb') or c.startswith('gin') or c.startswith('feat_')]
    if not embed_cols:
        # fallback: assume all numeric columns beyond ID/position belong to the embedding
        exclude = {ID_COL, order_col}
        embed_cols = [c for c in frame.columns if c not in exclude]
    if not embed_cols:
        sys.exit('Could not infer embedding columns in ginfinity output')
    return frame[embed_cols].astype(float).astype(str).agg(','.join, axis=1)

vec_series = infer_embedding_column(df)

out_df = pd.DataFrame({
    ID_COL: df[ID_COL].astype(str),
    'node_index': df[order_col].astype(int),
    'embedding_vector': vec_series,
})
out_df.to_csv('node_embeddings.tsv', sep='\t', index=False)
PY
    """
}
