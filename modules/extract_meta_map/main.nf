#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process EXTRACT_META_MAP {
    tag "extract_meta_map"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path orig_tsv

    output:
    path "id_meta.tsv", emit: meta_map

    script:
    """
    python3 - << 'PY'
import pandas as pd, sys
df = pd.read_csv('${orig_tsv}', sep='\t', dtype=str)
idc = '${params.id_column}'
# determine structure column: by name or index
scol_name = '${params.structure_column_name}'
scol_num = ${params.structure_column_num ?: 'None'}
if scol_name:
    scol = scol_name
elif scol_num is not None:
    try:
        scol = df.columns[int(scol_num)]
    except Exception:
        sys.exit(f'Specified structure_column_num {scol_num} out of range')
else:
    scol = ''
# detect sequence-like columns
seq_cols = [c for c in df.columns if 'sequence' in c.lower()]
if scol and scol in df.columns and scol not in seq_cols:
    seq_cols.append(scol)
if not seq_cols:
    sys.exit('No sequence columns found – expected a column containing "sequence"')
df[[idc] + seq_cols].to_csv('id_meta.tsv', sep='\\t', index=False)
PY
    """
}