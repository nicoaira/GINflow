#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_NULL_WINDOWS {
    tag { "generate_null_windows_batch_${ids.size()}" }

    label 'high_cpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginfinity:latest' :
        'nicoaira/ginfinity' }"

    // initialize BLAS thread limits before running the script
    beforeScript 'export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1'

    input:
    tuple val(ids), path(meta_tsv)

    output:
    tuple val(ids), path("windows_graphs.pt"), path("windows_metadata.tsv"), emit: window_files_null

    script:
    """
    ginfinity-generate-windows \\
      --input ${meta_tsv} \\
      --output-dir . \\
      --id-column ${params.id_column} \\
      --structure-column-name ${params.structure_column_name} \\
      --L ${params.L} \\
      ${params.keep_paired_neighbors ? '--keep-paired-neighbors' : ''} \\
      --mask-threshold ${params.mask_threshold} \\
      --num-workers ${task.cpus}

    export GINFLOW_LOW_COMPLEXITY_MASK=${params.low_complexity_mask == null ? '' : params.low_complexity_mask}
    export GINFLOW_STRUCTURE_COLUMN_NAME="${params.structure_column_name ?: ''}"
    export GINFLOW_STRUCTURE_COLUMN_NUM="${params.structure_column_num ?: ''}"
    export GINFLOW_ID_COLUMN="${params.id_column}"
    export GINFLOW_INPUT_TSV="${meta_tsv}"

    python - <<'PY'
import os
import sys
import csv
from pathlib import Path
from collections import OrderedDict

try:
    import torch
except ImportError as exc:
    print(f"[low-complexity] Failed to import torch: {exc}", file=sys.stderr)
    sys.exit(1)

threshold_raw = (os.environ.get('GINFLOW_LOW_COMPLEXITY_MASK') or '').strip()
if not threshold_raw:
    sys.exit(0)
try:
    threshold = float(threshold_raw)
except ValueError:
    print(f"[low-complexity] Invalid threshold '{threshold_raw}', skipping masking", file=sys.stderr)
    sys.exit(0)

if threshold >= 1.0:
    sys.exit(0)
if threshold < 0.0:
    print(f"[low-complexity] Threshold {threshold} < 0.0; clamping to 0.0", file=sys.stderr)
    threshold = 0.0

meta_path = Path('windows_metadata.tsv')
graph_path = Path('windows_graphs.pt')

if not meta_path.exists() or meta_path.stat().st_size == 0:
    sys.exit(0)

with meta_path.open(newline='') as meta_file:
    meta_reader = csv.DictReader(meta_file, delimiter='\t')
    meta_rows = list(meta_reader)
    meta_fieldnames = meta_reader.fieldnames or []

if not meta_rows:
    sys.exit(0)

input_path = os.environ.get('GINFLOW_INPUT_TSV')
if not input_path:
    print('[low-complexity] Missing GINFLOW_INPUT_TSV environment variable', file=sys.stderr)
    sys.exit(1)

with open(input_path, newline='') as input_file:
    input_reader = csv.DictReader(input_file, delimiter='\t')
    input_fieldnames = input_reader.fieldnames or []
    input_rows = list(input_reader)

if not input_rows:
    sys.exit(0)

structure_name = (os.environ.get('GINFLOW_STRUCTURE_COLUMN_NAME') or '').strip()
structure_num_raw = (os.environ.get('GINFLOW_STRUCTURE_COLUMN_NUM') or '').strip()
id_column = (os.environ.get('GINFLOW_ID_COLUMN') or '').strip()

if not id_column:
    print('[low-complexity] Missing ID column name', file=sys.stderr)
    sys.exit(1)

structure_col = None
if structure_name and structure_name in input_fieldnames:
    structure_col = structure_name
elif structure_num_raw:
    try:
        idx = int(structure_num_raw)
        structure_col = input_fieldnames[idx]
    except (ValueError, IndexError):
        print(f"[low-complexity] Invalid structure column index '{structure_num_raw}'", file=sys.stderr)
        sys.exit(1)
elif structure_name:
    print(f"[low-complexity] Structure column '{structure_name}' not found; skipping masking", file=sys.stderr)
    sys.exit(0)

if not structure_col:
    print('[low-complexity] Unable to resolve structure column; skipping masking', file=sys.stderr)
    sys.exit(0)

structures = {}
for row in input_rows:
    tid = row.get(id_column)
    if tid is None:
        continue
    structures[tid] = row.get(structure_col)

keep_rows = []
keep_ids = []
removed = 0
total = len(meta_rows)

for row in meta_rows:
    tid = row.get(id_column)
    structure = structures.get(tid)
    frac = None
    if isinstance(structure, str):
        structure = structure.strip()
        if structure:
            try:
                start = int(float(row['window_start']))
                end = int(float(row['window_end']))
            except (KeyError, ValueError, TypeError):
                start = None
                end = None
            if start is not None and end is not None and end >= start:
                start = max(start, 0)
                end = min(end, len(structure) - 1)
                if start < len(structure) and end >= start:
                    segment = structure[start:end + 1]
                    if segment:
                        frac = segment.count('.') / len(segment)
    if frac is not None and frac > threshold:
        removed += 1
        continue
    keep_rows.append(row)
    keep_ids.append(row.get('window_id'))

if removed == 0:
    sys.exit(0)

print(f"[low-complexity] Masking {removed} of {total} windows (threshold={threshold})")

with meta_path.open('w', newline='') as meta_out:
    writer = csv.DictWriter(meta_out, fieldnames=meta_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(keep_rows)

if not graph_path.exists():
    print('[low-complexity] Missing windows_graphs.pt; cannot filter graphs', file=sys.stderr)
    sys.exit(1)

graphs_obj = torch.load(graph_path, weights_only=False)
keep_set = {wid for wid in keep_ids if wid}

if isinstance(graphs_obj, dict):
    filtered = OrderedDict()
    for row in meta_rows:
        wid = row.get('window_id')
        if wid and wid in keep_set:
            if wid not in graphs_obj:
                print(f"[low-complexity] Warning: graph missing for window {wid}", file=sys.stderr)
                continue
            filtered[wid] = graphs_obj[wid]
elif isinstance(graphs_obj, (list, tuple)):
    filtered_list = []
    for row, graph in zip(meta_rows, graphs_obj):
        wid = row.get('window_id')
        if wid and wid in keep_set:
            filtered_list.append(graph)
    filtered = type(graphs_obj)(filtered_list)
else:
    print(f"[low-complexity] Unsupported graph serialization type: {type(graphs_obj)}", file=sys.stderr)
    sys.exit(1)

torch.save(filtered, graph_path)
PY
    """.stripIndent()
}
