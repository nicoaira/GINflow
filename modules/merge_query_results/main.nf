#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_QUERY_RESULTS {
    tag "merge_query_results"

    label 'lightweight'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    // Collect all per-query result files as values. Using `val` avoids Nextflow
    // staging (and renaming) the files in the work directory, which previously
    // produced bogus temporary names like `Script_***` and caused the process to
    // fail.  The files remain accessible at their original locations and are
    // read directly by the Python merging script below.
    val distances
    val scores_contigs
    val scores_windows
    val scores_agg

    output:
    path "distances.merged.sorted.tsv",                emit: distances
    path "pairs_scores_all_contigs.merged.tsv",        emit: scores
    path "pairs_scores_all_contigs.windows.merged.tsv", emit: scores_win
    path "pairs_scores_all_contigs.aggregated.merged.tsv", emit: scores_agg

    script:
    def merge_plan = [
        [distances.collect { it.toString() }, 'distances.merged.sorted.tsv'],
        [scores_contigs.collect { it.toString() }, 'pairs_scores_all_contigs.merged.tsv'],
        [scores_windows.collect { it.toString() }, 'pairs_scores_all_contigs.windows.merged.tsv'],
        [scores_agg.collect { it.toString() }, 'pairs_scores_all_contigs.aggregated.merged.tsv']
    ]
    def merge_json = groovy.json.JsonOutput.toJson(merge_plan)

    """
    python3 - <<'PY'
import csv
import json
from pathlib import Path

merge_plan = json.loads('''${merge_json}''')

def merge(files, out_path):
    files = [Path(p) for p in sorted(files)]
    if not files:
        return

    header = None
    with Path(out_path).open('w', encoding='utf-8', newline='') as out_fh:
        writer = None
        for path in files:
            if not path.exists():
                raise FileNotFoundError(f"Missing result file: {path}")
            with path.open('r', encoding='utf-8') as in_fh:
                reader = csv.reader(in_fh, delimiter='\t')
                try:
                    file_header = next(reader)
                except StopIteration:
                    continue
                if header is None:
                    header = file_header
                    writer = csv.writer(out_fh, delimiter='\t', lineterminator=chr(10))
                    writer.writerow(header)
                elif file_header != header:
                    raise ValueError(f"Header mismatch in {path}: {file_header} != {header}")
                for row in reader:
                    if row:
                        writer.writerow(row)

for files, output in merge_plan:
    merge(files, output)
PY
    """
}
