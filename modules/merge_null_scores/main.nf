#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_NULL_SCORES {
    tag { "merge_null_scores_${query_id}" }

    label 'lightweight'

    input:
    // Receive the list of null score paths as plain values; we read them from host
    // to avoid container volume issues.
    tuple val(query_id), val(score_files)

    output:
    tuple val(query_id), path('null_scores.tsv'), emit: null_scores

    script:
    def files_json = groovy.json.JsonOutput.toJson(score_files.collect { it.toString() })

    """
   python3 - << 'PY'
import csv
from pathlib import Path

files = sorted(Path(p) for p in ${files_json})
if not files:
    raise SystemExit('No null score files provided')

header = None
with Path('null_scores.tsv').open('w', encoding='utf-8', newline='') as out_fh:
    writer = None
    for idx, path in enumerate(files):
        if not path.exists():
            raise FileNotFoundError(f"Missing null score file: {path}")
        with path.open('r', encoding='utf-8') as in_fh:
            reader = csv.reader(in_fh, delimiter='\t')
            try:
                file_header = next(reader)
            except StopIteration:
                continue  # skip empty files
            if header is None:
                header = file_header
                writer = csv.writer(out_fh, delimiter='\t', lineterminator=chr(10))
                writer.writerow(header)
            elif file_header != header:
                raise ValueError(f"Header mismatch in {path}: {file_header} != {header}")
            for row in reader:
                if row:
                    writer.writerow(row)

print(f"Merged {len(files)} null score files -> null_scores.tsv")
PY
    """
}
