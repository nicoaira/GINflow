#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_NULL_SCORES {
    tag { "merge_null_scores_${query_id}" }

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-2.2.2:latest' :
        'amancevice/pandas:2.2.2' }"

    input:
    // Use val for score_files to avoid file name collisions when staging
    tuple val(query_id), val(score_files)

    output:
    tuple val(query_id), path('null_scores.tsv'), emit: null_scores

    script:
    """
    python3 - << 'PY'
import pandas as pd
files = ${groovy.json.JsonOutput.toJson(score_files.collect { it.toString() })}
# Defensive: keep deterministic file order
files = sorted(files)
dfs = [pd.read_csv(f, sep='\t') for f in files]
pd.concat(dfs, ignore_index=True).to_csv('null_scores.tsv', sep='\t', index=False)
print(f"Merged {len(files)} null score files -> null_scores.tsv")
PY
    """
}
