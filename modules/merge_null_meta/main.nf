#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_NULL_META {
    tag { "merge_null_meta_batch_${ids.size()}" }

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-2.2.2:latest' :
        'amancevice/pandas:2.2.2' }"

    input:
    tuple val(ids), path(meta_tsvs)

    output:
    tuple val(ids), path('batch.tsv'), emit: batch

    script:
    """
    python3 - << 'PY'
import pandas as pd
files = ${groovy.json.JsonOutput.toJson(meta_tsvs.collect { it.toString() })}
# Defensive: ensure deterministic order regardless of upstream ordering
files = sorted(files)
dfs = [pd.read_csv(f, sep='\t', dtype=str) for f in files]
pd.concat(dfs, ignore_index=True).to_csv('batch.tsv', sep='\t', index=False)
print(f"Merged {len(files)} files -> batch.tsv")
PY
    """
}
