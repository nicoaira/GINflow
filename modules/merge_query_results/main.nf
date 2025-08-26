#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_QUERY_RESULTS {
    tag "merge_query_results"

    label 'lightweight'

    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/amancevice-pandas-2.2.2:latest' :
        'amancevice/pandas:2.2.2' }"

    input:
    // Collect all per-query result files as values. Using `val` avoids Nextflow
    // staging (and renaming) the files in the work directory, which previously
    // produced bogus temporary names like `Script_***` and caused the process to
    // fail.  The files remain accessible at their original locations and are
    // read directly by the Python merging script below.
    val distances
    val scores_all
    val scores_unagg

    output:
    path "distances.merged.sorted.tsv"
    path "pairs_scores_all_contigs.merged.tsv"
    path "pairs_scores_all_contigs.unaggregated.merged.tsv"

    script:
    """
    python3 - <<'PY'
import pandas as pd

dist_files = ${groovy.json.JsonOutput.toJson(distances.collect { it.toString() })}
score_files = ${groovy.json.JsonOutput.toJson(scores_all.collect { it.toString() })}
unagg_files = ${groovy.json.JsonOutput.toJson(scores_unagg.collect { it.toString() })}

def merge(files, out):
    dfs = [pd.read_csv(f, sep='\t') for f in files]
    pd.concat(dfs, ignore_index=True).to_csv(out, sep='\t', index=False)

merge(dist_files, 'distances.merged.sorted.tsv')
merge(score_files, 'pairs_scores_all_contigs.merged.tsv')
merge(unagg_files, 'pairs_scores_all_contigs.unaggregated.merged.tsv')
PY
    """
}
