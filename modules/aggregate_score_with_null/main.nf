#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AGGREGATE_SCORE_WITH_NULL {
    tag { "aggregate_score_with_null_${query_id}" }

    label 'lightweight'

    publishDir "${params.outdir}/queries_results/${query_id}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-aggregate-score:latest' :
        'nicoaira/ginflow-aggregate-score:latest' }"

    input:
    tuple val(query_id), path(sorted_distances), path(null_scores)
    path  meta_map

    output:
    tuple val(query_id), path("pairs_scores_all_contigs.tsv"),            emit: contigs
    tuple val(query_id), path("pairs_scores_all_contigs.windows.tsv"),    emit: windows
    tuple val(query_id), path("pairs_scores_all_contigs.aggregated.tsv"), emit: enriched_agg

    script:
    """
    # 1) Aggregate scores from sorted distances
    python3 ${baseDir}/bin/aggregated_score.py \
      --input ${sorted_distances} \
      --id-column ${params.id_column} \
      --alpha1 ${params.alpha1} --alpha2 ${params.alpha2} \
      --beta1 ${params.beta1}   --beta2 ${params.beta2} \
      --gamma ${params.gamma}   --percentile ${params.percentile} \
      --mode contigs \
      --output raw_contigs.tsv \
      --output-unaggregated

    # 2) Enrich with meta map like the standard AGGREGATE_SCORE
    python3 - << 'PY'
import pandas as pd
idc      = '${params.id_column}'

meta     = pd.read_csv('id_meta.tsv', sep='\t', dtype=str)
meta_cols = [c for c in meta.columns if c != idc]

m1 = meta.rename(columns={idc:'query_id',   **{c:f'query_{c}'   for c in meta_cols}})
m2 = meta.rename(columns={idc:'subject_id', **{c:f'subject_{c}' for c in meta_cols}})

all_df   = pd.read_csv('raw_contigs.tsv', sep='\t').merge(m1,on='query_id').merge(m2,on='subject_id')
unagg_df = pd.read_csv('raw_contigs.unaggregated.tsv', sep='\t').merge(m1,on='query_id').merge(m2,on='subject_id')

all_df.to_csv('pairs_scores_all_contigs.tsv', sep='\t', index=False)
unagg_df.to_csv('pairs_scores_all_contigs.windows.tsv', sep='\t', index=False)
PY

    python3 ${baseDir}/bin/aggregate_structural_contigs.py \
      --input pairs_scores_all_contigs.tsv \
      --output pairs_scores_all_contigs.aggregated.tsv \
      --max-contig-overlap ${params.max_contig_overlap} \
      --structure-column-name ${params.structure_column_name}

    # 3) Insert z_score, p_value and q_value (FDR) next to 'score' using the provided null distribution
    python3 - << 'PY'
import pandas as pd, math

def norm_sf(z: float) -> float:
    # Numerically stable survival function for Normal(0,1)
    return 0.5 * math.erfc(z / math.sqrt(2.0))

scores = pd.read_csv('pairs_scores_all_contigs.aggregated.tsv', sep='\t')
null   = pd.read_csv('null_scores.tsv', sep='\t')

if 'score' not in scores.columns:
    raise SystemExit("'pairs_scores_all_contigs.tsv' missing 'score' column")
if 'score' not in null.columns:
    raise SystemExit("null distribution missing 'score' column")

mu = null['score'].mean()
sd = null['score'].std(ddof=0)
if sd == 0:
    # avoid div-by-zero; set z=0 and p=0.5 (non-informative)
    z = pd.Series([0.0] * len(scores))
    p = pd.Series([0.5] * len(scores))
else:
    z = (scores['score'] - mu) / sd
    p = z.map(norm_sf)

# Benjamini–Hochberg FDR adjustment
def fdr_bh(pvals: pd.Series) -> pd.Series:
    s = pvals.fillna(1.0)
    m = len(s)
    if m == 0:
        return s
    sorted_idx = s.sort_values(kind='mergesort').index
    sorted_p = s.loc[sorted_idx].tolist()
    adj = [min(p * m / (i + 1), 1.0) for i, p in enumerate(sorted_p)]
    for i in range(m - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    q_map = {idx: val for idx, val in zip(sorted_idx, adj)}
    return s.index.to_series().map(q_map)

# insert columns to the right of 'score'
cols = scores.columns.tolist()
try:
    idx = cols.index('score')
except ValueError:
    idx = len(cols) - 1
scores.insert(idx+1, 'z_score', z)
scores.insert(idx+2, 'p_value', p)
scores.insert(idx+3, 'q_value', fdr_bh(p))

scores.to_csv('pairs_scores_all_contigs.aggregated.tsv', sep='\t', index=False)
PY
    """
}
