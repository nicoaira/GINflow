#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ───────────────────────────────────────────────────────────
// User-tunable parameters (override on CLI)
// ───────────────────────────────────────────────────────────
params.input                 = "$baseDir/input.tsv"           // TSV with sequences
params.outdir                = "results"
params.model_path            = "$baseDir/model.pth"
params.structure_column_name = 'secondary_structure'
params.structure_column_num  = null
params.id_column             = 'transcript_id'

params.device                = 'cpu'
params.num_workers           = 4
params.subgraphs             = false
params.L                     = null
params.keep_paired_neighbors = false
params.retries               = 0
params.batch_size_embed      = 100       // sequences per embedding task

params.faiss_k               = 1000      // neighbors to retrieve in FAISS

// aggregation parameters
params.alpha1    = 0.25
params.alpha2    = 0.24
params.beta1     = 0.0057
params.beta2     = 1.15
params.gamma     = 0.41
params.percentile= 1
params.top_n     = 10

// plotting flags
params.plot_distances_distribution = true
params.hist_seed                   = 42
params.hist_frac                   = 0.001
params.hist_bins                   = 200

params.plot_score_distribution     = true
params.score_bins                  = 30

// ───────────────────────────────────────────────────────────
// Basic checks
// ───────────────────────────────────────────────────────────
if ( !file(params.input).exists() )
    error "Cannot find input TSV: ${params.input}"

// ───────────────────────────────────────────────────────────
// NEW: Generate windows if requested
// ───────────────────────────────────────────────────────────
process GENERATE_WINDOWS {
    tag "generate_windows"
    cpus 14

    // initialize BLAS thread limits before running the script
    beforeScript 'export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1'

    input:
    path orig_tsv

    output:
    tuple path("windows_graphs.pt"), path("windows_metadata.tsv"), emit: window_files

    script:
    """
    python3 ${baseDir}/modules/generate_windows.py \\
      --input ${orig_tsv} \\
      --output-dir . \\
      --id-column ${params.id_column} \\
      --structure-column-name ${params.structure_column_name} \\
      --L ${params.L} \\
      ${params.keep_paired_neighbors ? '--keep-paired-neighbors' : ''} \\
      --mask-threshold ${params.mask_threshold} \\
      --num-workers ${task.cpus}
    """
}


// ───────────────────────────────────────────────────────────
// 0)  Extract lightweight metadata (ID → sequences + structure)
//     Keeps bulky cols out of embeddings.tsv
// ───────────────────────────────────────────────────────────
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
df  = pd.read_csv('${orig_tsv}', sep='\\t', dtype=str)
idc = '${params.id_column}'
scol = '${params.structure_column_name}'.strip()
seq_cols = [c for c in df.columns if 'sequence' in c.lower()]
if scol and scol in df.columns and scol not in seq_cols:
    seq_cols.append(scol)
if not seq_cols:
    sys.exit('No sequence or structure columns found – expected at least one column containing \"sequence\" or the specified structure column')
df[[idc] + seq_cols].to_csv('id_meta.tsv', sep='\\t', index=False)
PY
    """
}

// ───────────────────────────────────────────────────────────
// Prepare transcript batches as TSV files
// ───────────────────────────────────────────────────────────
process PREP_BATCH {
    tag "prep_batch"
    input:
    val rows
    output:
    path "batch_${task.index}.tsv", emit: batch_file
    script:
    def header = rows[0].keySet().join("\t")
    def lines = rows.collect { r -> r.values().join("\t") }.join("\n")
    """
    cat > batch_${task.index}.tsv << 'EOF'
${header}
${lines}
EOF
    """
}

// 1) Generate embeddings per batch
process GENERATE_EMBEDDINGS {
    tag { "embeddings batch=${task.index}" }
    maxForks = 1

    input:
    each item // either a batch_tsv path or a tuple [graphs_pt, metadata_tsv]

    output:
    path "embeddings_batch_${task.index}.tsv", emit: batch_embeddings

    script:
    def common_args = """ \\
      --model-path ${params.model_path} \\
      --id-column ${params.id_column} \\
      --output embeddings_batch_${task.index}.tsv \\
      --device ${params.device} \\
      --num-workers ${params.num_workers} \\
      --batch-size ${params.batch_size}
    """

    def cmd
    if (params.subgraphs) {
        // 'item' is [graphs_pt, metadata_tsv]
        def graphs_pt_file = item[0]
        def metadata_tsv_file = item[1]
        cmd = """
        python3 ${baseDir}/modules/generate_embeddings.py \\
          --graph-pt ${graphs_pt_file} \\
          --meta-tsv ${metadata_tsv_file} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          ${common_args}
        """
    } else {
        // 'item' is batch_tsv_file (path)
        def batch_tsv_file = item
        cmd = """
        python3 ${baseDir}/modules/generate_embeddings.py \\
          --input ${batch_tsv_file} \\
          --structure-column-name ${params.structure_column_name} \\
          --keep-cols ${params.id_column},window_start,window_end \\
          ${common_args}
        """
    }
    """
    ${cmd}
    """
}

// 2) Merge per-batch embeddings
process MERGE_EMBEDDINGS {
    tag "merge_embeddings"
    container ''
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path batch_embeddings

    output:
    path "embeddings.tsv", emit: embeddings

    script:
    """
    head -n1 ${batch_embeddings[0]} > embeddings.tsv
    for f in ${batch_embeddings.join(' ')}; do
      tail -n +2 \$f >> embeddings.tsv
    done
    """
}

// 3) Build FAISS index
process BUILD_FAISS_INDEX {
    tag    "build_faiss_index"
    cpus   1
    memory '16 GB'
    publishDir "${params.outdir}/faiss_index", mode: 'copy'

    input:
    path embeddings

    output:
    path "faiss.index",       emit: faiss_idx
    path "faiss_mapping.tsv", emit: faiss_map

    script:
    """
    python3 ${baseDir}/modules/build_faiss_index.py \
      --input embeddings.tsv \
      --id-column ${params.id_column} \
      --query ${params.query} \
      --index-path faiss.index \
      --mapping-path faiss_mapping.tsv
    """
}

// 4) Query FAISS index
process QUERY_FAISS_INDEX {
    tag  "query_faiss_index"
    cpus params.num_workers

    input:
    path embeddings
    path faiss_idx
    path faiss_map

    output:
    path "distances.tsv", emit: distances

    script:
    """
    python3 ${baseDir}/modules/query_faiss_index.py \
      --input embeddings.tsv \
      --id-column ${params.id_column} \
      --query ${params.query} \
      --index-path faiss.index \
      --mapping-path faiss_mapping.tsv \
      --top-k ${params.faiss_k} \
      --output distances.tsv
    """
}

// 5) Sort distances
process SORT_DISTANCES {
    tag "sort_distances"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path distances

    output:
    path "distances.sorted.tsv", emit: sorted_distances

    script:
    """
    python3 - << 'PY'
import pandas as pd
df = pd.read_csv('distances.tsv', sep='\\t')
df.sort_values('distance', inplace=True)
df.to_csv('distances.sorted.tsv', sep='\\t', index=False)
PY
    """
}

// 6) Plot distance distribution
process PLOT_DISTANCES {
    when   { params.plot_distances_distribution }
    tag    "plot_distances"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path sorted_distances

    output:
    path "distance_distribution.png"

    script:
    """
    python3 - << 'PY'
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv('distances.sorted.tsv', sep='\\t')
sample = df.sample(frac=${params.hist_frac}, random_state=${params.hist_seed})
plt.hist(sample['distance'], bins=${params.hist_bins})
plt.xlabel('Distance')
plt.ylabel('Frequency')
plt.title('Distance Distribution')
plt.tight_layout()
plt.savefig('distance_distribution.png')
PY
    """
}

// 7) Aggregate and enrich (re-attach sequences *and* secondary structures)
process AGGREGATE_SCORE {
    tag "aggregate_score"
    cpus params.num_workers
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sorted_distances
    path meta_map

    output:
    path "pairs_scores_all_contigs.tsv",              emit: enriched_all
    path "pairs_scores_all_contigs.unaggregated.tsv", emit: enriched_unagg

    script:
    """
    python3 ${baseDir}/modules/aggregated_score.py \
      --input distances.sorted.tsv \
      --id-column ${params.id_column} \
      --alpha1 ${params.alpha1} --alpha2 ${params.alpha2} \
      --beta1 ${params.beta1}   --beta2 ${params.beta2} \
      --gamma ${params.gamma}   --percentile ${params.percentile} \
      --mode contigs \
      --output raw_contigs.tsv \
      --output-unaggregated

    python3 - << 'PY'
import pandas as pd
idc      = '${params.id_column}'

# meta_map has IDs + *all* sequence and structure columns
meta     = pd.read_csv('id_meta.tsv', sep='\\t', dtype=str)
meta_cols = [c for c in meta.columns if c != idc]

m1 = meta.rename(columns={idc:f'{idc}_1', **{c:f'{c}_1' for c in meta_cols}})
m2 = meta.rename(columns={idc:f'{idc}_2', **{c:f'{c}_2' for c in meta_cols}})

all_df   = pd.read_csv('raw_contigs.tsv', sep='\\t').merge(m1,on=f'{idc}_1').merge(m2,on=f'{idc}_2')
unagg_df = pd.read_csv('raw_contigs.unaggregated.tsv', sep='\\t').merge(m1,on=f'{idc}_1').merge(m2,on=f'{idc}_2')

all_df.to_csv('pairs_scores_all_contigs.tsv', sep='\\t', index=False)
unagg_df.to_csv('pairs_scores_all_contigs.unaggregated.tsv', sep='\\t', index=False)
PY
    """
}

// 8) Plot score distribution
process PLOT_SCORE {
    when   { params.plot_score_distribution }
    tag    "plot_score"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path enriched_all

    output:
    path "score_distribution.png"

    script:
    """
    python3 - << 'PY'
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv('pairs_scores_all_contigs.tsv', sep='\\t')
plt.hist(df['score'], bins=${params.score_bins})
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.title('Contig Score Distribution')
plt.tight_layout()
plt.savefig('score_distribution.png')
PY
    """
}

// 9) Filter top-N contigs
process FILTER_TOP_CONTIGS {
    tag "filter_top_contigs"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path enriched_all
    path enriched_unagg

    output:
    path "pairs_scores_top_contigs.tsv",              emit: top_contigs
    path "pairs_scores_top_contigs.unaggregated.tsv", emit: top_contigs_unagg

    script:
    """
    python3 - << 'PY'
import pandas as pd
top = ${params.top_n}
all_df = pd.read_csv('pairs_scores_all_contigs.tsv', sep='\\t')
un_df  = pd.read_csv('pairs_scores_all_contigs.unaggregated.tsv', sep='\\t')
sel    = all_df[all_df.contig_rank <= top]
ids    = sel.contig_id.unique()
sel.to_csv('pairs_scores_top_contigs.tsv', sep='\\t', index=False)
un_df[un_df.contig_id.isin(ids)].to_csv('pairs_scores_top_contigs.unaggregated.tsv', sep='\\t', index=False)
PY
    """
}

// 10) Draw contig-level SVGs/PNGs
process DRAW_CONTIG_SVGS {
    tag "draw_contig_svgs"
    cpus params.num_workers
    publishDir "${params.outdir}/drawings/contigs", mode: 'copy'

    input:
    path top_contigs_tsv

    output:
    path 'individual_svgs', emit: contig_individual

    script:
    """
    mkdir -p individual_svgs
    python3 /app/draw_pairs.py \
      --tsv ${top_contigs_tsv} --outdir individual_svgs \
      --width 500 --height 500 --highlight-colour "#00FF99" \
      --num-workers ${params.num_workers}
    """
}

// 11) Draw window-level SVGs/PNGs
process DRAW_UNAGG_SVGS {
    tag "draw_window_svgs"
    cpus params.num_workers
    publishDir "${params.outdir}/drawings/unagg_windows", mode: 'copy'

    input:
    path top_unagg_tsv

    output:
    path 'individual_svgs', emit: window_individual

    script:
    """
    mkdir -p individual_svgs
    python3 /app/draw_pairs.py \
      --tsv ${top_unagg_tsv} --outdir individual_svgs \
      --width 500 --height 500 --highlight-colour "#00FF99" \
      --num-workers ${params.num_workers}
    """
}

// 12) Generate contig-level HTML report
process GENERATE_AGGREGATED_REPORT {
    tag "gen_agg_report"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path top_contigs_tsv
    path contig_individual

    output:
    path "pairs_contigs_report.html"

    script:
    """
    python3 /app/generate_report.py \
      --pairs ${top_contigs_tsv} \
      --svg-dir ${contig_individual} \
      --id-column ${params.id_column} \
      --output pairs_contigs_report.html
    """
}

// 13) Generate window-level HTML report
process GENERATE_UNAGGREGATED_REPORT {
    tag "gen_unagg_report"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path top_contigs_unagg_tsv
    path window_individual

    output:
    path "pairs_contigs_report.unaggregated.html"

    script:
    """
    python3 /app/generate_report.py \
      --pairs ${top_contigs_unagg_tsv} \
      --svg-dir ${window_individual} \
      --id-column ${params.id_column} \
      --output pairs_contigs_report.unaggregated.html
    """
}

// ───────────────────────────────────────────────────────────
// Workflow wiring
// ───────────────────────────────────────────────────────────
workflow {
    // 0) Build the metadata map once
    def meta = EXTRACT_META_MAP(file(params.input))

    // Split transcripts into batches
    // Ensuring header and separator are correctly used for TSV
    def transcript_batches = Channel.fromPath(params.input)
        .splitCsv(header: params.header, sep:'\t', strip:true, by: params.batch_size_embed)
    def batch_files_ch = PREP_BATCH(transcript_batches).batch_file

    def gen
    if(params.subgraphs) {
        def window_files = GENERATE_WINDOWS(batch_files_ch)
        gen = GENERATE_EMBEDDINGS(window_files)
    } else {
        gen = GENERATE_EMBEDDINGS(batch_files_ch)
    }
    def merged = MERGE_EMBEDDINGS(gen.batch_embeddings)
    def idx     = BUILD_FAISS_INDEX(merged.embeddings)
    def dists   = QUERY_FAISS_INDEX(merged.embeddings, idx.faiss_idx, idx.faiss_map)
    def sorted  = SORT_DISTANCES(dists)
    PLOT_DISTANCES(sorted)

    // 7
    def (agg_all, agg_un) = AGGREGATE_SCORE(sorted, meta.meta_map)
    PLOT_SCORE(agg_all)

    // 8-9-10-11-12-13
    def (top_c, top_u)  = FILTER_TOP_CONTIGS(agg_all, agg_un)
    def contigs_draw    = DRAW_CONTIG_SVGS(top_c)
    def windows_draw    = DRAW_UNAGG_SVGS(top_u)
    GENERATE_AGGREGATED_REPORT(top_c, contigs_draw.contig_individual)
    GENERATE_UNAGGREGATED_REPORT(top_u, windows_draw.window_individual)
}
