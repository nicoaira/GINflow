#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PLOT_DISTANCES {
    tag "plot_distances_${query_id}"

    label 'lightweight'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/nicoaira/ginflow-plot-distances:latest' :
        'nicoaira/ginflow-plot-distances:latest' }"

    when   { params.plot_distances_distribution }
    publishDir "${params.outdir}/queries_results/${query_id}/plots", mode: 'copy'

    input:
    tuple val(query_id), path(sorted_distances)

    output:
    path "distance_distribution.png"

    script:
    """
    python3 - << 'PY'
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv('distances.sorted.tsv', sep='\t')
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
