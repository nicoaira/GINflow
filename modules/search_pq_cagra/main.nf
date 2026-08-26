process SEARCH_PQ_CAGRA {
    tag "${meta.id}"
    label 'process_medium'
    accelerator { (params.search_device == 'cpu' || params.cagra_to_hnsw) ? null : 1 }

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'oras://community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:e0102895f642d35a' : 'oras://community.wave.seqera.io/library/python_numpy_pq-cagra-adc-cpu:230079361cb63119') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:ed61007abe908991' : 'community.wave.seqera.io/library/python_numpy_pq-cagra-adc-cpu:2ea6d576a29716e1') }"

    input:
    tuple val(meta), path(windows), path(manifest)
    path database

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def device = (params.search_device == 'cpu' || params.cagra_to_hnsw) ? 'cpu' : 'cuda'
    def separate_rerank = BooleanParam.rerankEnabled(params.exact_rerank, params.hnswlib_rerank)
    def search_k = separate_rerank ? params.candidate_k : params.seed_k
    def min_sim = BooleanParam.annMinSimilarity(params.exact_rerank, params.hnswlib_rerank, params.quantize, params.seed_min_similarity)
    """
    mkdir -p query_windows
    cp ${windows} query_windows/
    cp ${manifest} query_windows/

    search_pq_cagra.py \\
        --windows-dir query_windows \\
        --database ${database} \\
        --output ${prefix}.seeds.tsv \\
        --k ${search_k} \\
        --min-similarity=${min_sim} \\
        --device ${device} \\
        --itopk-size ${params.cagra_itopk_size} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: stub
    END_VERSIONS
    """
}
