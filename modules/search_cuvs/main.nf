process SEARCH_CUVS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd' :
        'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4' }"

    input:
    tuple val(meta), path(windows), path(manifest), val(query_stats)
    path database
    val database_meta

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def gpu_flag = task.accelerator ? '--gpu' : ''
    def n_probes = params.cuvs_n_probes != null ? "--cuvs-n-probes ${params.cuvs_n_probes}" : ''
    def separate_rerank = params.exact_rerank || params.hnswlib_rerank
    def search_k = separate_rerank ? params.candidate_k : params.seed_k
    def min_sim = (separate_rerank || params.quantize.toString().toLowerCase() in ['pq', 'opq']) ? '-inf' : String.valueOf(params.seed_min_similarity)
    def search_device = task.accelerator ? 'gpu' : 'cpu'
    """
    search_faiss.py \\
        --windows ${windows} \\
        --manifest ${manifest} \\
        --database ${database} \\
        --output ${prefix}.seeds.tsv \\
        --k ${search_k} \\
        --min-similarity=${min_sim} \\
        --search-device ${search_device} \\
        ${n_probes} \\
        ${gpu_flag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuvs: \$(python3 -c "from importlib.metadata import version; print(version('cuvs'))")
        cupy: \$(python3 -c "import cupy; print(cupy.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuvs: stub
        cupy: stub
        numpy: stub
    END_VERSIONS
    """
}
