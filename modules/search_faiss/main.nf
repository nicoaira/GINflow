process SEARCH_FAISS {
    tag "${meta.id}"
    label 'process_low'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/405599a0a76073fa7c942f2773cad66490aabb99fb9825c79ba50ad1157e9a6a/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_numpy_faiss-gpu:7382ed4d7c6e6258' : 'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e') }"

    input:
    tuple val(meta), path(windows), path(manifest)
    path database

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def gpu_flag  = task.accelerator ? '--gpu' : ''
    def nprobe    = params.faiss_nprobe != null ? "--nprobe ${params.faiss_nprobe}" : ''
    def ef_search = params.faiss_hnsw_ef_search != null ? "--hnsw-ef-search ${params.faiss_hnsw_ef_search}" : ''
    def separate_rerank = params.exact_rerank || params.hnswlib_rerank
    def exact_index = params.faiss_index in ['flatip', 'flatl2']
    def search_k = (separate_rerank && !exact_index) ? params.candidate_k : params.seed_k
    def min_sim = (separate_rerank && !exact_index) || params.quantize.toString().toLowerCase() in ['pq', 'opq']
        ? '-inf'
        : String.valueOf(params.seed_min_similarity)
    """
    search_faiss.py \\
        --windows ${windows} \\
        --manifest ${manifest} \\
        --database ${database} \\
        --output ${prefix}.seeds.tsv \\
        --k ${search_k} \\
        --min-similarity=${min_sim} \\
        ${nprobe} \\
        ${ef_search} \\
        ${gpu_flag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
