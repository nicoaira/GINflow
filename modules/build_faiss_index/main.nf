process BUILD_FAISS_INDEX {
    tag "faiss"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/405599a0a76073fa7c942f2773cad66490aabb99fb9825c79ba50ad1157e9a6a/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_numpy_faiss-gpu:7382ed4d7c6e6258' : 'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e') }"

    input:
    path windows, stageAs: 'windows/*'
    path manifests, stageAs: 'manifests/*'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'
    path quantization, stageAs: 'quantization'
    val n_windows

    output:
    path "index", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def gpu_flag  = task.accelerator ? '--gpu' : ''
    def nlist     = params.faiss_nlist     != null ? "--nlist ${params.faiss_nlist}" : ''
    def nprobe    = params.faiss_nprobe    != null ? "--nprobe ${params.faiss_nprobe}" : ''
    def quantization_arg = params.quantize == 'sq' ? '--quantization quantization' : ''
    """
    build_faiss.py \\
        --windows windows/*.windows.npz \\
        --manifests manifests/*.windows.manifest.json \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir index \\
        --backend faiss \\
        --index-type ${params.faiss_index} \\
        --hnsw-m ${params.faiss_hnsw_m} \\
        --hnsw-ef-construction ${params.faiss_hnsw_ef_construction} \\
        --hnsw-ef-search ${params.faiss_hnsw_ef_search} \\
        ${nlist} \\
        ${nprobe} \\
        ${quantization_arg} \\
        ${gpu_flag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p index
    touch index/index.faiss
    touch index/windows.tsv
    touch index/embeddings.npz
    touch index/records.tsv
    if [ "${params.quantize}" = "sq" ]; then
        mkdir -p index/quantization
        cp -a quantization/. index/quantization/
    fi
    echo '{}' > index/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faiss: \$(python3 -c "import faiss; print(faiss.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
