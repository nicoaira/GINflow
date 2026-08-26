process GENERATE_QUANTIZED_WINDOWS {
    tag "quantized-windows"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path quantization, stageAs: 'quantization'
    val quantized_window_stats

    output:
    path "quantized_windows", emit: windows
    path "quantized_windows/*.windows.npz", emit: npz
    path "quantized_windows/*.windows.manifest.json", emit: manifests
    path "windows.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    generate_quantized_windows.py \\
        --input-dir quantization \\
        --output-dir quantized_windows \\
        --window-size ${params.window_size} \\
        --stride ${params.window_stride} \\
        ${args}

    cp quantized_windows/windows.json windows.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p quantized_windows
    python3 - <<'PY'
import numpy as np, json
from pathlib import Path
Path('quantized_windows').mkdir(exist_ok=True)
np.savez_compressed('quantized_windows/stub.windows.npz', stub=np.zeros((1, 8), dtype=np.uint8))
Path('quantized_windows/stub.windows.manifest.json').write_text('{"records":[],"window_size":11,"stride":1}')
Path('quantized_windows/windows.json').write_text('{}')
Path('windows.json').write_text('{"n_windows":0,"records":0,"window_size":11,"stride":1}')
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
