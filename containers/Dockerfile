# ----------------------------------------------------------------------------
# containers/Dockerfile.gpu  (CUDA 12.4 / cu124)
# ----------------------------------------------------------------------------

# 1) Start from NVIDIA's CUDA 12.4 runtime image
FROM nvidia/cuda:12.4.0-runtime-ubuntu22.04

# 2) Install Python and other OS deps
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
       python3 python3-pip python3-venv git curl && \
    rm -rf /var/lib/apt/lists/*

# 3) Create and activate a virtualenv at /venv
RUN python3 -m venv /venv
ENV PATH="/venv/bin:$PATH"

# 4) Install PyTorch 2.6.0 + CUDA 12.4 (cu124) and torchvision
RUN pip install --no-cache-dir \
       "torch==2.6.0+cu124" "torchvision==0.21.0+cu124" \
       --extra-index-url https://download.pytorch.org/whl/cu124

# 5) Install torch-geometric + matching CUDA-extensions for cu124
RUN pip install --no-cache-dir torch-geometric && \
    pip install --no-cache-dir \
      pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
      -f https://data.pyg.org/whl/torch-2.6.0+cu124.html

# 6) Install any other Python deps (faiss-cpu, pandas, matplotlib, etc.) and ginfinity
RUN pip install --no-cache-dir \
      pandas matplotlib faiss-cpu \
      dask distributed \
      git+https://github.com/nicoaira/GINFINITY.git#egg=ginfinity

# 7) Copy your pipeline's code into /app
WORKDIR /app
COPY main.nf nextflow.config ./
COPY workflows/ workflows/
COPY bin/ bin/
COPY conf/ conf/
COPY presets/ presets/
RUN chmod +x bin/*.py

# 8) No ENTRYPOINT or CMD listed here.  Nextflow's Docker executor will run
#    "bash -ue .command.sh" inside /workspace when it spins up this image.
# ----------------------------------------------------------------------------
