FROM eclipse-temurin:17-jre-jammy AS base

# Install OS packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-venv git curl \
    && rm -rf /var/lib/apt/lists/*

# Create a virtual environment and install Python dependencies
RUN python3 -m venv /venv
ENV PATH="/venv/bin:$PATH"
RUN pip install --no-cache-dir \
        pandas matplotlib torch torchvision faiss-cpu \
        git+https://github.com/nicoaira/GINFINITY.git#egg=ginfinity

# Install Nextflow (so that processes launched inside the container can also call Nextflow if they need to)
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow

# Set the working directory and copy in the entire pipeline (so that
# when Nextflow mounts /app as /workspace inside the container, all scripts are present)
WORKDIR /app
COPY main.nf nextflow.config ./
COPY workflows/ workflows/
COPY bin/ bin/
COPY conf/ conf/
COPY presets/ presets/

# Make helper scripts in bin/ executable
RUN chmod +x bin/*.py