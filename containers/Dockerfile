# Base image with Python and Java for Nextflow
FROM openjdk:11-jre-slim AS base

# Install OS packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip git curl && \
    rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip3 install --no-cache-dir \
    pandas matplotlib torch torchvision faiss-cpu

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/nextflow

# Copy pipeline scripts and set workdir
WORKDIR /app
COPY main.nf nextflow.config ./
COPY workflows/ workflows/
COPY bin/ bin/
COPY conf/ conf/
COPY presets/ presets/

# Make helper scripts executable
RUN chmod +x bin/*.py

# Default entrypoint
ENTRYPOINT ["nextflow"]
CMD ["run", "main.nf"]