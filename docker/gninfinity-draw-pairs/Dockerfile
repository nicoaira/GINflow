# Dockerfile for gin-draw-pairs
FROM python:3.9-slim

# Install system dependencies:
# - Git, Maven, Java for RNArtistCore
# - Cairo, Pango, GDK-Pixbuf for cairosvg rendering support
RUN apt-get update && DEBIAN_FRONTEND=noninteractive \
    apt-get install -y --no-install-recommends \
        git \
        maven \
        default-jdk-headless \
        libcairo2 \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libgdk-pixbuf2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Clone and build RNArtistCore
RUN git clone https://github.com/fjossinet/RNArtistCore.git /opt/RNArtistCore \
    && cd /opt/RNArtistCore \
    && mvn clean package

# Create a wrapper script for `rnartistcore`
RUN printf '#!/usr/bin/env bash\nJAR=$(ls /opt/RNArtistCore/target/rnartistcore-*-SNAPSHOT-jar-with-dependencies.jar)\nexec java -jar "$JAR" "$@"\n' > /usr/local/bin/rnartistcore \
    && chmod +x /usr/local/bin/rnartistcore

# Set working directory
WORKDIR /app

# Copy the draw_pairs script into the container
COPY nextflow/modules/draw_pairs.py ./draw_pairs.py

# Install required Python packages
RUN pip install cairosvg pandas