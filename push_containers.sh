#!/bin/bash
# Push all GINflow Docker containers to Docker Hub

set -e

echo "Pushing GINflow Docker containers to Docker Hub..."

# Login to Docker Hub
docker login

# Core modules
docker push nicoaira/ginflow-align-candidates:latest
docker push nicoaira/ginflow-build-faiss-index:latest
docker push nicoaira/ginflow-calculate-statistics:latest
docker push nicoaira/ginflow-cluster-seeds:latest
docker push nicoaira/ginflow-generate-window-vectors:latest
docker push nicoaira/ginflow-merge-embedding-chunks:latest
docker push nicoaira/ginflow-merge-window-chunks:latest
docker push nicoaira/ginflow-query-faiss-index:latest

# Drawing and reporting
docker push nicoaira/ginflow-draw-structures-rnartistcore:latest
docker push nicoaira/ginflow-draw-structures-r4rna:latest
docker push nicoaira/ginflow-generate-report:latest

# Extract meta map (uses pandas)
docker push nicoaira/ginflow-extract-meta-map:latest

# GINFINITY (node embeddings)
docker push nicoaira/ginfinity:latest

echo "All containers pushed successfully!"
