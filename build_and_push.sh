#!/bin/bash

# Build and push script for GINflow Docker images
# Usage: ./build_and_push.sh [username] [tag]

USERNAME=${1:-nicoaira}
TAG=${2:-latest}

# Change to the project root directory
cd /home/nicolas/programs/GINflow

echo "Building and pushing Docker images for GINflow..."
echo "Username: $USERNAME"
echo "Tag: $TAG"
echo ""

# List of all Dockerfiles and their corresponding image names
declare -A dockerfiles=(
    ["Dockerfile.aggregate_score"]="ginflow-aggregate-score"
    ["Dockerfile.build_faiss_index"]="ginflow-build-faiss-index"
    ["Dockerfile.draw_contig_svgs"]="ginflow-draw-contig-svgs"
    ["Dockerfile.draw_unagg_svgs"]="ginflow-draw-unagg-svgs"
    ["Dockerfile.extract_meta_map"]="ginflow-extract-meta-map"
    ["Dockerfile.filter_top_contigs"]="ginflow-filter-top-contigs"
    ["Dockerfile.generate_aggregated_report"]="ginflow-generate-aggregated-report"
    ["Dockerfile.generate_embeddings"]="ginflow-generate-embeddings"
    ["Dockerfile.generate_unaggregated_report"]="ginflow-generate-unaggregated-report"
    ["Dockerfile.generate_windows"]="ginflow-generate-windows"
    ["Dockerfile.plot_distances"]="ginflow-plot-distances"
    ["Dockerfile.plot_score"]="ginflow-plot-score"
    ["Dockerfile.query_faiss_index"]="ginflow-query-faiss-index"
    ["Dockerfile.sort_distances"]="ginflow-sort-distances"
)

# Build and push each image
for dockerfile in "${!dockerfiles[@]}"; do
    image_name="${dockerfiles[$dockerfile]}"
    full_tag="$USERNAME/$image_name:$TAG"
    
    echo "Building $image_name..."
    if docker build -f containers/$dockerfile -t $full_tag .; then
        echo "✓ Built $full_tag successfully"
        
        echo "Pushing $full_tag..."
        if docker push $full_tag; then
            echo "✓ Pushed $full_tag successfully"
        else
            echo "✗ Failed to push $full_tag"
        fi
    else
        echo "✗ Failed to build $full_tag"
    fi
    echo ""
done

echo "Build and push process completed!"