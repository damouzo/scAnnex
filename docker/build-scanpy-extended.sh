#!/bin/bash
# Build extended Scanpy container for scAnnex pipeline
# This container includes scanpy + scrublet + harmonypy + celltypist

set -e

IMAGE_NAME="scannex/scanpy-extended"
IMAGE_TAG="1.7.2"

echo "Building scAnnex extended Scanpy container..."
docker build -t ${IMAGE_NAME}:${IMAGE_TAG} -f docker/Dockerfile.scanpy-extended .

echo ""
echo "Build complete! Tagged as: ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""
echo "To use this container, update your module files with:"
echo "  container \"${IMAGE_NAME}:${IMAGE_TAG}\""
echo ""
echo "To push to a registry (optional):"
echo "  docker tag ${IMAGE_NAME}:${IMAGE_TAG} your-registry/${IMAGE_NAME}:${IMAGE_TAG}"
echo "  docker push your-registry/${IMAGE_NAME}:${IMAGE_TAG}"
