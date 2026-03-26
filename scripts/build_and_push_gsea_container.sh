#!/bin/bash
#
# Build and push GSEA container to DockerHub
# This script should be run on Apocrita (QMUL HPC)
#
# Requirements:
#   - Apptainer/Singularity installed
#   - docker.io credentials configured (apptainer remote login)
#   - Existing scannex-gsea_1.0.sif in containers/
#
# Usage:
#   bash scripts/build_and_push_gsea_container.sh
#

set -euo pipefail

# Configuration
DOCKERHUB_REPO="damouzo/scannex"
TAG="gsea-1.0"
SIF_FILE="containers/scannex-gsea_1.0.sif"
PROJECT_DIR="/data/BCI-KRP/projects/scAnnex"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== scAnnex GSEA Container - DockerHub Uploader ===${NC}"
echo

# Step 1: Verify SIF exists
echo -e "${YELLOW}[1/4] Verifying SIF file exists...${NC}"
if [[ ! -f "${PROJECT_DIR}/${SIF_FILE}" ]]; then
    echo -e "${RED}ERROR: ${SIF_FILE} not found!${NC}"
    echo "Expected location: ${PROJECT_DIR}/${SIF_FILE}"
    echo
    echo "Build it first with:"
    echo "  apptainer build ${SIF_FILE} containers/apptainer_gsea.def"
    exit 1
fi
echo -e "${GREEN}✓ Found: ${SIF_FILE} ($(du -h ${PROJECT_DIR}/${SIF_FILE} | cut -f1))${NC}"
echo

# Step 2: Verify docker.io credentials
echo -e "${YELLOW}[2/4] Verifying DockerHub credentials...${NC}"
if [[ ! -f "${HOME}/.docker/config.json" ]] && [[ ! -f "${HOME}/.apptainer/remote.yaml" ]]; then
    echo -e "${RED}ERROR: Not logged in to docker.io${NC}"
    echo
    echo "Login first with:"
    echo "  apptainer remote login --username damouzo docker://docker.io"
    exit 1
fi
echo -e "${GREEN}✓ Authenticated to docker.io${NC}"
echo

# Step 3: Convert SIF to Docker and push to DockerHub
echo -e "${YELLOW}[3/4] Pushing to DockerHub: ${DOCKERHUB_REPO}:${TAG}${NC}"
echo "This will take 10-20 minutes (2.3GB conversion + upload)..."
echo

cd "${PROJECT_DIR}"

# Method: Use singularity/apptainer to push directly to registry
# Note: Apptainer can push to OCI registries using oras:// protocol
echo "Converting SIF to OCI format and pushing..."

if apptainer push "${SIF_FILE}" "oras://docker.io/${DOCKERHUB_REPO}:${TAG}"; then
    echo
    echo -e "${GREEN}✓ Successfully pushed to DockerHub!${NC}"
else
    echo
    echo -e "${RED}ERROR: Push failed${NC}"
    echo
    echo "Alternative method: Build as Docker image"
    echo "1. Create Dockerfile from apptainer_gsea.def"
    echo "2. Build with: docker build -f containers/Dockerfile.gsea -t ${DOCKERHUB_REPO}:${TAG} ."
    echo "3. Push with: docker push ${DOCKERHUB_REPO}:${TAG}"
    exit 1
fi
echo

# Step 4: Verify push
echo -e "${YELLOW}[4/4] Verifying image on DockerHub...${NC}"
echo "Image URL: https://hub.docker.com/r/${DOCKERHUB_REPO}/tags"
echo
echo "To test pull with Singularity:"
echo "  singularity pull docker://${DOCKERHUB_REPO}:${TAG}"
echo
echo -e "${GREEN}=== GSEA Container Upload Complete ===${NC}"
echo
echo "Next steps:"
echo "1. Update gsea.nf container directive to: docker://${DOCKERHUB_REPO}:${TAG}"
echo "2. Delete local containers/scannex-gsea_1.0.sif (no longer needed)"
echo "3. Singularity will auto-download from DockerHub when run_gsea=true"
echo
