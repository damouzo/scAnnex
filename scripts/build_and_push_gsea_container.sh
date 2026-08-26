#!/bin/bash
#
# Build and push GSEA container to GHCR (ghcr.io/damouzo/scannex/gsea)
# Normally built automatically on push by .github/workflows/build-containers.yml
# from containers/gsea/Dockerfile. This script is the manual fallback (Apocrita).
#
# Requirements:
#   - Apptainer/Singularity installed (build) and docker CLI (push)
#   - ghcr.io credentials configured (`docker login ghcr.io`)
#
# Usage:
#   bash scripts/build_and_push_gsea_container.sh
#

set -euo pipefail

# Configuration
GHCR_REPO="ghcr.io/damouzo/scannex/gsea"
TAG="latest"
SIF_FILE="containers/scannex-gsea_latest.sif"
PROJECT_DIR="/data/BCI-KRP/projects/scAnnex"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== scAnnex GSEA Container - GHCR Uploader ===${NC}"
echo

# Step 1: Verify SIF exists
echo -e "${YELLOW}[1/4] Verifying / building SIF file...${NC}"
if [[ ! -f "${PROJECT_DIR}/${SIF_FILE}" ]]; then
    echo "SIF not found. Building from containers/gsea/Dockerfile..."
    echo "  (recommended: build the Dockerfile in CI; this builds from apptainer def)"
    apptainer build "${PROJECT_DIR}/${SIF_FILE}" "${PROJECT_DIR}/containers/apptainer_gsea.def"
fi
echo -e "${GREEN}✓ Found: ${SIF_FILE} ($(du -h ${PROJECT_DIR}/${SIF_FILE} | cut -f1))${NC}"
echo

# Step 2: Verify ghcr.io credentials
echo -e "${YELLOW}[2/4] Verifying GHCR credentials...${NC}"
if [[ ! -f "${HOME}/.docker/config.json" ]]; then
    echo -e "${RED}ERROR: Not logged in to ghcr.io${NC}"
    echo
    echo "Login first with:"
    echo "  docker login ghcr.io"
    exit 1
fi
echo -e "${GREEN}✓ Authenticated to ghcr.io${NC}"
echo

# Step 3: Convert SIF to Docker and push to GHCR
echo -e "${YELLOW}[3/4] Pushing to GHCR: ${GHCR_REPO}:${TAG}${NC}"
echo "This will take 10-20 minutes (2.3GB conversion + upload)..."
echo

cd "${PROJECT_DIR}"

echo "Converting SIF to OCI format and pushing (oras:// protocol)..."

if apptainer push "${SIF_FILE}" "oras://${GHCR_REPO}:${TAG}"; then
    echo
    echo -e "${GREEN}✓ Successfully pushed to GHCR!${NC}"
else
    echo
    echo -e "${RED}ERROR: Push failed${NC}"
    echo
    echo "Alternative method: Build as Docker image"
    echo "1. Build with: docker build -f containers/gsea/Dockerfile -t ${GHCR_REPO}:${TAG} containers/gsea"
    echo "2. Push with:  docker push ${GHCR_REPO}:${TAG}"
    exit 1
fi
echo

# Step 4: Verify push
echo -e "${YELLOW}[4/4] Verifying image on GHCR...${NC}"
echo "Image URL: https://github.com/damouzo/scannex/pkgs/container/scannex%2Fgsea"
echo
echo "To test pull with Singularity:"
echo "  singularity pull docker://${GHCR_REPO}:${TAG}"
echo
echo -e "${GREEN}=== GSEA Container Upload Complete ===${NC}"
echo
echo "Next steps:"
echo "1. gsea.nf already references oras://${GHCR_REPO}:${TAG}"
echo "2. Singularity will auto-download from GHCR when run_gsea=true"
echo
