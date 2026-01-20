#!/bin/bash
################################################################################
# scAnnex SLC Pipeline Launch Script
# 
# This script sets up the environment and launches the full SLC pipeline
# from UNIFY_INPUT to AUTO_ANNOT_CELLTYPIST
#
# Hardware: 16GB RAM, 4 CPUs (WSL2)
# Environment: Conda (scannex-minimal)
################################################################################

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}  scAnnex SLC Pipeline - Production Launch${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo

# Set paths
export PATH="/home/damo/miniforge3/bin:$PATH"
PROJECT_DIR="/home/damo/scAnnex"
cd "${PROJECT_DIR}"

# Activate conda environment
echo -e "${GREEN}[1/4] Activating conda environment (scannex)...${NC}"
source /home/damo/miniforge3/bin/activate scannex

# Verify environment
echo -e "${GREEN}[2/4] Verifying environment...${NC}"
python -c "import scanpy as sc; print(f'  ✓ Scanpy {sc.__version__}')"
python -c "import anndata as ad; print(f'  ✓ AnnData {ad.__version__}')"
python -c "import celltypist as ct; print(f'  ✓ CellTypist {ct.__version__}')"
nextflow -version | head -3
echo

# Default parameters (can be overridden via command line)
INPUT_SAMPLESHEET="${1:-test_data/samplesheet_slc_test.csv}"
OUTPUT_DIR="${2:-results_slc_$(date +%Y%m%d_%H%M%S)}"
RESUME="${3:-}"

echo -e "${GREEN}[3/4] Pipeline Configuration:${NC}"
echo "  • Input: ${INPUT_SAMPLESHEET}"
echo "  • Output: ${OUTPUT_DIR}"
echo "  • Profile: conda,laptop (8GB RAM, 4 CPUs)"
echo "  • Doublet Detection: ENABLED"
echo "  • CellTypist Annotation: ENABLED"
echo "  • Batch Integration: DISABLED (single sample)"
echo

# Display sample information
echo -e "${YELLOW}Sample Information:${NC}"
cat "${INPUT_SAMPLESHEET}"
echo

# Confirm launch
echo -e "${YELLOW}Ready to launch pipeline. Press ENTER to continue or Ctrl+C to cancel...${NC}"
read -r

# Launch pipeline
echo -e "${GREEN}[4/4] Launching scAnnex SLC Pipeline...${NC}"
echo

RESUME_FLAG=""
if [ -n "${RESUME}" ]; then
    RESUME_FLAG="-resume"
    echo "  → Resume mode: ENABLED"
fi

nextflow run main.nf \
    -profile conda,laptop \
    --input "${INPUT_SAMPLESHEET}" \
    --outdir "${OUTPUT_DIR}" \
    --run_doublet_detection true \
    --doublet_removal true \
    --run_auto_annotation true \
    --annotation_method celltypist \
    --celltypist_model 'Immune_All_Low.pkl' \
    --run_integration false \
    -with-report "${OUTPUT_DIR}/pipeline_report.html" \
    -with-timeline "${OUTPUT_DIR}/timeline.html" \
    -with-dag "${OUTPUT_DIR}/dag.html" \
    ${RESUME_FLAG}

EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo
    echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}  ✓ Pipeline completed successfully!${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}"
    echo
    echo -e "${BLUE}Output Location:${NC}"
    echo "  ${OUTPUT_DIR}/"
    echo
    echo -e "${BLUE}Key Outputs:${NC}"
    echo "  • Unified H5AD: ${OUTPUT_DIR}/unified/"
    echo "  • QC Results: ${OUTPUT_DIR}/qc/"
    echo "  • Doublet Scores: ${OUTPUT_DIR}/doublets/"
    echo "  • Processed Data: ${OUTPUT_DIR}/standard_processing_results/"
    echo "  • Cell Type Annotations: ${OUTPUT_DIR}/celltypist/"
    echo "  • Pipeline Report: ${OUTPUT_DIR}/pipeline_report.html"
    echo
    echo -e "${BLUE}Next Steps:${NC}"
    echo "  1. Review QC metrics and doublet detection results"
    echo "  2. Inspect cell type annotations"
    echo "  3. Launch dashboard for interactive visualization:"
    echo "     cd dashboard && docker run -p 3838:3838 -v \$(pwd)/../${OUTPUT_DIR}:/data scannex-dashboard"
else
    echo
    echo -e "${YELLOW}═══════════════════════════════════════════════════════${NC}"
    echo -e "${YELLOW}  Pipeline finished with errors (exit code: ${EXIT_CODE})${NC}"
    echo -e "${YELLOW}═══════════════════════════════════════════════════════${NC}"
    echo
    echo "Check the error logs in:"
    echo "  ${OUTPUT_DIR}/pipeline_info/"
    echo
    echo "To resume the pipeline from where it stopped:"
    echo "  ./run_slc_pipeline.sh ${INPUT_SAMPLESHEET} ${OUTPUT_DIR} resume"
fi

exit ${EXIT_CODE}
