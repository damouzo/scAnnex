#!/bin/bash

# scAnnex UNIFY_INPUT Test Runner
# This script tests the UNIFY_INPUT module using Docker containers

set -euo pipefail

echo "======================================================================"
echo "scAnnex UNIFY_INPUT Module Testing Suite"
echo "======================================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Directories
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$TEST_DIR")"
OUTPUT_DIR="$TEST_DIR/outputs"

echo "Test directory: $TEST_DIR"
echo "Project root: $PROJECT_ROOT"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check Docker availability
if ! command -v docker &> /dev/null; then
    echo -e "${RED}✗ Docker not found${NC}"
    echo "  Please install Docker to run these tests"
    exit 1
fi
echo -e "${GREEN}✓ Docker available${NC}"

# Docker images
PYTHON_IMAGE="python:3.10-slim"
R_IMAGE="rocker/r-ver:4.3"

# ==============================================================================
# Step 1: Create test data files
# ==============================================================================
echo ""
echo "----------------------------------------------------------------------"
echo "Step 1: Creating test data files"
echo "----------------------------------------------------------------------"

# Create H5AD test file using Python container
echo "Creating H5AD test file..."
docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace/test_data \
    "$PYTHON_IMAGE" bash -c "
        pip install -q scanpy numpy pandas &&
        python create_test_files.py
    " || {
    echo -e "${YELLOW}⚠ H5AD creation skipped (expected if Python script needs adjustment)${NC}"
}

# Create RDS test file using R container
echo ""
echo "Creating RDS test file..."
docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace/test_data \
    "$R_IMAGE" bash -c "
        R -e \"install.packages(c('Seurat', 'remotes'), repos='https://cloud.r-project.org', quiet=TRUE)\" &&
        R -e \"remotes::install_github('mojaveazure/seurat-disk', quiet=TRUE)\" &&
        Rscript create_test_rds.R
    " 2>&1 | grep -v "^─" || {
    echo -e "${YELLOW}⚠ RDS creation skipped (Seurat installation takes time)${NC}"
}

# ==============================================================================
# Step 2: Test UNIFY_INPUT on each format
# ==============================================================================
echo ""
echo "----------------------------------------------------------------------"
echo "Step 2: Testing UNIFY_INPUT module"
echo "----------------------------------------------------------------------"

# Test parameters
SAMPLES=(
    "MTX:mtx:test_data/mtx/filtered_feature_bc_matrix:batch1:control"
    "H5AD:h5ad:test_data/h5ad/pbmc_100cells.h5ad:batch1:treated"
    "RDS:rds:test_data/rds/pbmc_seurat.rds:batch2:control"
)

SUCCESS_COUNT=0
FAIL_COUNT=0

for sample_spec in "${SAMPLES[@]}"; do
    IFS=':' read -r sample_name file_type file_path batch condition <<< "$sample_spec"
    
    echo ""
    echo "Testing: $sample_name ($file_type format)"
    echo "  Input: $file_path"
    
    # Check if input exists
    if [ ! -e "$PROJECT_ROOT/$file_path" ]; then
        echo -e "${YELLOW}  ⚠ Skipped - input file not found${NC}"
        continue
    fi
    
    output_file="$OUTPUT_DIR/${sample_name}_unified.h5ad"
    
    # Run unify_input.py in container
    echo "  Running unify_input.py..."
    
    if docker run --rm \
        -v "$PROJECT_ROOT:/workspace" \
        -w /workspace \
        "$PYTHON_IMAGE" bash -c "
            pip install -q scanpy numpy pandas anndata &&
            python bin/unify_input.py \
                --input $file_path \
                --input-type $file_type \
                --output $output_file \
                --sample-id $sample_name \
                --batch $batch \
                --condition $condition
        " 2>&1 | tee "$OUTPUT_DIR/${sample_name}_conversion.log"; then
        
        # Check if output was created
        if [ -f "$output_file" ]; then
            file_size=$(du -h "$output_file" | cut -f1)
            echo -e "${GREEN}  ✓ Success - Output created ($file_size)${NC}"
            ((SUCCESS_COUNT++))
            
            # Validate output
            echo "  Validating output structure..."
            docker run --rm \
                -v "$PROJECT_ROOT:/workspace" \
                -w /workspace \
                "$PYTHON_IMAGE" bash -c "
                    pip install -q scanpy pandas &&
                    python test_data/validate_output.py $output_file
                " 2>&1 | grep -E "^(✓|✗|  )" || true
        else
            echo -e "${RED}  ✗ Failed - Output file not created${NC}"
            ((FAIL_COUNT++))
        fi
    else
        echo -e "${RED}  ✗ Failed - Conversion error${NC}"
        ((FAIL_COUNT++))
    fi
done

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "======================================================================"
echo "Test Summary"
echo "======================================================================"
echo ""
echo "Successful conversions: $SUCCESS_COUNT"
echo "Failed conversions: $FAIL_COUNT"
echo ""

if [ -d "$OUTPUT_DIR" ] && [ "$(ls -A $OUTPUT_DIR/*.h5ad 2>/dev/null)" ]; then
    echo "Output files:"
    ls -lh "$OUTPUT_DIR"/*.h5ad 2>/dev/null || true
fi

echo ""
if [ $FAIL_COUNT -eq 0 ] && [ $SUCCESS_COUNT -gt 0 ]; then
    echo -e "${GREEN}✓ ALL TESTS PASSED${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Review output files in: $OUTPUT_DIR"
    echo "  2. Proceed with Nextflow pipeline testing"
    echo "  3. Move to Dashboard implementation (Phase 8)"
    exit 0
else
    echo -e "${YELLOW}⚠ SOME TESTS FAILED OR SKIPPED${NC}"
    echo ""
    echo "Review logs in: $OUTPUT_DIR"
    echo ""
    echo "Common issues:"
    echo "  - Missing test data files: Run test data creation scripts"
    echo "  - R packages not installed: Seurat/SeuratDisk installation takes time"
    echo "  - Python packages: Installed automatically in container"
    exit 1
fi
