#!/bin/bash

# Quick validation test for UNIFY_INPUT module
# Tests MTX format conversion which requires no additional setup

set -e

echo "========================================================================"
echo "scAnnex UNIFY_INPUT Quick Validation Test"
echo "========================================================================"
echo ""

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "✗ Docker not found. Please install Docker or use manual testing."
    echo "  See test_data/README.md for manual testing instructions."
    exit 1
fi

echo "✓ Docker available"
echo ""

# Setup paths
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$TEST_DIR")"
MTX_INPUT="$TEST_DIR/mtx/filtered_feature_bc_matrix"
OUTPUT_DIR="$TEST_DIR/outputs"
OUTPUT_FILE="$OUTPUT_DIR/PBMC_MTX_quick_test.h5ad"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check MTX input exists
if [ ! -d "$MTX_INPUT" ]; then
    echo "✗ MTX test data not found: $MTX_INPUT"
    echo "  Please run the download script first"
    exit 1
fi

echo "Test Configuration:"
echo "  Input: $MTX_INPUT"
echo "  Output: $OUTPUT_FILE"
echo "  Container: python:3.10-slim"
echo ""

echo "----------------------------------------------------------------------"
echo "Step 1: Running UNIFY_INPUT conversion (MTX → H5AD)"
echo "----------------------------------------------------------------------"
echo ""

# Run conversion in Docker container
docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace \
    python:3.10-slim \
    bash -c "
        set -e
        echo 'Installing Python dependencies...'
        pip install -q scanpy numpy pandas anndata
        echo ''
        echo 'Running unify_input.py...'
        python bin/unify_input.py \
            --input test_data/mtx/filtered_feature_bc_matrix \
            --input-type mtx \
            --output test_data/outputs/PBMC_MTX_quick_test.h5ad \
            --sample-id PBMC_MTX_TEST \
            --batch batch1 \
            --condition control
    "

echo ""
echo "----------------------------------------------------------------------"
echo "Step 2: Validating output structure"
echo "----------------------------------------------------------------------"
echo ""

# Validate output
docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace \
    python:3.10-slim \
    bash -c "
        set -e
        pip install -q scanpy pandas
        python test_data/validate_output.py test_data/outputs/PBMC_MTX_quick_test.h5ad
    "

echo ""
echo "========================================================================"
echo "✓ QUICK TEST PASSED"
echo "========================================================================"
echo ""
echo "Output file: $OUTPUT_FILE"
if [ -f "$OUTPUT_FILE" ]; then
    echo "File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
fi
echo ""
echo "Next steps:"
echo "  1. Review output in: $OUTPUT_DIR"
echo "  2. Test H5AD and RDS formats (see README.md)"
echo "  3. Run full pipeline test with Nextflow"
echo "  4. Proceed to Dashboard implementation"
