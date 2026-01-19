#!/bin/bash

# scAnnex Analytical Core Test Runner
# Tests QC and Integration modules with PBMC test data

set -euo pipefail

echo "======================================================================"
echo "scAnnex Analytical Core Testing (Phases 2 & 4)"
echo "======================================================================"
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Directories
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$TEST_DIR")"
INPUT_FILE="$TEST_DIR/outputs/PBMC_MTX_quick_test.h5ad"
OUTPUT_DIR="$TEST_DIR/analytical_core_results"

echo "Test directory: $TEST_DIR"
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check Docker
if ! command -v docker &> /dev/null; then
    echo -e "${YELLOW}✗ Docker not found - cannot run tests${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Docker available${NC}"

# Check input file
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${YELLOW}✗ Input file not found: $INPUT_FILE${NC}"
    echo "  Please run test_data/quick_test.sh first to generate test data"
    exit 1
fi
echo -e "${GREEN}✓ Input file exists${NC}"
echo ""

# Docker image
PYTHON_IMAGE="python:3.10-slim"

# ==============================================================================
# PHASE 2: Quality Control Testing
# ==============================================================================
echo "======================================================================"
echo "Phase 2: Testing Enhanced QC Module"
echo "======================================================================"
echo ""

QC_OUTPUT="$OUTPUT_DIR/qc_filtered.h5ad"
QC_DIR="$OUTPUT_DIR/qc_results"

echo "Running quality_control.py with MAD thresholds..."
docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace \
    "$PYTHON_IMAGE" bash -c "
        pip install -q scanpy numpy pandas matplotlib seaborn scikit-learn &&
        python bin/quality_control.py \
            --input test_data/outputs/PBMC_MTX_quick_test.h5ad \
            --output test_data/analytical_core_results/qc_filtered.h5ad \
            --qc-dir test_data/analytical_core_results/qc_results \
            --use-mad-thresholds \
            --mad-threshold 5.0
    " 2>&1 | tee "$OUTPUT_DIR/qc_test.log"

if [ -f "$QC_OUTPUT" ]; then
    echo -e "\n${GREEN}✓ QC module test PASSED${NC}"
    echo "  Output: $QC_OUTPUT"
    echo "  QC plots: $QC_DIR/"
else
    echo -e "\n${YELLOW}✗ QC module test FAILED${NC}"
    exit 1
fi

# ==============================================================================
# PHASE 4: Normalization & Integration Testing  
# ==============================================================================
echo ""
echo "======================================================================"
echo "Phase 4: Testing Normalization & Integration Module"
echo "======================================================================"
echo ""

NORM_OUTPUT="$OUTPUT_DIR/normalized_integrated.h5ad"
INT_DIR="$OUTPUT_DIR/integration_results"

echo "Running normalize_integrate.py..."
docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace \
    "$PYTHON_IMAGE" bash -c "
        pip install -q scanpy numpy pandas matplotlib seaborn harmonypy scikit-learn &&
        python bin/normalize_integrate.py \
            --input test_data/analytical_core_results/qc_filtered.h5ad \
            --output test_data/analytical_core_results/normalized_integrated.h5ad \
            --output-dir test_data/analytical_core_results/integration_results \
            --run-integration \
            --batch-key batch \
            --n-top-genes 2000 \
            --n-pcs 30 \
            --harmony-theta 2.0
    " 2>&1 | tee "$OUTPUT_DIR/integration_test.log"

if [ -f "$NORM_OUTPUT" ]; then
    echo -e "\n${GREEN}✓ Integration module test PASSED${NC}"
    echo "  Output: $NORM_OUTPUT"
    echo "  Integration plots: $INT_DIR/"
else
    echo -e "\n${YELLOW}✗ Integration module test FAILED${NC}"
    exit 1
fi

# ==============================================================================
# Validation
# ==============================================================================
echo ""
echo "======================================================================"
echo "Validation: Checking Output Structure"
echo "======================================================================"
echo ""

docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace \
    "$PYTHON_IMAGE" bash -c "
        pip install -q scanpy pandas &&
        python -c \"
import scanpy as sc
import sys

print('Loading final output...')
adata = sc.read_h5ad('test_data/analytical_core_results/normalized_integrated.h5ad')

print(f'\\n✓ Dimensions: {adata.n_obs} cells × {adata.n_vars} genes')
print(f'\\n.obsm keys (should include PCA and Harmony):')
for key in adata.obsm.keys():
    print(f'  - {key}: {adata.obsm[key].shape}')

print(f'\\n.layers keys (should include counts):')
for key in adata.layers.keys():
    print(f'  - {key}')

print(f'\\n.obs columns (should include QC metrics):')
qc_cols = [col for col in adata.obs.columns if 'pct_' in col or 'n_genes' in col or 'total_counts' in col]
for col in qc_cols[:5]:
    print(f'  - {col}')

# Validate required fields
required_obsm = ['X_pca', 'X_pca_harmony', 'X_umap']
missing = [k for k in required_obsm if k not in adata.obsm.keys()]

if missing:
    print(f'\\n✗ Missing required .obsm keys: {missing}')
    sys.exit(1)
else:
    print(f'\\n✓ All required .obsm keys present')

print(f'\\n✓ VALIDATION PASSED')
        \"
    "

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "======================================================================"
echo "Test Summary"
echo "======================================================================"
echo ""
echo "Test Results:"
echo "  ✓ Phase 2 (QC): PASSED"
echo "  ✓ Phase 4 (Integration): PASSED"
echo "  ✓ Output validation: PASSED"
echo ""
echo "Generated Files:"
echo "  1. QC filtered data: $QC_OUTPUT"
echo "  2. QC plots: $QC_DIR/"
echo "  3. Integrated data: $NORM_OUTPUT"
echo "  4. Integration plots: $INT_DIR/"
echo ""
echo "Test Logs:"
echo "  - QC: $OUTPUT_DIR/qc_test.log"
echo "  - Integration: $OUTPUT_DIR/integration_test.log"
echo ""
echo -e "${GREEN}✓ ALL ANALYTICAL CORE TESTS PASSED${NC}"
echo ""
echo "Next Steps:"
echo "  1. Review QC and integration plots"
echo "  2. Examine integration metrics in $INT_DIR/integration_report.json"
echo "  3. Proceed with clustering and annotation (Phase 6)"
echo "  4. Implement Dashboard (Phase 8)"
