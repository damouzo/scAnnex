#!/usr/bin/env bash

#==============================================================================
# test_integration_quick.sh
#==============================================================================
# Quick integration test with optimized parameters for small test dataset
# (1,049 cells × 14,859 genes)
#
# Optimizations:
#   - Reduced HVGs: 2000 → 1000 (sufficient for small dataset)
#   - Reduced PCs: 50 → 20 (appropriate for 1k cells)
#   - Reduced neighbors: 15 → 10 (better for small datasets)
#   - Reduced Harmony iterations: 10 → 5 (faster convergence)
#==============================================================================

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/analytical_core_results"
INTEGRATION_DIR="$OUTPUT_DIR/integration_results"

# Docker image
PYTHON_IMAGE="python:3.11-slim"

# Input/Output files
INPUT_FILE="$OUTPUT_DIR/qc_filtered.h5ad"
OUTPUT_FILE="$OUTPUT_DIR/normalized_integrated.h5ad"
LOG_FILE="$OUTPUT_DIR/integration_quick_test.log"

#==============================================================================
# Main
#==============================================================================

log_info "Starting quick integration test with optimized parameters..."
echo "Project root: $PROJECT_ROOT"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Check input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    log_error "Input file not found: $INPUT_FILE"
    log_error "Run the QC module first: cd test_data && bash test_analytical_core.sh"
    exit 1
fi

# Create output directory
mkdir -p "$INTEGRATION_DIR"

# Display input file info
log_info "Input file size: $(du -h "$INPUT_FILE" | cut -f1)"

# Run integration with optimized parameters
log_info "Running integration with optimized parameters..."
echo "  - HVGs: 1000 (reduced from 2000)"
echo "  - PCs: 20 (reduced from 50)"
echo "  - Neighbors: 10 (reduced from 15)"
echo "  - Harmony iterations: 5 (reduced from 10)"
echo ""

START_TIME=$(date +%s)

docker run --rm \
    -v "$PROJECT_ROOT:/workspace" \
    -w /workspace \
    "$PYTHON_IMAGE" bash -c "
        echo 'Installing Python packages...'
        pip install -q scanpy numpy pandas matplotlib seaborn harmonypy scikit-learn 2>&1 | grep -v 'Requirement already satisfied' || true
        echo 'Running normalize_integrate.py...'
        python bin/normalize_integrate.py \
            --input test_data/analytical_core_results/qc_filtered.h5ad \
            --output test_data/analytical_core_results/normalized_integrated.h5ad \
            --output-dir test_data/analytical_core_results/integration_results \
            --run-integration \
            --batch-key batch \
            --n-top-genes 1000 \
            --n-pcs 20 \
            --n-neighbors 10 \
            --harmony-theta 2.0 \
            --harmony-max-iter 5
    " 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""

if [[ $EXIT_CODE -eq 0 ]]; then
    log_success "Integration completed successfully in ${DURATION}s"
    
    # Verify output file
    if [[ -f "$OUTPUT_FILE" ]]; then
        OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        log_success "Output file created: $OUTPUT_FILE ($OUTPUT_SIZE)"
        
        # Verify output structure using Python
        log_info "Verifying output structure..."
        docker run --rm \
            -v "$PROJECT_ROOT:/workspace" \
            -w /workspace \
            "$PYTHON_IMAGE" bash -c "
                pip install -q scanpy 2>&1 | grep -v 'Requirement already satisfied' || true
                python -c \"
import scanpy as sc
import sys

try:
    adata = sc.read_h5ad('test_data/analytical_core_results/normalized_integrated.h5ad')
    
    print(f'Shape: {adata.shape[0]} cells × {adata.shape[1]} genes')
    print(f'Observations (metadata): {list(adata.obs.columns)}')
    print(f'Embeddings: {list(adata.obsm.keys())}')
    
    # Critical checks
    assert 'X_pca' in adata.obsm, 'Missing X_pca!'
    assert 'X_pca_harmony' in adata.obsm, 'Missing X_pca_harmony!'
    assert 'X_umap' in adata.obsm, 'Missing X_umap!'
    
    print(f'\\n✓ X_pca shape: {adata.obsm[\"X_pca\"].shape}')
    print(f'✓ X_pca_harmony shape: {adata.obsm[\"X_pca_harmony\"].shape}')
    print(f'✓ X_umap shape: {adata.obsm[\"X_umap\"].shape}')
    print('\\n✓ All required embeddings present!')
    sys.exit(0)
    
except Exception as e:
    print(f'ERROR: {e}', file=sys.stderr)
    sys.exit(1)
\"
            "
        
        if [[ $? -eq 0 ]]; then
            log_success "Output structure verified successfully"
            echo ""
            log_info "Integration results saved to: $INTEGRATION_DIR"
            log_info "Next steps:"
            echo "  1. Review integration plots in: $INTEGRATION_DIR"
            echo "  2. Launch dashboard: cd dashboard && ./run_dashboard.sh run"
            echo "  3. Point dashboard to: /srv/shiny-server/data/normalized_integrated.h5ad"
        else
            log_error "Output structure verification failed"
            exit 1
        fi
    else
        log_error "Output file not created: $OUTPUT_FILE"
        exit 1
    fi
else
    log_error "Integration failed with exit code: $EXIT_CODE"
    log_error "Check log file: $LOG_FILE"
    exit 1
fi

echo ""
log_success "Quick integration test completed successfully!"
