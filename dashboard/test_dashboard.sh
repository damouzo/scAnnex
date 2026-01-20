#!/bin/bash
#
# Quick test of dashboard functionality
# Checks that all components work before launching full server
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-${SCRIPT_DIR}/../results_slc_first_run}"

echo "Testing scAnnex Dashboard Components..."
echo ""

# Activate conda if available
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    if conda env list | grep -q "scannex-dashboard"; then
        echo "✓ Activating scannex-dashboard environment"
        conda activate scannex-dashboard
    fi
fi

# Test 1: Check R packages
echo "→ Testing R packages..."
R --quiet --no-save << 'EOF'
packages <- c("shiny", "shinydashboard", "shinyWidgets", "reticulate", 
              "plotly", "DT", "ggplot2", "viridis", "data.table", "jsonlite")

missing <- c()
for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing <- c(missing, pkg)
    }
}

if (length(missing) > 0) {
    cat("✗ Missing R packages:", paste(missing, collapse=", "), "\n")
    quit(status = 1)
} else {
    cat("✓ All R packages found\n")
}
EOF

# Test 2: Check Python packages
echo "→ Testing Python packages..."
python3 << 'EOF'
try:
    import scanpy
    import anndata
    import numpy
    import pandas
    print("✓ All Python packages found")
except ImportError as e:
    print(f"✗ Missing Python package: {e}")
    exit(1)
EOF

# Test 3: Check data file
echo "→ Checking for data file..."
if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "✗ Results directory not found: $RESULTS_DIR"
    exit 1
fi

H5AD_FILE=$(find "$RESULTS_DIR" -name "*annotated*.h5ad" -type f | head -1)
if [[ -z "$H5AD_FILE" ]]; then
    echo "✗ No annotated .h5ad file found in: $RESULTS_DIR"
    exit 1
fi

echo "✓ Found data file: $(basename $H5AD_FILE)"

# Test 4: Quick data load test
echo "→ Testing data loading..."
python3 << EOF
import anndata as ad
import sys

try:
    adata = ad.read_h5ad("$H5AD_FILE")
    print(f"✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Check for required fields
    if 'X_umap' not in adata.obsm:
        print("⚠ Warning: No UMAP coordinates found")
    else:
        print("✓ UMAP coordinates present")
    
    if 'predicted_labels' in adata.obs.columns:
        print("✓ Cell type annotations present")
        n_types = adata.obs['predicted_labels'].nunique()
        print(f"  Found {n_types} cell types")
    else:
        print("⚠ Warning: No predicted_labels in metadata")
        
except Exception as e:
    print(f"✗ Error loading data: {e}")
    sys.exit(1)
EOF

# Test 5: R + Python integration (reticulate)
echo "→ Testing R-Python integration..."
R --quiet --no-save << EOF
# Set Python path BEFORE loading reticulate (same as global.R)
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix != "") {
    python_path <- file.path(conda_prefix, "bin", "python3")
} else {
    python_path <- "/usr/bin/python3"
}
Sys.setenv(RETICULATE_PYTHON = python_path)

# Now load reticulate
library(reticulate)
use_python(python_path, required = TRUE)

tryCatch({
    ad <- import("anndata")
    adata <- ad\$read_h5ad("$H5AD_FILE")
    cat("✓ R can access Python and read H5AD files\n")
    cat(sprintf("  %d cells loaded via reticulate\n", adata\$n_obs))
}, error = function(e) {
    cat("✗ R-Python integration failed:", e\$message, "\n")
    quit(status = 1)
})
EOF

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "✓ All tests passed! Dashboard is ready to launch"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "To start the dashboard:"
echo "  ./launch_dashboard.sh $RESULTS_DIR"
echo ""
