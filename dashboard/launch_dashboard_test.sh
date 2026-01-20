#!/bin/bash
# Launch scAnnex Dashboard with detailed error reporting

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║              scAnnex Dashboard - Launch Script                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: conda not found in PATH"
    echo "   Please install Miniconda or Anaconda first"
    exit 1
fi

# Initialize conda
echo "🔧 Initializing conda..."
eval "$(conda shell.bash hook 2>/dev/null)" || {
    echo "❌ ERROR: Failed to initialize conda"
    exit 1
}

# Activate environment
ENV_NAME="scannex-dashboard"
echo "🔧 Activating environment: $ENV_NAME"

if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "❌ ERROR: Conda environment '${ENV_NAME}' not found"
    echo ""
    echo "Please create it first:"
    echo "  ./setup_dashboard.sh"
    exit 1
fi

conda activate "$ENV_NAME" || {
    echo "❌ ERROR: Failed to activate conda environment"
    exit 1
}

# Set Python path
export RETICULATE_PYTHON="$(which python3)"
echo "✓ Python path: $RETICULATE_PYTHON"

# Check required files
echo ""
echo "📋 Checking required files..."
REQUIRED_FILES=("app.R" "global.R" "server.R" "ui.R")
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "  ❌ Missing: $file"
        exit 1
    else
        echo "  ✓ Found: $file"
    fi
done

# Check test data
TEST_DATA="/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad"
if [ -f "$TEST_DATA" ]; then
    echo "  ✓ Test data: $TEST_DATA"
else
    echo "  ⚠️  Test data not found: $TEST_DATA"
    echo "     (You'll need to specify a different file in the dashboard)"
fi

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "  Launching Dashboard..."
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "  URL: http://localhost:8888"
echo "  Press Ctrl+C to stop"
echo ""
echo "═══════════════════════════════════════════════════════════════"
echo ""

# Launch with R
R --quiet --no-save -e "shiny::runApp('.', host='127.0.0.1', port=8888, launch.browser=FALSE)" 2>&1 | while IFS= read -r line; do
    # Filter out the numpy warning
    if [[ ! "$line" =~ "UserWarning: Signature" ]] && \
       [[ ! "$line" =~ "falling back to type probe" ]] && \
       [[ ! "$line" =~ "machar = _get_machar" ]] && \
       [[ ! "$line" =~ "getlimits.py" ]]; then
        echo "$line"
    fi
done
