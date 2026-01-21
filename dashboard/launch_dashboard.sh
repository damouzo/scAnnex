#!/bin/bash
#
# scAnnex Dashboard Launcher
# Auto-setup and launch the interactive dashboard in one command
#
# Usage: ./launch_dashboard.sh [path/to/results_directory]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-${SCRIPT_DIR}/../results}"

# Convert relative path to absolute path
if [[ ! "${RESULTS_DIR}" = /* ]]; then
    RESULTS_DIR="$(cd "${SCRIPT_DIR}" && cd "${RESULTS_DIR}" && pwd)"
fi

# Colors
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

print_header() {
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Interactive Dashboard${NC}                           ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
}

print_success() { echo -e "${GREEN}✓${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }
print_info() { echo -e "${YELLOW}ℹ${NC} $1"; }

# Check if conda is available
if ! command -v conda &> /dev/null; then
    print_header
    echo ""
    print_error "Conda not found!"
    echo ""
    echo "Please install Conda/Mamba first:"
    echo "  https://github.com/conda-forge/miniforge"
    echo ""
    exit 1
fi

# Validate results directory
if [[ ! -d "$RESULTS_DIR" ]]; then
    print_header
    echo ""
    print_error "Results directory not found: $RESULTS_DIR"
    echo ""
    echo "Usage: $0 [path/to/results_directory]"
    exit 1
fi

# Find h5ad file
H5AD_FILE=$(find "$RESULTS_DIR" -name "*annotated*.h5ad" -type f | head -1)
if [[ -z "$H5AD_FILE" ]]; then
    print_header
    echo ""
    print_error "No annotated .h5ad file found in: $RESULTS_DIR"
    echo ""
    echo "Please run the scAnnex pipeline first to generate results."
    exit 1
fi

print_header
echo ""
print_success "Found results: $RESULTS_DIR"
print_success "Found data: $(basename $H5AD_FILE)"
echo ""

ENV_NAME="scannex-dashboard"
ENV_FILE="${SCRIPT_DIR}/environment_dashboard.yml"

# Check if environment exists, create if not
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo ""
    print_info "Dashboard environment not detected."
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}  First-time setup: Creating dashboard environment${NC}"
    echo -e "${YELLOW}  This will take approximately 5-10 minutes (one-time only)${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    conda env create -f "$ENV_FILE" -n "$ENV_NAME" || {
        echo ""
        print_error "Failed to create environment"
        echo ""
        echo "Try manually:"
        echo "  cd ${SCRIPT_DIR}"
        echo "  conda env create -f environment_dashboard.yml"
        exit 1
    }
    
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    print_success "Environment setup complete!"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    sleep 1
else
    print_success "Using existing dashboard environment"
fi

# Kill any existing dashboard on port 3838
lsof -ti:3838 | xargs kill -9 2>/dev/null || true
sleep 1

echo ""
print_info "Starting dashboard server..."
echo ""

export SCANNEX_DATA_PATH="$RESULTS_DIR"
cd "${SCRIPT_DIR}"

echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║${NC}  ${GREEN}Dashboard Starting${NC}                                         ${BLUE}║${NC}"
echo -e "${BLUE}╠════════════════════════════════════════════════════════════════╣${NC}"
echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
echo -e "${BLUE}║${NC}  ${CYAN}Access the dashboard at:${NC}                                   ${BLUE}║${NC}"
echo -e "${BLUE}║${NC}  ${YELLOW}http://localhost:3838${NC}                                      ${BLUE}║${NC}"
echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
echo -e "${BLUE}║${NC}  Press ${RED}Ctrl+C${NC} to stop                                        ${BLUE}║${NC}"
echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}Loading R libraries and data (may take 10-20 seconds)...${NC}"
echo ""

# WSL2 fix: Source conda and activate directly instead of using 'conda run'
# This avoids httpuv socket binding issues
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Run dashboard directly in activated environment
# Run dashboard directly in activated environment
# Redirect R messages to log file to keep output clean
LOG_FILE="${SCRIPT_DIR}/dashboard.log"
> "$LOG_FILE"  # Clear old log

Rscript --quiet -e "
suppressPackageStartupMessages({
  library(shiny)
})
cat('\n✓ Loading packages...\n')
cat('✓ Starting dashboard server...\n\n')
shiny::runApp('.', host='0.0.0.0', port=3838, launch.browser=FALSE, quiet=TRUE)
" 2>> "$LOG_FILE" | grep -v -E "Loading required package|Attaching package|The following|masked from|^[[:space:]]*$|^$" || true

