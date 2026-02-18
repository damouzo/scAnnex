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

# Detect execution environment
detect_environment() {
    # Check if running on HPC cluster
    if [[ -n "${SLURM_CLUSTER_NAME:-}" ]] || \
       [[ -n "${SLURM_JOB_ID:-}" ]] || \
       command -v srun &> /dev/null; then
        echo "HPC"
    else
        echo "LOCAL"
    fi
}

ENVIRONMENT=$(detect_environment)

# Colors
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

print_header() {
    if [[ "$ENVIRONMENT" == "HPC" ]]; then
        echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
        echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Dashboard - HPC Environment Detected${NC}           ${BLUE}║${NC}"
        echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
    else
        echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
        echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Interactive Dashboard${NC}                           ${BLUE}║${NC}"
        echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
    fi
}

print_success() { echo -e "${GREEN}✓${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }
print_info() { echo -e "${YELLOW}ℹ${NC} $1"; }
print_warning() { echo -e "${YELLOW}⚠${NC} $1"; }

# Function to find a free port
find_free_port() {
    local PORT=${1:-3838}
    while ss -tln 2>/dev/null | grep -q ":${PORT} " || netstat -tln 2>/dev/null | grep -q ":${PORT} "; do
        PORT=$((PORT + 1))
    done
    echo $PORT
}

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

# HPC Environment Warning and Recommendation
if [[ "$ENVIRONMENT" == "HPC" ]]; then
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}  HPC Environment Detected${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    print_warning "You are running on an HPC cluster"
    echo ""
    echo "For production use on HPC, it is recommended to use:"
    echo ""
    echo -e "  ${GREEN}bash launch_dashboard_hpc.sh${NC}"
    echo ""
    echo "This will:"
    echo "  - Request a SLURM interactive job on a compute node"
    echo "  - Allocate dedicated resources (configurable CPU/RAM)"
    echo "  - Provide SSH tunnel instructions for access"
    echo "  - Follow HPC best practices"
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo "Current script will run on this node (login or compute)."
    echo "This is acceptable for:"
    echo "  - Quick testing"
    echo "  - Small datasets (<10k cells)"
    echo "  - Short sessions"
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    read -p "Continue with current node? [y/N]: " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo ""
        print_info "Exiting. To use HPC mode:"
        echo ""
        echo "  cd ${SCRIPT_DIR}"
        echo "  bash launch_dashboard_hpc.sh $RESULTS_DIR"
        echo ""
        exit 0
    fi
    echo ""
fi

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

# Find available port
DEFAULT_PORT=3838
DASHBOARD_PORT=$(find_free_port $DEFAULT_PORT)

if [[ $DASHBOARD_PORT -ne $DEFAULT_PORT ]]; then
    echo ""
    print_info "Port $DEFAULT_PORT is in use, using port $DASHBOARD_PORT instead"
fi

# Kill any existing dashboard on the selected port (shouldn't happen after port detection)
lsof -ti:$DASHBOARD_PORT | xargs kill -9 2>/dev/null || true
sleep 1

echo ""
print_info "Starting dashboard server..."
echo ""

export SCANNEX_DATA_PATH="$RESULTS_DIR"
cd "${SCRIPT_DIR}"

# Display different messages based on environment
if [[ "$ENVIRONMENT" == "HPC" ]]; then
    HOSTNAME=$(hostname -f)
    
    # Try to detect login node for SSH tunnel instructions
    if [[ "$HOSTNAME" == *"apocrita"* ]] || [[ "$HOSTNAME" == *"qmul"* ]]; then
        LOGIN_NODE="login.hpc.qmul.ac.uk"
    else
        # Generic HPC - try to infer login node
        LOGIN_NODE=$(echo "$HOSTNAME" | sed 's/^[^.]*\./login./')
    fi
    
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  ${GREEN}Dashboard Starting on HPC${NC}                                 ${BLUE}║${NC}"
    echo -e "${BLUE}╠════════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
    echo -e "${BLUE}║${NC}  ${CYAN}Node:${NC} $HOSTNAME"
    echo -e "${BLUE}║${NC}  ${CYAN}Port:${NC} $DASHBOARD_PORT"
    echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}  SSH TUNNEL SETUP (Required)${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo "Open a NEW terminal on your LOCAL machine and run:"
    echo ""
    echo -e "  ${YELLOW}ssh -N -L ${DASHBOARD_PORT}:${HOSTNAME}:${DASHBOARD_PORT} ${USER}@${LOGIN_NODE}${NC}"
    echo ""
    echo "Keep this tunnel running while using the dashboard."
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}  ACCESS DASHBOARD${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo "After setting up the SSH tunnel, open your browser:"
    echo ""
    echo -e "  ${GREEN}http://localhost:${DASHBOARD_PORT}${NC}"
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo -e "${RED}Press Ctrl+C in THIS terminal to stop the dashboard${NC}"
    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
else
    # Local environment - simpler message
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  ${GREEN}Dashboard Starting${NC}                                         ${BLUE}║${NC}"
    echo -e "${BLUE}╠════════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
    echo -e "${BLUE}║${NC}  ${CYAN}Access the dashboard at:${NC}                                   ${BLUE}║${NC}"
    echo -e "${BLUE}║${NC}  ${YELLOW}http://localhost:${DASHBOARD_PORT}${NC}"
    echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
    echo -e "${BLUE}║${NC}  Press ${RED}Ctrl+C${NC} to stop                                        ${BLUE}║${NC}"
    echo -e "${BLUE}║${NC}                                                                ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
fi
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
shiny::runApp('.', host='0.0.0.0', port=${DASHBOARD_PORT}, launch.browser=FALSE, quiet=TRUE)
" 2>> "$LOG_FILE" | grep -v -E "Loading required package|Attaching package|The following|masked from|^[[:space:]]*$|^$" || true

