#!/bin/bash
#
# scAnnex Dashboard Launcher
# One-command solution to launch the interactive dashboard
#
# Usage: ./launch_dashboard.sh [path/to/results_directory]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-${SCRIPT_DIR}/../results_slc_first_run}"

# Convert relative path to absolute path
if [[ ! "${RESULTS_DIR}" = /* ]]; then
    RESULTS_DIR="$(cd "${SCRIPT_DIR}" && cd "${RESULTS_DIR}" && pwd)"
fi

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

print_header() {
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Interactive Dashboard Launcher${NC}                    ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
}

print_success() { echo -e "${GREEN}✓${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }
print_info() { echo -e "${YELLOW}ℹ${NC} $1"; }
print_step() { echo -e "${CYAN}→${NC} $1"; }

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

#==============================================================================
# Detection: Find the best available method
#==============================================================================

detect_method() {
    # Priority order: Conda+R > Apptainer > Singularity > Docker > Manual install
    # (Conda first because it's lightweight and likely already set up)
    
    if command_exists conda && [[ -f "${SCRIPT_DIR}/environment_dashboard.yml" ]]; then
        # Check if environment exists
        if conda env list | grep -q "scannex-dashboard"; then
            echo "conda"
            return
        fi
    fi
    
    if command_exists apptainer; then
        echo "apptainer"
    elif command_exists singularity; then
        echo "singularity"
    elif command_exists docker; then
        echo "docker"
    else
        echo "manual"
    fi
}

#==============================================================================
# Method 1: Apptainer (HPC-friendly, no sudo needed)
#==============================================================================

launch_apptainer() {
    local TOOL="$1"
    print_step "Launching with ${TOOL}..."
    
    # Container URL (will be hosted on GitHub Container Registry or Sylabs)
    CONTAINER_URL="docker://ghcr.io/yourusername/scannex-dashboard:latest"
    CONTAINER_SIF="${SCRIPT_DIR}/scannex-dashboard.sif"
    
    # Check if container exists locally
    if [[ ! -f "$CONTAINER_SIF" ]]; then
        print_info "Container not found locally. Pulling from registry..."
        print_info "This is a one-time download (~500MB, may take a few minutes)"
        
        ${TOOL} pull "$CONTAINER_SIF" "$CONTAINER_URL" || {
            print_error "Failed to pull container. Building locally instead..."
            if [[ -f "${SCRIPT_DIR}/scannex-dashboard.def" ]]; then
                ${TOOL} build --fakeroot "$CONTAINER_SIF" "${SCRIPT_DIR}/scannex-dashboard.def" || {
                    print_error "Build failed. Try: sudo ${TOOL} build ..."
                    return 1
                }
            else
                print_error "No container definition found"
                return 1
            fi
        }
    else
        print_success "Found existing container: $CONTAINER_SIF"
    fi
    
    # Find available port
    PORT=3838
    while ss -tln 2>/dev/null | grep -q ":${PORT} " || netstat -tln 2>/dev/null | grep -q ":${PORT} "; do
        PORT=$((PORT + 1))
    done
    
    print_success "Using port: ${PORT}"
    print_info "Starting dashboard..."
    echo ""
    print_header
    echo -e "${GREEN}Dashboard is starting!${NC}"
    echo ""
    echo -e "Access it at: ${CYAN}http://localhost:${PORT}${NC}"
    echo ""
    echo "Press Ctrl+C to stop"
    echo ""
    
    # Launch with instance (runs in background)
    ${TOOL} instance start \
        --bind "${RESULTS_DIR}:/data" \
        --env SCANNEX_DATA_PATH=/data \
        "$CONTAINER_SIF" \
        scannex-dashboard-instance
    
    # Run shiny server in the instance
    ${TOOL} exec instance://scannex-dashboard-instance \
        R -e "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=${PORT})" &
    
    SHINY_PID=$!
    
    # Wait for user interrupt
    trap "cleanup_apptainer ${TOOL}" EXIT INT TERM
    wait $SHINY_PID
}

cleanup_apptainer() {
    local TOOL="$1"
    echo ""
    print_info "Stopping dashboard..."
    ${TOOL} instance stop scannex-dashboard-instance 2>/dev/null || true
    print_success "Dashboard stopped"
}

#==============================================================================
# Method 2: Docker
#==============================================================================

launch_docker() {
    print_step "Launching with Docker..."
    
    CONTAINER_NAME="scannex-dashboard"
    IMAGE_NAME="ghcr.io/yourusername/scannex-dashboard:latest"
    
    # Try to pull pre-built image
    if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
        print_info "Pulling pre-built container..."
        docker pull "$IMAGE_NAME" || {
            print_info "Pull failed. Building locally..."
            docker build -t "$IMAGE_NAME" "${SCRIPT_DIR}" || return 1
        }
    fi
    
    # Stop existing container if running
    docker stop "$CONTAINER_NAME" 2>/dev/null || true
    docker rm "$CONTAINER_NAME" 2>/dev/null || true
    
    PORT=3838
    
    print_success "Starting dashboard..."
    docker run -d \
        --name "$CONTAINER_NAME" \
        -p "${PORT}:3838" \
        -v "${RESULTS_DIR}:/data" \
        "$IMAGE_NAME"
    
    echo ""
    print_header
    echo -e "${GREEN}Dashboard is running!${NC}"
    echo ""
    echo -e "Access it at: ${CYAN}http://localhost:${PORT}${NC}"
    echo ""
    echo "To stop: docker stop ${CONTAINER_NAME}"
    echo ""
}

#==============================================================================
# Method 3: Conda environment (lightweight, no container)
#==============================================================================

launch_conda() {
    print_step "Launching with Conda environment..."
    
    ENV_NAME="scannex-dashboard"
    ENV_FILE="${SCRIPT_DIR}/environment_dashboard.yml"
    
    # Check if environment exists
    if ! conda env list | grep -q "^${ENV_NAME} "; then
        print_info "Creating Conda environment (one-time setup)..."
        conda env create -f "$ENV_FILE" -n "$ENV_NAME" || return 1
    fi
    
    print_success "Activating environment..."
    
    PORT=3838
    
    # Check if port is already in use
    if ss -tln 2>/dev/null | grep -q ":${PORT} " || netstat -tln 2>/dev/null | grep -q ":${PORT} "; then
        print_info "Port ${PORT} is in use, trying next port..."
        PORT=$((PORT + 1))
    fi
    
    print_success "Starting dashboard on port ${PORT}..."
    echo ""
    print_header
    echo -e "${GREEN}Dashboard is starting!${NC}"
    echo ""
    echo -e "Access it at: ${CYAN}http://localhost:${PORT}${NC}"
    echo ""
    echo "Press Ctrl+C to stop"
    echo ""
    echo -e "${YELLOW}Loading R libraries and data... (this may take 10-20 seconds)${NC}"
    echo ""
    
    export SCANNEX_DATA_PATH="$RESULTS_DIR"
    
    # Activate and run using conda run
    cd "${SCRIPT_DIR}"
    
    # Filter output to show only important messages
    {
        conda run -n "$ENV_NAME" --no-capture-output R --quiet --no-save -e "shiny::runApp('.', host='0.0.0.0', port=${PORT})" 2>&1
    } | while IFS= read -r line; do
        # Skip known harmless warnings
        if echo "$line" | grep -q "UserWarning: Signature"; then
            continue
        elif echo "$line" | grep -q "This warnings indicates broken support"; then
            continue
        # Show important messages
        elif echo "$line" | grep -q "Listening on"; then
            echo ""
            echo -e "${GREEN}✓ Dashboard ready!${NC}"
            echo "$line"
            echo ""
        elif echo "$line" | grep -qE "(Loading required package|scAnnex Dashboard|Auto-detected|Default data path|initialized)"; then
            echo "  $line"
        elif echo "$line" | grep -qE "(Error|error|Failed|failed)" && ! echo "$line" | grep -q "UserWarning"; then
            echo -e "${RED}✗${NC} $line"
        fi
    done
}

#==============================================================================
# Method 4: Manual installation guide
#==============================================================================

show_manual_instructions() {
    print_error "No container system or Conda environment found"
    echo ""
    echo "Please install one of the following:"
    echo ""
    echo "1. ${CYAN}Apptainer${NC} (Recommended for HPC):"
    echo "   https://apptainer.org/docs/admin/main/installation.html"
    echo ""
    echo "2. ${CYAN}Docker${NC} (Recommended for local development):"
    echo "   https://docs.docker.com/engine/install/"
    echo ""
    echo "3. ${CYAN}Conda/Mamba${NC} (Lightweight, no container):"
    echo "   https://github.com/conda-forge/miniforge"
    echo "   Then run: conda env create -f dashboard/environment_dashboard.yml"
    echo ""
    echo "After installation, run this script again!"
}

#==============================================================================
# Main
#==============================================================================

main() {
    print_header
    echo ""
    
    # Validate results directory
    if [[ ! -d "$RESULTS_DIR" ]]; then
        print_error "Results directory not found: $RESULTS_DIR"
        echo "Usage: $0 [path/to/results_directory]"
        exit 1
    fi
    
    # Find h5ad file
    H5AD_FILE=$(find "$RESULTS_DIR" -name "*annotated*.h5ad" -type f | head -1)
    if [[ -z "$H5AD_FILE" ]]; then
        print_error "No annotated .h5ad file found in: $RESULTS_DIR"
        exit 1
    fi
    
    print_success "Found results: $RESULTS_DIR"
    print_success "Found data: $(basename $H5AD_FILE)"
    echo ""
    
    # Detect and launch
    print_step "Detecting available container/environment systems..."
    METHOD=$(detect_method)
    
    case "$METHOD" in
        apptainer)
            print_success "Using Apptainer (HPC-optimized)"
            launch_apptainer apptainer
            ;;
        singularity)
            print_success "Using Singularity"
            launch_apptainer singularity
            ;;
        docker)
            print_success "Using Docker"
            launch_docker
            ;;
        conda)
            print_success "Using Conda environment"
            launch_conda
            ;;
        manual)
            show_manual_instructions
            exit 1
            ;;
    esac
}

main "$@"
