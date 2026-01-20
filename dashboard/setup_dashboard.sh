#!/bin/bash
#
# scAnnex Dashboard - One-Time Setup
# This script sets up the dashboard environment (run once after cloning)
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

print_header() {
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Dashboard Setup${NC}                                  ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
}

print_success() { echo -e "${GREEN}✓${NC} $1"; }
print_info() { echo -e "${YELLOW}ℹ${NC} $1"; }
print_step() { echo -e "${CYAN}→${NC} $1"; }

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_header
echo ""

#==============================================================================
# Option 1: Pull pre-built container (fastest, no build needed)
#==============================================================================

setup_prebuilt_container() {
    local tool="$1"
    print_step "Downloading pre-built container..."
    
    # GitHub Container Registry URL (update with your username)
    CONTAINER_URL="docker://ghcr.io/YOURUSERNAME/scannex-dashboard:latest"
    CONTAINER_SIF="${SCRIPT_DIR}/scannex-dashboard.sif"
    
    print_info "This is a one-time download (~500MB)"
    
    ${tool} pull "$CONTAINER_SIF" "$CONTAINER_URL" && {
        print_success "Container ready!"
        echo ""
        print_info "Next: Run ./launch_dashboard.sh to start"
        return 0
    }
    
    print_info "Pre-built container not available yet. Building locally..."
    return 1
}

#==============================================================================
# Option 2: Build container locally (requires build once)
#==============================================================================

build_container() {
    local tool="$1"
    print_step "Building container locally..."
    
    CONTAINER_SIF="${SCRIPT_DIR}/scannex-dashboard.sif"
    
    print_info "This may take 10-15 minutes (one-time build)"
    
    ${tool} build --fakeroot "$CONTAINER_SIF" "${SCRIPT_DIR}/scannex-dashboard.def" || {
        print_info "Fakeroot build failed. Trying with sudo..."
        sudo ${tool} build "$CONTAINER_SIF" "${SCRIPT_DIR}/scannex-dashboard.def" || {
            print_error "Build failed"
            return 1
        }
    }
    
    print_success "Container built successfully!"
    echo ""
    print_info "Next: Run ./launch_dashboard.sh to start"
}

#==============================================================================
# Option 3: Conda environment (lightweight, no container)
#==============================================================================

setup_conda_env() {
    print_step "Setting up Conda environment..."
    
    ENV_NAME="scannex-dashboard"
    ENV_FILE="${SCRIPT_DIR}/environment_dashboard.yml"
    
    if conda env list | grep -q "^${ENV_NAME} "; then
        print_info "Environment already exists"
        read -p "Recreate? [y/N] " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            conda env remove -n "$ENV_NAME" -y
        else
            print_success "Using existing environment"
            echo ""
            print_info "Next: Run ./launch_dashboard.sh to start"
            return 0
        fi
    fi
    
    print_info "Creating environment (may take 5-10 minutes)..."
    conda env create -f "$ENV_FILE" -n "$ENV_NAME"
    
    print_success "Conda environment ready!"
    echo ""
    print_info "Next: Run ./launch_dashboard.sh to start"
}

#==============================================================================
# Main Setup Flow
#==============================================================================

main() {
    # Detect what's available
    if command_exists apptainer; then
        print_success "Found Apptainer"
        echo ""
        echo "Choose setup method:"
        echo "  1) Download pre-built container (fastest, ~500MB)"
        echo "  2) Build container locally (slower, but works offline)"
        echo ""
        read -p "Choice [1]: " choice
        choice=${choice:-1}
        
        case "$choice" in
            1)
                setup_prebuilt_container apptainer || build_container apptainer
                ;;
            2)
                build_container apptainer
                ;;
        esac
        
    elif command_exists singularity; then
        print_success "Found Singularity"
        setup_prebuilt_container singularity || build_container singularity
        
    elif command_exists docker; then
        print_success "Found Docker"
        print_info "Docker will pull pre-built image automatically on first launch"
        print_success "Setup complete!"
        echo ""
        print_info "Next: Run ./launch_dashboard.sh to start"
        
    elif command_exists conda; then
        print_success "Found Conda"
        setup_conda_env
        
    else
        print_info "No container system or Conda found"
        echo ""
        echo "Please install one of:"
        echo ""
        echo "  1. ${CYAN}Conda/Mamba${NC} (Recommended - easiest, no sudo needed):"
        echo "     curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
        echo "     bash Miniforge3-Linux-x86_64.sh"
        echo ""
        echo "  2. ${CYAN}Apptainer${NC} (HPC environments):"
        echo "     https://apptainer.org/docs/admin/main/installation.html"
        echo ""
        echo "  3. ${CYAN}Docker${NC} (Local development):"
        echo "     https://docs.docker.com/engine/install/"
        echo ""
        echo "After installation, run this script again!"
        exit 1
    fi
}

main "$@"
