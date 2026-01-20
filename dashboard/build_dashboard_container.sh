#!/bin/bash
#
# Build script for scAnnex Dashboard Containers
# Supports both Docker and Apptainer/Singularity for HPC environments
#
# Usage:
#   ./build_dashboard_container.sh [docker|apptainer|both]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print functions
print_header() {
    echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}ℹ $1${NC}"
}

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Build Docker container
build_docker() {
    print_header "Building Docker Container"
    
    if ! command_exists docker; then
        print_error "Docker not found. Install Docker first:"
        echo "    https://docs.docker.com/engine/install/"
        return 1
    fi
    
    cd "$SCRIPT_DIR"
    
    print_info "Building image: scannex-dashboard:latest"
    docker build -t scannex-dashboard:latest .
    
    print_success "Docker container built successfully!"
    echo ""
    print_info "Run with:"
    echo "    docker run -d -p 3838:3838 \\"
    echo "      --name scannex-dashboard \\"
    echo "      -v \$(pwd)/../results_slc_first_run:/data \\"
    echo "      scannex-dashboard:latest"
    echo ""
    echo "    Access at: http://localhost:3838"
}

# Build Apptainer/Singularity container
build_apptainer() {
    print_header "Building Apptainer/Singularity Container"
    
    if command_exists apptainer; then
        BUILDER="apptainer"
    elif command_exists singularity; then
        BUILDER="singularity"
    else
        print_error "Neither Apptainer nor Singularity found. Install one of them:"
        echo ""
        echo "  Apptainer (recommended for newer systems):"
        echo "    https://apptainer.org/docs/admin/main/installation.html"
        echo ""
        echo "  Singularity (older systems):"
        echo "    https://sylabs.io/guides/latest/user-guide/quick_start.html"
        return 1
    fi
    
    cd "$SCRIPT_DIR"
    
    print_info "Building with: $BUILDER"
    print_info "Output: scannex-dashboard.sif"
    
    # Check if we can build (requires sudo or --fakeroot)
    if [[ $EUID -eq 0 ]]; then
        # Root user
        $BUILDER build scannex-dashboard.sif scannex-dashboard.def
    else
        # Try with --fakeroot first, fallback to sudo
        print_info "Attempting build with --fakeroot..."
        if $BUILDER build --fakeroot scannex-dashboard.sif scannex-dashboard.def 2>/dev/null; then
            print_success "Built with --fakeroot"
        else
            print_info "Fakeroot failed, trying with sudo..."
            sudo $BUILDER build scannex-dashboard.sif scannex-dashboard.def
        fi
    fi
    
    print_success "Apptainer/Singularity container built successfully!"
    echo ""
    print_info "Run with:"
    echo "    $BUILDER run \\"
    echo "      --bind \$(pwd)/../results_slc_first_run:/data \\"
    echo "      scannex-dashboard.sif"
    echo ""
    echo "    Access at: http://localhost:3838"
    echo ""
    print_info "For HPC with SLURM:"
    echo "    sbatch launch_dashboard.slurm"
}

# Show usage
show_usage() {
    echo "Usage: $0 [docker|apptainer|both]"
    echo ""
    echo "Build containers for the scAnnex interactive dashboard"
    echo ""
    echo "Options:"
    echo "  docker      Build Docker container only"
    echo "  apptainer   Build Apptainer/Singularity container only"
    echo "  both        Build both containers (default)"
    echo ""
    echo "Examples:"
    echo "  $0 docker           # Build Docker container"
    echo "  $0 apptainer        # Build Apptainer container"
    echo "  $0                  # Build both"
}

# Main
main() {
    local mode="${1:-both}"
    
    case "$mode" in
        docker)
            build_docker
            ;;
        apptainer|singularity)
            build_apptainer
            ;;
        both)
            build_docker || print_error "Docker build failed"
            echo ""
            build_apptainer || print_error "Apptainer build failed"
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $mode"
            show_usage
            exit 1
            ;;
    esac
    
    echo ""
    print_header "Build Complete"
    print_info "See dashboard/README.md for detailed usage instructions"
}

main "$@"
