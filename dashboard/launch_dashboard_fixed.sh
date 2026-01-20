#!/bin/bash
#
# Dashboard Launcher - Versión con troubleshooting automático
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-${SCRIPT_DIR}/../results_slc_first_run}"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
RED='\033[0;31m'
NC='\033[0m'

print_header() {
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Interactive Dashboard${NC}                             ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
}

print_success() { echo -e "${GREEN}✓${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }
print_info() { echo -e "${YELLOW}ℹ${NC} $1"; }
print_step() { echo -e "${CYAN}→${NC} $1"; }

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Kill any existing dashboard processes
cleanup_previous() {
    print_step "Checking for previous dashboard instances..."
    
    # Check if port is in use
    if ss -tln 2>/dev/null | grep -q ":3838 " || netstat -tln 2>/dev/null | grep -q ":3838 "; then
        print_info "Port 3838 is in use. Attempting to free it..."
        pkill -f "shiny::runApp" 2>/dev/null || true
        sleep 2
    fi
    
    if ss -tln 2>/dev/null | grep -q ":3838 " || netstat -tln 2>/dev/null | grep -q ":3838 "; then
        print_error "Port 3838 still in use. Run: pkill -f shiny"
        return 1
    fi
    
    print_success "Port 3838 is available"
}

# Detect WSL2 and provide instructions
detect_wsl() {
    if grep -qi microsoft /proc/version 2>/dev/null; then
        WSL_IP=$(hostname -I | awk '{print $1}')
        print_info "WSL2 detected. IP: ${WSL_IP}"
        echo ""
        echo -e "  ${CYAN}Access options:${NC}"
        echo -e "    1. http://localhost:3838 ${GREEN}(try first)${NC}"
        echo -e "    2. http://${WSL_IP}:3838 ${YELLOW}(if localhost fails)${NC}"
        echo ""
        return 0
    fi
    return 1
}

# Main launch
main() {
    print_header
    echo ""
    
    # Check conda environment
    if [[ -z "${CONDA_PREFIX:-}" ]] || [[ "${CONDA_PREFIX}" != *"scannex-dashboard"* ]]; then
        print_step "Activating scannex-dashboard environment..."
        
        # Try to activate conda
        if command_exists conda; then
            eval "$(conda shell.bash hook)"
            conda activate scannex-dashboard || {
                print_error "Failed to activate environment"
                echo ""
                echo "Run setup first:"
                echo "  ./setup_dashboard.sh"
                exit 1
            }
        else
            print_error "Conda not found. Run setup first: ./setup_dashboard.sh"
            exit 1
        fi
    fi
    
    print_success "Environment: scannex-dashboard"
    
    # Check data directory
    if [[ ! -d "$RESULTS_DIR" ]]; then
        print_error "Results directory not found: $RESULTS_DIR"
        exit 1
    fi
    
    H5AD_FILE=$(find "$RESULTS_DIR" -name "*annotated*.h5ad" -type f | head -1)
    if [[ -z "$H5AD_FILE" ]]; then
        print_error "No annotated .h5ad file found in: $RESULTS_DIR"
        exit 1
    fi
    
    print_success "Data found: $(basename $H5AD_FILE)"
    echo ""
    
    # Cleanup previous instances
    cleanup_previous || exit 1
    echo ""
    
    # Set environment variables
    export RETICULATE_PYTHON="$(which python3)"
    export SCANNEX_DATA_PATH="$RESULTS_DIR"
    
    print_success "Configuration complete"
    echo ""
    
    # Detect WSL and show access instructions
    IS_WSL=false
    detect_wsl && IS_WSL=true
    
    print_header
    echo -e "${GREEN}Starting Dashboard...${NC}"
    echo ""
    print_info "The dashboard will start in a few seconds"
    print_info "Press Ctrl+C to stop"
    echo ""
    echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
    
    # Launch in background and monitor
    cd "$SCRIPT_DIR"
    R --vanilla -e "shiny::runApp('.', host='0.0.0.0', port=3838, launch.browser=FALSE)" &
    RPID=$!
    
    # Wait for server to start
    print_step "Waiting for server to start..."
    sleep 5
    
    # Test if accessible
    if curl -s http://localhost:3838 > /dev/null 2>&1; then
        echo ""
        print_success "Dashboard is running!"
        echo ""
        echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
        echo -e "  ${GREEN}Access the dashboard at:${NC}"
        echo -e "  ${CYAN}http://localhost:3838${NC}"
        if $IS_WSL; then
            WSL_IP=$(hostname -I | awk '{print $1}')
            echo -e "  ${CYAN}http://${WSL_IP}:3838${NC} ${YELLOW}(if localhost fails)${NC}"
        fi
        echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
        echo ""
    else
        echo ""
        print_error "Server started but not accessible via localhost"
        echo ""
        print_info "Troubleshooting:"
        echo "  1. Check if R process is running: ps aux | grep shiny"
        echo "  2. Try WSL IP instead of localhost (see above)"
        echo "  3. Check firewall settings"
        echo "  4. Try manual launch (see MANUAL_LAUNCH.md)"
    fi
    
    # Wait for user interrupt
    trap "cleanup" EXIT INT TERM
    wait $RPID
}

cleanup() {
    echo ""
    print_info "Stopping dashboard..."
    kill $RPID 2>/dev/null || true
    pkill -f "shiny::runApp" 2>/dev/null || true
    sleep 1
    print_success "Dashboard stopped"
}

main "$@"
