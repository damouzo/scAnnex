#!/bin/bash
#
# scAnnex Dashboard HPC Launcher
# Launches dashboard on HPC compute nodes via SLURM interactive session
#
# Usage: ./launch_dashboard_hpc.sh [OPTIONS] [results_directory]
#
# Options:
#   --cpus N          Number of CPUs (default: 4)
#   --mem SIZE        Memory allocation (default: 8G)
#   --time DURATION   Time limit (default: 4:00:00)
#   --port PORT       Dashboard port (default: 3838, auto-increment if busy)
#   --partition NAME  SLURM partition (default: compute)
#   --help            Show this help message
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
DEFAULT_CPUS=4
DEFAULT_MEM="8G"
DEFAULT_TIME="4:00:00"
DEFAULT_PORT=3838
DEFAULT_PARTITION="compute"

# Parse arguments
CPUS=$DEFAULT_CPUS
MEM=$DEFAULT_MEM
TIME_LIMIT=$DEFAULT_TIME
PORT=$DEFAULT_PORT
PARTITION=$DEFAULT_PARTITION
RESULTS_DIR=""

show_help() {
    cat << EOF
scAnnex Dashboard HPC Launcher

Launches the scAnnex interactive dashboard on HPC compute nodes
via SLURM interactive session.

USAGE:
    $0 [OPTIONS] [results_directory]

OPTIONS:
    --cpus N          Number of CPUs (default: $DEFAULT_CPUS)
    --mem SIZE        Memory allocation (default: $DEFAULT_MEM)
                      Examples: 8G, 16G, 32G
    --time DURATION   Time limit (default: $DEFAULT_TIME)
                      Format: HH:MM:SS or D-HH:MM:SS
                      Examples: 2:00:00, 8:00:00, 1-00:00:00
    --port PORT       Dashboard port (default: $DEFAULT_PORT)
                      Will auto-increment if port is busy
    --partition NAME  SLURM partition (default: $DEFAULT_PARTITION)
                      Options: compute, computeshort, highmem, etc.
    --help            Show this help message

EXAMPLES:
    # Basic usage with default settings (4 CPU, 8G RAM, 4 hours)
    $0 ../results

    # Large dataset with more resources
    $0 --cpus 8 --mem 16G --time 8:00:00 ../results

    # Quick session on short partition (max 1 hour)
    $0 --partition computeshort --time 1:00:00 ../results

    # Custom port
    $0 --port 3850 ../results

RESOURCE RECOMMENDATIONS:
    Dataset Size        CPUs    Memory    Time
    -------------------------------------------------
    < 10k cells         2-4     4-8G      2-4h
    10k-50k cells       4-8     8-16G     4-8h
    50k-100k cells      8-16    16-32G    8-12h
    > 100k cells        16+     32G+      12-24h

NOTES:
    - Dashboard will run on a dedicated compute node
    - You will need to set up an SSH tunnel to access it
    - Press Ctrl+C to stop the dashboard and release resources
    - Session will automatically end after time limit expires

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --time)
            TIME_LIMIT="$2"
            shift 2
            ;;
        --port)
            PORT="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --help)
            show_help
            exit 0
            ;;
        -*)
            echo "Error: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
        *)
            RESULTS_DIR="$1"
            shift
            ;;
    esac
done

# Set default results directory if not provided
if [[ -z "$RESULTS_DIR" ]]; then
    RESULTS_DIR="${SCRIPT_DIR}/../results"
fi

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
    echo -e "${BLUE}║${NC}  ${CYAN}scAnnex Dashboard - HPC Mode${NC}                            ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
}

print_success() { echo -e "${GREEN}✓${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }
print_info() { echo -e "${YELLOW}ℹ${NC} $1"; }

# Validate environment
if ! command -v srun &> /dev/null; then
    print_header
    echo ""
    print_error "SLURM not detected (srun command not found)"
    echo ""
    echo "This script is designed for HPC clusters with SLURM scheduler."
    echo "For local execution, use: bash launch_dashboard.sh"
    echo ""
    exit 1
fi

if ! command -v conda &> /dev/null; then
    print_header
    echo ""
    print_error "Conda not found"
    echo ""
    echo "Please ensure conda/mamba is available on compute nodes."
    exit 1
fi

# Validate results directory
if [[ ! -d "$RESULTS_DIR" ]]; then
    print_header
    echo ""
    print_error "Results directory not found: $RESULTS_DIR"
    echo ""
    echo "Usage: $0 [OPTIONS] <results_directory>"
    echo "Use --help for more information"
    exit 1
fi

# Find h5ad file
H5AD_FILE=$(find "$RESULTS_DIR" -name "*annotated*.h5ad" -type f 2>/dev/null | head -1)
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
print_success "Results directory: $RESULTS_DIR"
print_success "H5AD file: $(basename $H5AD_FILE)"
echo ""

ENV_NAME="scannex-dashboard"
ENV_FILE="${SCRIPT_DIR}/environment_dashboard.yml"

# Check if environment exists, create if not
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo ""
    print_info "Dashboard environment not detected"
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
    print_success "Environment setup complete"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
else
    print_success "Dashboard environment ready"
fi

# Display resource request
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}  Requesting SLURM Interactive Session${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "  Partition:  $PARTITION"
echo "  CPUs:       $CPUS"
echo "  Memory:     $MEM"
echo "  Time limit: $TIME_LIMIT"
echo "  Port:       $PORT (will auto-increment if busy)"
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
print_info "Waiting for compute node allocation..."
print_info "This may take a few moments depending on cluster load"
echo ""

# Launch interactive session with dashboard
# Note: We use srun --pty to get an interactive shell on compute node
srun \
    --partition="$PARTITION" \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME_LIMIT" \
    --pty \
    bash -c "
    set -euo pipefail
    
    # Define colors for use inside the job
    GREEN='\033[0;32m'
    CYAN='\033[0;36m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    RED='\033[0;31m'
    NC='\033[0m'
    
    # Get compute node hostname
    COMPUTE_NODE=\$(hostname -f)
    
    # Function to find free port
    find_free_port() {
        local PORT=$PORT
        while ss -tln 2>/dev/null | grep -q \":\${PORT} \" || netstat -tln 2>/dev/null | grep -q \":\${PORT} \"; do
            PORT=\$((PORT + 1))
        done
        echo \$PORT
    }
    
    # Find available port
    DASHBOARD_PORT=\$(find_free_port)
    
    echo \"\"
    echo -e \"\${BLUE}╔════════════════════════════════════════════════════════════════╗\${NC}\"
    echo -e \"\${BLUE}║\${NC}  \${GREEN}Compute Node Allocated\${NC}                                     \${BLUE}║\${NC}\"
    echo -e \"\${BLUE}╚════════════════════════════════════════════════════════════════╝\${NC}\"
    echo \"\"
    echo -e \"\${GREEN}✓\${NC} Node:       \${COMPUTE_NODE}\"
    echo -e \"\${GREEN}✓\${NC} Port:       \${DASHBOARD_PORT}\"
    echo -e \"\${GREEN}✓\${NC} Resources:  ${CPUS} CPUs, ${MEM} RAM\"
    echo -e \"\${GREEN}✓\${NC} Time limit: ${TIME_LIMIT}\"
    echo \"\"
    
    # Determine login node (HPC-specific)
    if [[ \"\${COMPUTE_NODE}\" == *\"apocrita\"* ]] || [[ \"\${COMPUTE_NODE}\" == *\"qmul\"* ]]; then
        LOGIN_NODE=\"login.hpc.qmul.ac.uk\"
    else
        # Generic HPC - try to infer login node
        LOGIN_NODE=\$(echo \"\${COMPUTE_NODE}\" | sed 's/^[^.]*\./login./')
    fi
    
    echo -e \"\${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\${NC}\"
    echo -e \"\${YELLOW}SSH Tunnel Required\${NC}\"
    echo -e \"\${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\${NC}\"
    echo \"\"
    echo -e \"\${YELLOW}VSCode Remote SSH:\${NC}\"
    echo \"  Open new VSCode terminal and run:\"
    echo -e \"  \${GREEN}ssh -N -L \${DASHBOARD_PORT}:localhost:\${DASHBOARD_PORT} \${COMPUTE_NODE}\${NC}\"
    echo \"\"
    echo -e \"\${YELLOW}Local terminal (MobaXterm/PowerShell/Mac):\${NC}\"
    echo \"  Open new terminal on your local machine and run:\"
    echo -e \"  \${GREEN}ssh -N -L \${DASHBOARD_PORT}:\${COMPUTE_NODE}:\${DASHBOARD_PORT} \\\$USER@\${LOGIN_NODE}\${NC}\"
    echo \"\"
    echo \"Then open browser: http://localhost:\${DASHBOARD_PORT}\"
    echo \"Keep tunnel running. Press Ctrl+C here to stop dashboard.\"
    echo \"\"
    echo -e \"\${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\${NC}\"
    echo \"\"
    
    # Source conda
    source \"\$(conda info --base)/etc/profile.d/conda.sh\"
    conda activate \"$ENV_NAME\"
    
    # Set environment variable for data path
    export SCANNEX_DATA_PATH=\"$RESULTS_DIR\"
    
    # Change to dashboard directory
    cd \"${SCRIPT_DIR}\"
    
    # Launch dashboard
    echo -e \"\${YELLOW}Starting dashboard...\${NC}\"
    echo \"\"
    
    Rscript --quiet -e \"
    suppressPackageStartupMessages({
      library(shiny)
    })
    cat('✓ Ready\\\\n\\\\n')
    shiny::runApp('.', host='0.0.0.0', port=\${DASHBOARD_PORT}, launch.browser=FALSE, quiet=FALSE)
    \"
    
    # This will only execute if dashboard exits cleanly
    echo \"\"
    echo -e \"\${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\${NC}\"
    echo -e \"\${GREEN}Dashboard stopped\${NC}\"
    echo -e \"\${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\${NC}\"
    echo \"\"
"

# This will only execute if srun exits (normally or due to error)
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}SLURM interactive session ended${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "Session ended. Possible reasons:"
echo "  - Time limit reached (${TIME_LIMIT})"
echo "  - You pressed Ctrl+C"
echo "  - SLURM preempted the job"
echo ""
echo "To relaunch the dashboard, run:"
echo "  cd ${SCRIPT_DIR}"
echo "  bash launch_dashboard_hpc.sh $RESULTS_DIR"
echo ""
