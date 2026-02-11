#!/bin/bash
#
# scAnnex Dashboard - Background Launch Wrapper
# Launches dashboard in true daemon mode, detached from parent process
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-}"

if [[ -z "$RESULTS_DIR" ]]; then
    echo "Error: Results directory required"
    echo "Usage: $0 <results_directory>"
    exit 1
fi

# Convert to absolute path
if [[ ! "${RESULTS_DIR}" = /* ]]; then
    RESULTS_DIR="$(cd "${SCRIPT_DIR}" && cd "${RESULTS_DIR}" && pwd)"
fi

LOG_FILE="${RESULTS_DIR}/dashboard_launch.log"
PID_FILE="${RESULTS_DIR}/dashboard.pid"

# Kill any existing dashboard
lsof -ti:3838 | xargs kill -9 2>/dev/null || true
sleep 1

# Launch dashboard with proper daemonization
# Use setsid to create new session, completely detached
setsid bash -c "
    cd '${SCRIPT_DIR}'
    export SCANNEX_DATA_PATH='$RESULTS_DIR'
    
    # Source conda
    source \"\$(conda info --base)/etc/profile.d/conda.sh\"
    conda activate scannex-dashboard
    
    # Run dashboard server with explicit keep-alive
    exec Rscript --vanilla -e \"
    suppressPackageStartupMessages({
      library(shiny)
    })
    cat('Starting dashboard server for:', Sys.getenv('SCANNEX_DATA_PATH'), '\\\n')
    shiny::runApp('.', host='0.0.0.0', port=3838, launch.browser=FALSE, quiet=FALSE)
    \" < /dev/null
" > "$LOG_FILE" 2>&1 &

DASH_PID=$!

# Save PID
echo $DASH_PID > "$PID_FILE"

# Wait a moment to see if it starts
sleep 3

if ps -p $DASH_PID > /dev/null 2>&1; then
    echo "Dashboard launched successfully (PID: $DASH_PID)"
    echo "Log file: $LOG_FILE"
    exit 0
else
    echo "Error: Dashboard failed to start"
    echo "Check log: $LOG_FILE"
    exit 1
fi
