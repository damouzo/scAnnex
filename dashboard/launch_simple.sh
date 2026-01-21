#!/bin/bash
# Simple dashboard launcher for testing

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-${SCRIPT_DIR}/../results}"

# Colors
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m'

ENV_NAME="scannex-dashboard"

# Check if environment exists, create if not
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${CYAN}Creating dashboard environment (first time setup, ~5-10 min)...${NC}"
    conda env create -f "${SCRIPT_DIR}/environment_dashboard.yml" -n "$ENV_NAME" || exit 1
    echo -e "${GREEN}✓ Environment created${NC}"
fi

# Kill any existing dashboard on port 3838
lsof -ti:3838 | xargs kill -9 2>/dev/null || true

export SCANNEX_DATA_PATH="$RESULTS_DIR"
cd "${SCRIPT_DIR}"

echo ""
echo -e "${GREEN}Starting dashboard...${NC}"
echo -e "Access at: ${CYAN}http://localhost:3838${NC}"
echo "Press Ctrl+C to stop"
echo ""

# Run dashboard
conda run -n "$ENV_NAME" R --quiet -e "shiny::runApp('.', host='0.0.0.0', port=3838, launch.browser=FALSE)"
