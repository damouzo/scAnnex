#!/bin/bash
#
# Simple scAnnex Dashboard Launcher (Conda only)
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-${SCRIPT_DIR}/../results}"

# Convert relative path to absolute path
if [[ ! "${RESULTS_DIR}" = /* ]]; then
    RESULTS_DIR="$(cd "${SCRIPT_DIR}" && cd "${RESULTS_DIR}" && pwd)"
fi

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m'

echo ""
echo -e "${CYAN}═══════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}  scAnnex Interactive Dashboard${NC}"
echo -e "${CYAN}═══════════════════════════════════════════════════════════${NC}"
echo ""

# Validate results directory
if [[ ! -d "$RESULTS_DIR" ]]; then
    echo -e "${RED}✗ Error:${NC} Results directory not found: $RESULTS_DIR"
    echo "Usage: $0 [path/to/results_directory]"
    exit 1
fi

# Find h5ad file
H5AD_FILE=$(find "$RESULTS_DIR" -name "*annotated*.h5ad" -type f | head -1)
if [[ -z "$H5AD_FILE" ]]; then
    echo -e "${RED}✗ Error:${NC} No annotated .h5ad file found in: $RESULTS_DIR"
    exit 1
fi

echo -e "${GREEN}✓${NC} Found results: $RESULTS_DIR"
echo -e "${GREEN}✓${NC} Found data: $(basename $H5AD_FILE)"
echo ""

# Check port
PORT=3838
if ss -tln 2>/dev/null | grep -q ":${PORT} " || netstat -tln 2>/dev/null | grep -q ":${PORT} "; then
    echo -e "${CYAN}ℹ${NC} Port ${PORT} is in use, trying next port..."
    PORT=$((PORT + 1))
fi

echo -e "${GREEN}✓${NC} Starting dashboard on port ${PORT}..."
echo ""
echo -e "  ${CYAN}Dashboard URL:${NC} http://localhost:${PORT}"
echo ""
echo -e "  Press ${CYAN}Ctrl+C${NC} to stop the dashboard"
echo ""
echo -e "  ${CYAN}⏳ Loading R libraries and data...${NC} (this may take 10-20 seconds)"
echo ""
echo -e "${CYAN}═══════════════════════════════════════════════════════════${NC}"
echo ""

# Export data path
export SCANNEX_DATA_PATH="$RESULTS_DIR"

# Launch with conda
cd "${SCRIPT_DIR}"

# Create a simple wrapper to filter output
{
    conda run -n scannex-dashboard --no-capture-output R --quiet --no-save -e "shiny::runApp('.', host='0.0.0.0', port=${PORT})" 2>&1
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
