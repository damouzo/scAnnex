#!/bin/bash
################################################################################
# Environment Verification Script
# 
# Tests that the conda environment is correctly set up and ready for pipeline
################################################################################

set -euo pipefail

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}  scAnnex Environment Verification${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo

# Set paths
export PATH="/home/damo/miniforge3/bin:$PATH"
PROJECT_DIR="/home/damo/scAnnex"
cd "${PROJECT_DIR}"

PASS=0
FAIL=0

# Test 1: Check conda/mamba installation
echo -e "${BLUE}[1/5] Checking Conda/Mamba installation...${NC}"
if command -v mamba &> /dev/null; then
    VERSION=$(mamba --version | head -1)
    echo -e "  ${GREEN}✓ Mamba found: ${VERSION}${NC}"
    ((PASS++))
else
    echo -e "  ${RED}✗ Mamba not found${NC}"
    ((FAIL++))
fi

# Test 2: Check if scannex environment exists
echo -e "${BLUE}[2/5] Checking scannex environment...${NC}"
if conda env list | grep -q "^scannex "; then
    echo -e "  ${GREEN}✓ Environment 'scannex' exists${NC}"
    ((PASS++))
else
    echo -e "  ${RED}✗ Environment 'scannex' not found${NC}"
    echo -e "  ${YELLOW}  Run: mamba env create -f env/scanpy.yml${NC}"
    ((FAIL++))
fi

# Test 3: Test environment activation and package imports
echo -e "${BLUE}[3/5] Testing package imports...${NC}"
source /home/damo/miniforge3/bin/activate scannex

python << 'PYEOF'
import sys
errors = []

try:
    import scanpy as sc
    print(f"  ✓ Scanpy {sc.__version__}")
except ImportError as e:
    errors.append(f"  ✗ Scanpy: {e}")

try:
    import anndata as ad
    print(f"  ✓ AnnData {ad.__version__}")
except ImportError as e:
    errors.append(f"  ✗ AnnData: {e}")

try:
    import celltypist as ct
    print(f"  ✓ CellTypist {ct.__version__}")
except ImportError as e:
    errors.append(f"  ✗ CellTypist: {e}")

try:
    import scrublet
    print(f"  ✓ Scrublet imported")
except ImportError as e:
    errors.append(f"  ✗ Scrublet: {e}")

try:
    import harmonypy
    print(f"  ✓ HarmonyPy imported")
except ImportError as e:
    errors.append(f"  ✗ HarmonyPy: {e}")

if errors:
    for err in errors:
        print(err)
    sys.exit(1)
PYEOF

if [ $? -eq 0 ]; then
    echo -e "  ${GREEN}✓ All packages imported successfully${NC}"
    ((PASS++))
else
    echo -e "  ${RED}✗ Package import failed${NC}"
    ((FAIL++))
fi

# Test 4: Check Nextflow installation
echo -e "${BLUE}[4/5] Checking Nextflow...${NC}"
if command -v nextflow &> /dev/null; then
    NF_VERSION=$(nextflow -version 2>&1 | grep "version" | awk '{print $3}')
    echo -e "  ${GREEN}✓ Nextflow found: ${NF_VERSION}${NC}"
    ((PASS++))
else
    echo -e "  ${RED}✗ Nextflow not found${NC}"
    ((FAIL++))
fi

# Test 5: Check demo data
echo -e "${BLUE}[5/5] Checking demo data...${NC}"
if [ -f "data_demo/H5AD/samplesheet.csv" ]; then
    echo -e "  ${GREEN}✓ Samplesheet found${NC}"
    if [ -d "data_demo/10xMTX/filtered_feature_bc_matrix" ]; then
        MATRIX=$(ls data_demo/10xMTX/filtered_feature_bc_matrix/matrix.mtx.gz 2>/dev/null)
        if [ -n "$MATRIX" ]; then
            SIZE=$(du -h "$MATRIX" | cut -f1)
            echo -e "  ${GREEN}✓ MTX data found (${SIZE})${NC}"
            ((PASS++))
        else
            echo -e "  ${RED}✗ MTX matrix not found${NC}"
            ((FAIL++))
        fi
    else
        echo -e "  ${RED}✗ MTX directory not found${NC}"
        ((FAIL++))
    fi
else
    echo -e "  ${RED}✗ Samplesheet not found${NC}"
    ((FAIL++))
fi

# Summary
echo
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}  Verification Summary${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo -e "  ${GREEN}Passed: ${PASS}/5${NC}"
if [ ${FAIL} -gt 0 ]; then
    echo -e "  ${RED}Failed: ${FAIL}/5${NC}"
fi
echo

if [ ${FAIL} -eq 0 ]; then
    echo -e "${GREEN}✓✓✓ ALL CHECKS PASSED - READY TO RUN PIPELINE ✓✓✓${NC}"
    echo
    echo -e "${BLUE}Next step:${NC}"
    echo "  ./run_slc_pipeline.sh"
    echo
    exit 0
else
    echo -e "${RED}✗✗✗ SOME CHECKS FAILED - FIX ISSUES BEFORE RUNNING ✗✗✗${NC}"
    echo
    echo -e "${YELLOW}Common fixes:${NC}"
    echo "  1. Create environment: mamba env create -f env/scanpy.yml"
    echo "  2. Activate environment: conda activate scannex"
    echo "  3. Install Nextflow: curl -s https://get.nextflow.io | bash"
    echo
    exit 1
fi
