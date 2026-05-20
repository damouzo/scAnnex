#!/usr/bin/env bash
# Run scAnnex with the synthetic PBMC demo dataset.
#
# Prerequisites:
#   1. Generate the H5AD files (once):
#        python data_demo/generate_demo.py
#
# Usage (from project root, interactive):
#        bash data_demo/run_command.sh
#
# For HPC/SLURM execution use submit.sh instead.

set -euo pipefail

# Navigate to the project root (parent of this script's directory).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."

# ── Data check ────────────────────────────────────────────────────────────────
if [[ ! -f data_demo/Ctrl_PBMC_C1_filtered_matrix.h5ad ]]; then
    echo "ERROR: Demo H5AD files not found. Generate them first:"
    echo "  python data_demo/generate_demo.py"
    exit 1
fi

# ── Run ───────────────────────────────────────────────────────────────────────
nextflow run . \
    --input              data_demo/samplesheet.csv \
    --contrasts_file     data_demo/contrasts.csv \
    --outdir             results_demo \
    --run_integration    true \
    --batch_key          sample_id \
    --dge_groupby        condition \
    --dge_reference      Control \
    --run_gsea           true \
    --organism           human \
    --celltypist_models  "Immune_All_Low.pkl" \
    -profile conda \
    --max_memory         8.GB \
    -resume
