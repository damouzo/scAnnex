#!/usr/bin/env bash
#SBATCH --job-name=scannex-demo
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=data_demo/logs/scannex_demo_%j.out
#SBATCH --error=data_demo/logs/scannex_demo_%j.err
#
# SLURM submission script for the scAnnex demo dataset.
# Run from the project root:
#   sbatch data_demo/submit.sh
#
# Prerequisites (run once from an interactive node):
#   python data_demo/generate_demo.py   # also creates data_demo/logs/

set -euo pipefail

# SLURM copies the script to a temp location before executing it, so
# BASH_SOURCE[0] would resolve to /var/spool/slurmd/. Use SLURM_SUBMIT_DIR
# instead, which SLURM sets to the directory where 'sbatch' was called.
# This requires submitting from the project root (see usage comment above).
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

# ── Paths ─────────────────────────────────────────────────────────────────────
WORKDIR="/gpfs/scratch/${USER}/scannex_demo_work"

# ── Environment ───────────────────────────────────────────────────────────────
module unload openjdk 2>/dev/null || true
module load nextflow

# Run Nextflow from WORKDIR so .nextflow/ state files land on scratch (writable).
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# Samplesheet paths are relative to the project root; since Nextflow launches
# from WORKDIR, rewrite them to absolute paths in a temporary copy.
TMP_SAMPLESHEET="${WORKDIR}/samplesheet_abs.csv"
sed "s|data_demo/|${PROJECT_DIR}/data_demo/|g" \
    "${PROJECT_DIR}/data_demo/samplesheet.csv" > "${TMP_SAMPLESHEET}"

# ── Run ───────────────────────────────────────────────────────────────────────
nextflow run "${PROJECT_DIR}" \
    -profile singularity,apocrita \
    -w "${WORKDIR}" \
    --input              "${TMP_SAMPLESHEET}" \
    --contrasts_file     "${PROJECT_DIR}/data_demo/contrasts.csv" \
    --outdir             "${PROJECT_DIR}/results_demo" \
    --run_integration    true \
    --batch_key          sample_id \
    --dge_groupby        condition \
    --dge_reference      Control \
    --run_gsea           true \
    --organism           human \
    --celltypist_models  "Immune_All_Low.pkl" \
    -resume
