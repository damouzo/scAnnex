#!/bin/bash
#SBATCH --job-name=scannex-dashboard
#SBATCH --partition=interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=dashboard_%j.log
#SBATCH --error=dashboard_%j.err

#==============================================================================
# SLURM Job Script for scAnnex Dashboard on HPC
#==============================================================================
# This script launches the scAnnex interactive dashboard using Apptainer
# on an HPC cluster with SLURM scheduler.
#
# USAGE:
#   1. Edit the USER_DATA_DIR variable below to point to your results
#   2. Submit job: sbatch launch_dashboard.slurm
#   3. Wait for job to start (check with: squeue -u $USER)
#   4. Find the compute node and port in the log file
#   5. Set up SSH tunnel to access from your local machine
#
# REQUIREMENTS:
#   - Apptainer/Singularity installed on compute nodes
#   - scannex-dashboard.sif container built
#   - Results directory with .h5ad file
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# USER CONFIGURATION
#------------------------------------------------------------------------------

# Path to your results directory (will be mounted at /data in container)
USER_DATA_DIR="${SLURM_SUBMIT_DIR}/../results_slc_first_run"

# Path to container image
CONTAINER_IMAGE="${SLURM_SUBMIT_DIR}/scannex-dashboard.sif"

# Port for Shiny server (pick a unique port between 8000-9999)
SHINY_PORT=8888

#------------------------------------------------------------------------------
# SETUP
#------------------------------------------------------------------------------

echo "================================================================"
echo "  scAnnex Dashboard - SLURM Job"
echo "================================================================"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Node:         ${HOSTNAME}"
echo "Start time:   $(date)"
echo "Data dir:     ${USER_DATA_DIR}"
echo "Container:    ${CONTAINER_IMAGE}"
echo "Port:         ${SHINY_PORT}"
echo "================================================================"
echo ""

# Check if container exists
if [[ ! -f "$CONTAINER_IMAGE" ]]; then
    echo "ERROR: Container not found: $CONTAINER_IMAGE"
    echo "Build it first with: ./build_dashboard_container.sh apptainer"
    exit 1
fi

# Check if data directory exists
if [[ ! -d "$USER_DATA_DIR" ]]; then
    echo "ERROR: Data directory not found: $USER_DATA_DIR"
    exit 1
fi

# Find an available port (in case default is taken)
while ss -tln | grep -q ":${SHINY_PORT} "; do
    echo "Port ${SHINY_PORT} is in use, trying next port..."
    SHINY_PORT=$((SHINY_PORT + 1))
done

echo "Using port: ${SHINY_PORT}"
echo ""

#------------------------------------------------------------------------------
# SSH TUNNEL INSTRUCTIONS
#------------------------------------------------------------------------------

echo "================================================================"
echo "  ACCESS INSTRUCTIONS"
echo "================================================================"
echo ""
echo "To access the dashboard from your local machine:"
echo ""
echo "1. Open a new terminal on your LOCAL machine"
echo ""
echo "2. Set up SSH tunnel (keep this running):"
echo ""
echo "   ssh -N -L ${SHINY_PORT}:${HOSTNAME}:${SHINY_PORT} ${USER}@<your-hpc-login-node>"
echo ""
echo "3. Open your browser and navigate to:"
echo ""
echo "   http://localhost:${SHINY_PORT}"
echo ""
echo "================================================================"
echo ""
echo "Dashboard is starting..."
echo ""

#------------------------------------------------------------------------------
# LAUNCH DASHBOARD
#------------------------------------------------------------------------------

# Set environment variables
export SCANNEX_DATA_PATH=/data
export SHINY_PORT=${SHINY_PORT}

# Launch the dashboard with Apptainer
apptainer exec \
    --bind "${USER_DATA_DIR}:/data" \
    --env SCANNEX_DATA_PATH=/data \
    "${CONTAINER_IMAGE}" \
    R -e "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=${SHINY_PORT})"

echo ""
echo "================================================================"
echo "Dashboard stopped at: $(date)"
echo "================================================================"
