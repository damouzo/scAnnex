process LAUNCH_DASHBOARD {
    tag "Dashboard Ready"
    label 'process_low'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path(h5ad_file)
    val(dashboard_dir)
    val(results_path)

    output:
    path("dashboard_info.txt"), emit: info

    when:
    params.enable_dashboard

    script:
    def port = params.dashboard_port ?: 3838
    def host = params.dashboard_host ?: 'localhost'
    """
    # Detect execution environment
    if [[ -n "\${SLURM_CLUSTER_NAME:-}" ]] || command -v srun &> /dev/null; then
        ENVIRONMENT="HPC"
    else
        ENVIRONMENT="LOCAL"
    fi
    
    # Determine conda environment location
    # Priority: 1) User-specified param, 2) HPC scratch, 3) Home directory
    if [[ -n "${params.dashboard_conda_dir}" ]]; then
        CONDA_BASE_DIR="${params.dashboard_conda_dir}"
    elif [[ -d "/gpfs/scratch/\${USER}" ]]; then
        # HPC with scratch (e.g., Apocrita)
        CONDA_BASE_DIR="/gpfs/scratch/\${USER}/conda_envs"
    else
        # Local system or HPC without scratch
        CONDA_BASE_DIR="\${HOME}/.conda/envs"
    fi
    
    ENV_NAME="scannex-dashboard"
    ENV_PATH="\${CONDA_BASE_DIR}/scannex-dashboard"
    ENV_FILE="${dashboard_dir}/environment_dashboard.yml"
    
    # Check if environment exists (by path, not name)
    if [[ -d "\${ENV_PATH}" ]]; then
        # Environment exists - check age and warn if old
        if command -v stat &> /dev/null; then
            ENV_AGE_SECONDS=\$(( \$(date +%s) - \$(stat -c %Y "\${ENV_PATH}" 2>/dev/null || stat -f %m "\${ENV_PATH}" 2>/dev/null || echo 0) ))
            ENV_AGE_DAYS=\$(( ENV_AGE_SECONDS / 86400 ))
            
            if [[ \${ENV_AGE_DAYS} -gt 60 ]] && [[ "\${CONDA_BASE_DIR}" == *"/gpfs/scratch/"* ]]; then
                echo ""
                echo "⚠️  Warning: Dashboard environment is \${ENV_AGE_DAYS} days old"
                echo "    (scratch files auto-delete after 65 days)"
                echo ""
                echo "    If dashboard fails, recreate environment:"
                echo "      rm -rf \${ENV_PATH}"
                echo "      bash launch_dashboard_hpc.sh ${results_path}"
                echo ""
            fi
        fi
    else
        # Environment does not exist - attempt auto-creation
        echo ""
        echo "════════════════════════════════════════════════════════════════"
        echo " First-time setup: Creating dashboard environment"
        echo "════════════════════════════════════════════════════════════════"
        echo ""
        echo "Location: \${ENV_PATH}"
        echo "This will take approximately 2-5 minutes (one-time only)"
        echo ""
        
        # Create parent directory if needed
        mkdir -p "\${CONDA_BASE_DIR}" || {
            echo "⚠️  Failed to create directory: \${CONDA_BASE_DIR}"
            echo "    Please check permissions or specify custom location:"
            echo "      --dashboard_conda_dir /path/to/writable/directory"
            echo ""
            exit 0  # Non-critical failure, skip dashboard
        }
        
        # Try to load miniforge module (HPC) or use system mamba
        if command -v module &> /dev/null; then
            module load miniforge 2>/dev/null || true
        fi
        
        # Prefer mamba over conda (10-20x faster dependency resolution)
        if command -v mamba &> /dev/null; then
            CONDA_CMD="mamba"
            echo "✓ Using mamba for faster environment creation"
        else
            CONDA_CMD="conda"
            echo "ℹ Using conda (mamba not available, will be slower)"
        fi
        echo ""
        
        # Attempt creation with 10-minute timeout
        # Use -p (path) instead of -n (name) to specify custom location
        timeout 600 \${CONDA_CMD} env create -f "\${ENV_FILE}" -p "\${ENV_PATH}" 2>&1 || {
            EXIT_CODE=\$?
            echo ""
            
            if [ \${EXIT_CODE} -eq 124 ]; then
                echo "⚠️  Timeout: Environment creation exceeded 10 minutes"
            elif [ \${EXIT_CODE} -eq 127 ]; then
                echo "⚠️  timeout command not available (old system)"
                echo "    Attempting without timeout..."
                \${CONDA_CMD} env create -f "\${ENV_FILE}" -p "\${ENV_PATH}" 2>&1 || {
                    echo "⚠️  Failed to create dashboard environment"
                }
            else
                echo "⚠️  Failed to create dashboard environment automatically"
            fi
            
            # Only show manual instructions if creation actually failed
            if [[ ! -d "\${ENV_PATH}" ]]; then
                echo ""
                echo "════════════════════════════════════════════════════════════════"
                echo " Manual Setup Instructions"
                echo "════════════════════════════════════════════════════════════════"
                echo ""
                echo "Option 1: HPC compute node (recommended for HPC)"
                echo "  srun --cpus-per-task=4 --mem=8G --time=30:00 --pty bash"
                echo "  module load miniforge  # if on HPC"
                echo "  mkdir -p \${CONDA_BASE_DIR}"
                echo "  cd ${dashboard_dir}"
                echo "  mamba env create -f environment_dashboard.yml -p \${ENV_PATH}"
                echo "  exit"
                echo ""
                echo "Option 2: Current node (if allowed)"
                echo "  cd ${dashboard_dir}"
                echo "  conda env create -f environment_dashboard.yml -p \${ENV_PATH}"
                echo ""
                echo "Then relaunch the dashboard:"
                echo "  bash launch_dashboard_hpc.sh ${results_path}"
                echo ""
                echo "════════════════════════════════════════════════════════════════"
                echo ""
            fi
        }
    fi
    
    # Save dashboard info (environment-aware)
    if [[ "\$ENVIRONMENT" == "HPC" ]]; then
        cat > dashboard_info.txt << 'EOF_HPC'
Dashboard Configuration (HPC Mode)
===================================
Results Directory: ${results_path}
Dashboard Port: ${port} (default, may auto-increment if busy)
Environment: scannex-dashboard
Setup completed: \$(date)

RECOMMENDED: HPC Production Mode
---------------------------------
For production use on HPC with dedicated compute resources:

  cd ${dashboard_dir}
  bash launch_dashboard_hpc.sh ${results_path}

This will:
  - Request SLURM interactive job on compute node
  - Allocate dedicated resources (configurable)
  - Provide SSH tunnel instructions
  - Follow HPC best practices

Resource recommendations:
  Dataset Size        Command
  -------------------------------------------------
  < 10k cells         bash launch_dashboard_hpc.sh
  10k-50k cells       bash launch_dashboard_hpc.sh --cpus 4 --mem 8G
  50k-100k cells      bash launch_dashboard_hpc.sh --cpus 8 --mem 16G
  > 100k cells        bash launch_dashboard_hpc.sh --cpus 16 --mem 32G

ALTERNATIVE: Quick Mode (current node)
--------------------------------------
For quick testing or small datasets:

  cd ${dashboard_dir}
  bash launch_dashboard.sh ${results_path}

Note: This runs on current node (login/compute).
Use only for small datasets or quick testing.

For detailed instructions:
  See docs/internal/HPC_DASHBOARD_GUIDE.md
EOF_HPC
    else
        cat > dashboard_info.txt << 'EOF_LOCAL'
Dashboard Configuration (Local Mode)
=====================================
Results Directory: ${results_path}
Dashboard Port: ${port} (auto-increment if busy)
Dashboard Host: ${host}
Environment: scannex-dashboard
Setup completed: \$(date)

Launch Dashboard:
-----------------
  cd ${dashboard_dir}
  bash launch_dashboard.sh ${results_path}

Once started, access the dashboard at:
  http://localhost:${port}

Control commands:
  Stop dashboard:  Press Ctrl+C in the terminal
  View logs:       cat ${dashboard_dir}/dashboard.log
EOF_LOCAL
    fi
    
    # Print completion message (environment-aware)
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo " scAnnex Pipeline Complete - Dashboard Ready"
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    
    if [[ "\$ENVIRONMENT" == "HPC" ]]; then
        echo "Environment: HPC Cluster Detected"
        echo ""
        echo "RECOMMENDED: Production Mode (dedicated compute resources)"
        echo ""
        echo "  cd ${dashboard_dir}"
        echo "  bash launch_dashboard_hpc.sh ${results_path}"
        echo ""
        echo "This will request a SLURM interactive job and provide"
        echo "SSH tunnel instructions for accessing the dashboard."
        echo ""
        echo "For custom resources (large datasets):"
        echo "  bash launch_dashboard_hpc.sh --cpus 8 --mem 16G ${results_path}"
        echo ""
        echo "────────────────────────────────────────────────────────────────"
        echo "ALTERNATIVE: Quick Mode (for testing/small datasets)"
        echo ""
        echo "  cd ${dashboard_dir}"
        echo "  bash launch_dashboard.sh ${results_path}"
        echo ""
        echo "Full guide: docs/internal/HPC_DASHBOARD_GUIDE.md"
    else
        echo "Environment: Local System"
        echo ""
        echo "Launch Dashboard:"
        echo ""
        echo "  cd ${dashboard_dir}"
        echo "  bash launch_dashboard.sh ${results_path}"
        echo ""
        echo "Once started, access at: http://localhost:${port}"
    fi
    
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    """

    stub:
    """
    touch dashboard_info.txt
    """
}
