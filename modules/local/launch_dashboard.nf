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
    
    # Check if conda environment exists, create if needed
    ENV_NAME="scannex-dashboard"
    ENV_FILE="${dashboard_dir}/environment_dashboard.yml"
    
    if ! conda env list | grep -q "^\${ENV_NAME} "; then
        echo ""
        echo "════════════════════════════════════════════════════════════════"
        echo " First-time setup: Creating dashboard environment"
        echo "════════════════════════════════════════════════════════════════"
        echo ""
        echo "This will take approximately 5-10 minutes (one-time only)"
        echo ""
        
        conda env create -f "\${ENV_FILE}" -n "\${ENV_NAME}" || {
            echo ""
            echo "Warning: Failed to create dashboard environment automatically"
            echo ""
            echo "Please create manually:"
            echo "  cd ${dashboard_dir}"
            echo "  conda env create -f environment_dashboard.yml"
            echo ""
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
