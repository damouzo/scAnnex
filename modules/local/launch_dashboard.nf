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
    
    # Save dashboard info
    cat > dashboard_info.txt << EOF
Dashboard Configuration
======================
Results Directory: ${results_path}
Dashboard Port: ${port}
Dashboard Host: ${host}
Dashboard URL: http://${host}:${port}

Environment: \${ENV_NAME}
Setup completed: \$(date)

Launch command:
  cd ${dashboard_dir}
  bash launch_dashboard.sh ${results_path}

Alternative launch methods:
  bash ${dashboard_dir}/auto_launch_dashboard.sh ${results_path}

Control commands:
  Stop dashboard:  lsof -ti:${port} | xargs kill
  View logs:       tail -f dashboard_launch.log
EOF

    # Print completion message
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo " Interactive Dashboard Available"
    echo "════════════════════════════════════════════════════════════════"
     echo ""
    echo "  cd ${dashboard_dir}"
    echo "  bash launch_dashboard.sh ${results_path}"
    echo ""
    echo "Once started, access the dashboard at:"
    echo "  http://${host}:${port}"
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    """

    stub:
    """
    touch dashboard_info.txt
    """
}
