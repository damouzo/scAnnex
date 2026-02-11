process LAUNCH_DASHBOARD {
    tag "Dashboard Auto-Launch"
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
    # Save dashboard info
    cat > dashboard_info.txt << EOF
Dashboard Configuration
======================
Results Directory: ${results_path}
Dashboard Port: ${port}
Dashboard Host: ${host}
Dashboard URL: http://${host}:${port}

Auto-launched: Yes
Launch time: \$(date)

Manual control commands:
  Stop dashboard:  lsof -ti:${port} | xargs kill
  Restart:         cd ${dashboard_dir} && bash launch_dashboard.sh ${results_path}
  View logs:       tail -f ${results_path}/dashboard_launch.log
EOF

    # Launch dashboard using background launcher
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo " Pipeline Completed - Launching Dashboard"
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    
    if bash ${dashboard_dir}/background_launch.sh ${results_path}; then
        echo ""
        echo "════════════════════════════════════════════════════════════════"
        echo "  Dashboard URL:"
        echo "  http://${host}:${port}"
        echo "════════════════════════════════════════════════════════════════"
        echo ""
        echo "Dashboard is running in the background"
        echo ""
        echo "Useful commands:"
        echo "  View logs:  tail -f ${results_path}/dashboard_launch.log"
        echo "  Stop:       lsof -ti:${port} | xargs kill"
        echo ""
    else
        echo ""
        echo "Failed to auto-launch dashboard"
        echo ""
        echo "To launch manually:"
        echo "  cd ${dashboard_dir}"
        echo "  bash launch_dashboard.sh ${results_path}"
        echo ""
        echo "════════════════════════════════════════════════════════════════"
        echo "  Dashboard URL (after manual launch):"
        echo "  http://${host}:${port}"
        echo "════════════════════════════════════════════════════════════════"
        echo ""
    fi
    """

    stub:
    """
    touch dashboard_info.txt
    """
}
