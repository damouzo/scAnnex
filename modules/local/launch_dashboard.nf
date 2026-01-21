process LAUNCH_DASHBOARD {
    tag "Dashboard Ready"
    label 'process_low'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path(h5ad_file)

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
Results Directory: ${params.outdir}
Dashboard Port: ${port}
Dashboard Host: ${host}
Dashboard URL: http://${host}:${port}

Launch Command:
  cd ${projectDir}/dashboard
  bash launch_dashboard.sh ${params.outdir}

Or with R directly:
  cd ${projectDir}/dashboard
  Rscript -e "shiny::runApp(port=${port}, host='${host}')"

Documentation:
  ${projectDir}/dashboard/README.md
  ${projectDir}/dashboard/QUICKSTART.md
EOF

    # Print message to console (will appear in Nextflow output)
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  🎉 Pipeline Completed Successfully!"
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    echo "📊 Your interactive dashboard is ready to launch!"
    echo ""
    echo "📂 Results location:"
    echo "   ${params.outdir}"
    echo ""
    echo "🚀 To launch the dashboard, run:"
    echo ""
    echo "   cd ${projectDir}/dashboard"
    echo "   bash launch_dashboard.sh ${params.outdir}"
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  🌐 Dashboard URL (after launch):"
    echo "  http://${host}:${port}"
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    echo "💡 Tips:"
    echo "   - Click or copy the URL to open in your browser"
    echo "   - Use Ctrl+C to stop the dashboard"
    echo "   - Dashboard info saved to: ${params.outdir}/dashboard_info.txt"
    echo ""
    """

    stub:
    """
    touch dashboard_info.txt
    """
}
