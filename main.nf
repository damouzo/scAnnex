#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scAnnex: Single-Cell RNA-seq Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/your-org/scannex
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCANNEX } from './workflows/scannex'

//
// WORKFLOW: Run main scAnnex analysis pipeline
//
workflow {
    // Print help message if needed
    if (params.help) {
        log.info paramsHelp("nextflow run main.nf --input samplesheet.csv")
        System.exit(0)
    }

    // Validate input parameters
    validateParameters()
    
    // Run main workflow
    SCANNEX ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success && params.enable_dashboard) {
        // Calculate absolute results path
        def results_path = params.outdir.startsWith('/') ? params.outdir : "${workflow.launchDir}/${params.outdir}"
        def dashboard_dir = "${workflow.projectDir}/dashboard"
        def port = params.dashboard_port ?: 3838
        def host = params.dashboard_host ?: 'localhost'
        
        log.info ""
        log.info "════════════════════════════════════════════════════════════════"
        log.info " Pipeline Completed Successfully"
        log.info "════════════════════════════════════════════════════════════════"
        log.info "Results saved to: ${results_path}"
        log.info "Interactive Dashboard Available. To launch it, run:"
        log.info "     cd ${dashboard_dir}"
        log.info "     bash launch_dashboard.sh ${results_path}"
        log.info ""
        log.info "Once started, access the dashboard at:"
        log.info "     http://${host}:${port}"
        log.info "════════════════════════════════════════════════════════════════"
        log.info ""
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
