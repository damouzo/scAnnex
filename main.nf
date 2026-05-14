#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scAnnex: Single-Cell RNA-seq Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/damouzo/scAnnex
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
    SCANNEX()

    workflow.onComplete = {
        if (workflow.success && params.enable_dashboard) {
            // Calculate absolute results path
            def results_path = params.outdir.startsWith('/') ? params.outdir : "${workflow.launchDir}/${params.outdir}"
            def dashboard_dir = "${workflow.projectDir}/dashboard"
            def port = params.dashboard_port ?: 3838
            def host = params.dashboard_host ?: 'localhost'

            // Detect HPC environment
            def isHPC = System.getenv('SLURM_CLUSTER_NAME') != null ||
                        'srun'.execute().text != null ||
                        workflow.executor == 'slurm'

            log.info ""
            log.info "════════════════════════════════════════════════════════════════"
            log.info " Pipeline Completed Successfully"
            log.info "════════════════════════════════════════════════════════════════"
            log.info "Results saved to: ${results_path}"

            if (isHPC) {
                // Estimate memory recommendation from the h5ad file size.
                // The dashboard loads the file fully into memory; in-memory
                // footprint is typically ~2x the on-disk size (HDF5 with sparse
                // count layers + dense log-norm + obsm + annotations).
                // Round up to the next standard size (8/16/32/64/128 G).
                def h5ad_file = new File("${results_path}/auto_annot/auto_annotated_global.h5ad")
                def rec_mem  = "16G"   // safe default
                def rec_cpus = 4
                if (h5ad_file.exists()) {
                    def size_gb   = h5ad_file.length() / 1e9
                    def raw_need  = size_gb * 2.0               // 2x expansion factor
                    def mem_tiers = [8, 16, 32, 64, 128]
                    def tier      = mem_tiers.find { it >= raw_need } ?: 128
                    rec_mem  = "${tier}G"
                    rec_cpus = tier >= 64 ? 16 : tier >= 32 ? 8 : 4
                }

                log.info "Interactive Dashboard Available"
                log.info ""
                log.info "To launch the dashboard, run:"
                log.info "     cd ${dashboard_dir}"
                log.info "     bash launch_dashboard_hpc.sh --mem ${rec_mem} --cpus ${rec_cpus} ${results_path}"
                log.info ""
                log.info "Memory recommendation: ${rec_mem} (estimated from output h5ad size)."
                log.info "This will request a SLURM interactive job with dedicated"
                log.info "compute resources and provide SSH tunnel instructions."
            } else {
                log.info "Interactive Dashboard Available. To launch it, run:"
                log.info "     cd ${dashboard_dir}"
                log.info "     bash launch_dashboard.sh ${results_path}"
                log.info ""
                log.info "Once started, access the dashboard at:"
                log.info "     http://${host}:${port}"
            }

            log.info "════════════════════════════════════════════════════════════════"
            log.info ""
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/