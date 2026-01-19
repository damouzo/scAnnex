/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Utility functions for scAnnex pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { fromSamplesheet        } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from 'plugin/nf-validation'

/*
========================================================================================
    SUBWORKFLOW TO INITIALIZE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Validate parameters
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:
    def ch_versions = channel.empty()

    //
    // Print version and exit if required
    //
    if (version) {
        log.info "${workflow.manifest.name} ${workflow.manifest.version}"
        System.exit(0)
    }

    //
    // Print help and exit if required
    //
    if (help) {
        log.info """
        ${workflow.manifest.name} v${workflow.manifest.version}
        
        ${workflow.manifest.description}
        
        Usage:
            nextflow run ${workflow.manifest.name} --input samplesheet.csv --outdir results -profile docker
        
        Mandatory arguments:
            --input                       Path to comma-separated file containing information about the samples
            --outdir                      The output directory where the results will be saved
        
        Optional arguments:
            --run_doublet_detection       Run Scrublet for doublet detection (default: true)
            --run_harmony                 Run Harmony for batch correction (default: true)
            --leiden_resolution           Leiden clustering resolution(s) (default: '0.5,1.0')
            --marker_list                 Path to CSV file with marker genes for auto-annotation
        
        Profiles:
            docker                        Use Docker containers
            singularity                   Use Singularity containers
            test                          Run with test dataset
        """.stripIndent()
        System.exit(0)
    }

    //
    // Validate parameters
    //
    if (validate_params) {
        validateParameters()
    }

    //
    // Create channel from input file provided through params.input
    //
    def ch_input = channel.fromSamplesheet('input')

    emit:
    samplesheet = ch_input
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    def summary_params = paramsSummaryMap(workflow, parameters_schema: 'nextflow_schema.json')

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            // TODO: Implement email notification
        }

        // Print summary
        log.info """
        ============================================
        Pipeline completed!
        ============================================
        Status:     ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Work dir:   ${workflow.workDir}
        Results:    ${outdir}
        Duration:   ${workflow.duration}
        ============================================
        """.stripIndent()
    }

    workflow.onError {
        log.error "Pipeline execution stopped with the following error: ${workflow.errorMessage}"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Check and validate pipeline parameters
//
def validateParameters() {
    if (!params.input) {
        error "Input samplesheet not specified with --input"
    }
}
