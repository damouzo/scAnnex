/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNIFY_INPUT             } from '../modules/local/unify_input'
include { QUALITY_CONTROL         } from '../modules/local/quality_control'
include { DOUBLET_DETECTION       } from '../modules/local/doublet_detection'
include { NORMALIZE_INTEGRATE     } from '../modules/local/normalize_integrate'
include { DIMENSIONALITY_REDUCTION } from '../modules/local/dimensionality_reduction'
include { AUTO_ANNOTATION         } from '../modules/local/auto_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCANNEX {
    main:
    //
    // STEP 1: Unify input format
    //
    def input_ch = channel.of(
        [
            [id: 'sample'],
            file(params.input, checkIfExists: true)
        ]
    )
    
    UNIFY_INPUT (
        input_ch,
        params.input_type
    )
    
    //
    // STEP 2: Quality control
    //
    QUALITY_CONTROL (
        UNIFY_INPUT.out.h5ad
    )
    
    //
    // STEP 3: Doublet detection (optional)
    //
    def qc_output = QUALITY_CONTROL.out.h5ad
    
    if (params.run_doublet_detection) {
        DOUBLET_DETECTION (
            qc_output
        )
        qc_output = DOUBLET_DETECTION.out.h5ad
    }
    
    //
    // STEP 4: Normalization and integration
    //
    NORMALIZE_INTEGRATE (
        qc_output
    )
    
    //
    // STEP 5: Dimensionality reduction and clustering
    //
    DIMENSIONALITY_REDUCTION (
        NORMALIZE_INTEGRATE.out.h5ad
    )
    
    //
    // STEP 6: Auto annotation (optional)
    //
    def final_output = DIMENSIONALITY_REDUCTION.out.h5ad
    
    if (params.run_auto_annotation && params.marker_list) {
        AUTO_ANNOTATION (
            final_output,
            file(params.marker_list, checkIfExists: true)
        )
        final_output = AUTO_ANNOTATION.out.h5ad
    }
    
    //
    // Emit final output
    //
    emit:
    h5ad = final_output
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
