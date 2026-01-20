/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS (SLC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNIFY_INPUT             } from '../modules/local/unify_input'
include { QUALITY_CONTROL         } from '../modules/local/quality_control'
include { DOUBLET_DETECTION       } from '../modules/local/doublet_detection'
include { STANDARD_PROCESSING     } from '../modules/local/standard_processing'
include { AUTO_ANNOT_CELLTYPIST   } from '../modules/local/auto_annot_celltypist'
include { NORMALIZE_INTEGRATE     } from '../modules/local/normalize_integrate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCANNEX {
    main:
    //
    // STEP 1: Parse samplesheet and create input channel
    //
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                file_type: row.file_type,
                batch: row.batch,
                condition: row.condition
            ]
            [ meta, file(row.file_path, checkIfExists: true) ]
        }
        .set { input_ch }
    
    //
    // STEP 2: Unify input format (per sample)
    //
    UNIFY_INPUT (
        input_ch
    )
    
    //
    // STEP 3: Quality Control (SLC with Quantile filtering & Attrition Log)
    //
    QUALITY_CONTROL (
        UNIFY_INPUT.out.h5ad
    )
    
    //
    // STEP 3: Doublet Detection (SLC with optional removal)
    //
    def processing_input = QUALITY_CONTROL.out.h5ad
    
    if (params.run_doublet_detection) {
        DOUBLET_DETECTION (
            processing_input
        )
        processing_input = DOUBLET_DETECTION.out.h5ad
    }
    
    //
    // STEP 4: Standard Processing (SLC: Normalize → PCA → UMAP → Multi-res Clustering)
    //
    STANDARD_PROCESSING (
        processing_input
    )
    
    //
    // STEP 5: Auto-Annotation with CellTypist (SLC)
    //
    def annotated_output = STANDARD_PROCESSING.out.h5ad
    
    if (params.run_auto_annotation) {
        AUTO_ANNOT_CELLTYPIST (
            annotated_output
        )
        annotated_output = AUTO_ANNOT_CELLTYPIST.out.h5ad
    }
    
    //
    // STEP 6: Integration (Optional - for multi-batch datasets)
    //
    def final_output = annotated_output
    
    if (params.run_integration && params.batch_key) {
        NORMALIZE_INTEGRATE (
            annotated_output
        )
        final_output = NORMALIZE_INTEGRATE.out.h5ad
    }
    
    //
    // Emit final outputs
    //
    emit:
    h5ad = final_output
    qc_results = QUALITY_CONTROL.out.qc_dir
    standard_results = STANDARD_PROCESSING.out.results_dir
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
