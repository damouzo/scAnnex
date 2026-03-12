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
include { AUTO_ANNOT_H5AD_TO_RDS  } from '../modules/local/auto_annot_h5ad_to_rds'
include { AUTO_ANNOT_SCTYPE       } from '../modules/local/auto_annot_sctype'
include { AUTO_ANNOT_AZIMUTH      } from '../modules/local/auto_annot_azimuth'
include { AUTO_ANNOT_SINGLER      } from '../modules/local/auto_annot_singler'
include { AUTO_ANNOT_SUMMARIZE    } from '../modules/local/auto_annot_summarize'
include { INTEGRATE_HARMONYPY     } from '../modules/local/integrate_harmonypy'
include { MERGE_SAMPLES           } from '../modules/local/merge_samples'
include { DIFFERENTIAL_EXPRESSION } from '../modules/local/differential_expression'
include { LAUNCH_DASHBOARD        } from '../modules/local/launch_dashboard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCANNEX {
    main:
    //
    // STEP 1: Parse input (supports both single files AND CSV samplesheets)
    //
    def input_ch
    if (params.input.endsWith('.csv') || params.input.endsWith('.tsv')) {
        // CSV/TSV samplesheet with multiple samples
        input_ch = channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, sep: params.input.endsWith('.tsv') ? '\t' : ',')
            .map { row ->
                def meta = [
                    id: row.sample_id ?: file(row.file_path).simpleName,
                    file_type: row.file_type ?: params.input_type,
                    batch: row.batch ?: 'batch1',
                    condition: row.condition ?: 'default'
                ]
                [ meta, file(row.file_path, checkIfExists: true) ]
            }
    } else {
        // Single file input (h5ad, rds, or mtx directory)
        input_ch = channel
            .fromPath(params.input, checkIfExists: true)
            .map { file ->
                def meta = [
                    id: file.simpleName,
                    file_type: params.input_type,
                    batch: 'batch1',
                    condition: 'default'
                ]
                [ meta, file ]
            }
    }
    
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
    // STEP 5: Integration (Optional - global multi-sample Harmony)
    //
    def base_annotation_h5ad

    if (params.run_integration && params.batch_key) {
        def integration_inputs = STANDARD_PROCESSING.out.h5ad.map { _meta, h5ad -> h5ad }.collect()
        def split_script = file("${projectDir}/bin/split_integrated_by_sample.py", checkIfExists: true)
        INTEGRATE_HARMONYPY (
            integration_inputs,
            split_script
        )
        base_annotation_h5ad = INTEGRATE_HARMONYPY.out.h5ad
    } else {
        def h5ad_files = STANDARD_PROCESSING.out.h5ad.map { _meta, h5ad -> h5ad }.collect()
        MERGE_SAMPLES (
            h5ad_files
        )
        base_annotation_h5ad = MERGE_SAMPLES.out.h5ad
    }

    //
    // STEP 6: Auto-Annotation (Global object, multi-tool in parallel)
    //
    def annotated_h5ad = base_annotation_h5ad

    if (params.run_auto_annotation) {
        def empty_annotations = Channel.fromPath("${projectDir}/assets/empty_annotations.csv", checkIfExists: true)
        def empty_status_json = Channel.fromPath("${projectDir}/assets/empty_status.json", checkIfExists: true)

        def celltypist_annotations_ch
        def celltypist_status_ch
        def azimuth_annotations_ch
        def azimuth_status_ch
        def singler_annotations_ch
        def singler_status_ch
        def sctype_annotations_ch
        def sctype_status_ch
        def auto_annot_rds_ch

        def run_any_r_annotator = params.azimuth_enable || params.singler_enable || params.sctype_enable
        if (run_any_r_annotator) {
            AUTO_ANNOT_H5AD_TO_RDS (
                base_annotation_h5ad
            )
            auto_annot_rds_ch = AUTO_ANNOT_H5AD_TO_RDS.out.rds
        }

        if (params.celltypist_enable) {
            AUTO_ANNOT_CELLTYPIST (
                base_annotation_h5ad.map { h5ad -> h5ad }
            )
            celltypist_annotations_ch = AUTO_ANNOT_CELLTYPIST.out.annotations
            celltypist_status_ch = AUTO_ANNOT_CELLTYPIST.out.status_json
        } else {
            celltypist_annotations_ch = empty_annotations
            celltypist_status_ch = empty_status_json
        }

        if (params.azimuth_enable) {
            AUTO_ANNOT_AZIMUTH (
                auto_annot_rds_ch
            )
            azimuth_annotations_ch = AUTO_ANNOT_AZIMUTH.out.annotations
            azimuth_status_ch = AUTO_ANNOT_AZIMUTH.out.status_json
        } else {
            azimuth_annotations_ch = empty_annotations
            azimuth_status_ch = empty_status_json
        }

        if (params.singler_enable) {
            AUTO_ANNOT_SINGLER (
                auto_annot_rds_ch
            )
            singler_annotations_ch = AUTO_ANNOT_SINGLER.out.annotations
            singler_status_ch = AUTO_ANNOT_SINGLER.out.status_json
        } else {
            singler_annotations_ch = empty_annotations
            singler_status_ch = empty_status_json
        }

        if (params.sctype_enable) {
            AUTO_ANNOT_SCTYPE (
                auto_annot_rds_ch
            )
            sctype_annotations_ch = AUTO_ANNOT_SCTYPE.out.annotations
            sctype_status_ch = AUTO_ANNOT_SCTYPE.out.status_json
        } else {
            sctype_annotations_ch = empty_annotations
            sctype_status_ch = empty_status_json
        }

        AUTO_ANNOT_SUMMARIZE (
            base_annotation_h5ad.map { h5ad -> h5ad },
            celltypist_annotations_ch,
            celltypist_status_ch,
            sctype_annotations_ch,
            sctype_status_ch,
            azimuth_annotations_ch,
            azimuth_status_ch,
            singler_annotations_ch,
            singler_status_ch
        )

        annotated_h5ad = AUTO_ANNOT_SUMMARIZE.out.h5ad
    }

    //
    // STEP 7: Differential Expression Analysis
    //
    if (params.run_dge) {
        DIFFERENTIAL_EXPRESSION (
            annotated_h5ad,
            params.contrasts_file ? 
                channel.fromPath(params.contrasts_file, checkIfExists: true) :
                channel.value(file('NO_FILE'))
        )
    }

    def final_output = annotated_h5ad.map { h5ad -> [[:], h5ad] }
    
    //
    // STEP 8: Launch Interactive Dashboard (Optional - enabled by default)
    //
    if (params.enable_dashboard) {
        // Calculate absolute results path
        def results_path = params.outdir.startsWith('/') ? params.outdir : "${workflow.launchDir}/${params.outdir}"
        
        LAUNCH_DASHBOARD (
            final_output.map { _meta, h5ad -> h5ad },
            "${projectDir}/dashboard",
            results_path
        )
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
