# Polyglot Annotation Suite - Implementation Guide

This document contains all remaining components for the modular multi-tool annotation suite.

---

## COMPONENT 1: Azimuth R Script

**File:** `bin/auto_annot_azimuth.R`

```r
#!/usr/bin/env Rscript

#==============================================================================
# Azimuth Automated Cell-Type Annotation
#==============================================================================
# Uses Azimuth reference-based annotation for Seurat objects
#
# Requirements:
#   - Azimuth >= 0.5.0
#   - Seurat >= 5.0.0
#   - SeuratObject
#
# Usage:
#   Rscript auto_annot_azimuth.R \\
#       --input seurat.rds \\
#       --output azimuth_annotations.csv \\
#       --reference pbmc
#==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(Azimuth)
})

# Command line arguments
option_list <- list(
    make_option(c("-i", "--input"), type = "character", required = TRUE,
                help = "Input RDS file (Seurat object)"),
    make_option(c("-o", "--output"), type = "character", required = TRUE,
                help = "Output CSV file (cell_id, label, score, tool)"),
    make_option(c("-r", "--reference"), type = "character", default = "pbmc",
                help = "Azimuth reference [pbmc, lung, kidney, etc.]"),
    make_option(c("--annotation-level"), type = "character", default = "predicted.celltype.l2",
                help = "Annotation level to extract"),
    make_option(c("--verbose"), action = "store_true", default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))

log_msg <- function(msg) {
    if (opt$verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

log_msg("Loading Seurat object...")
seurat_obj <- readRDS(opt$input)
log_msg(sprintf("Loaded: %d cells x %d genes", ncol(seurat_obj), nrow(seurat_obj)))

log_msg(sprintf("Running Azimuth with reference: %s", opt$reference))
seurat_obj <- RunAzimuth(seurat_obj, reference = opt$reference)

log_msg("Extracting annotations...")
cell_ids <- colnames(seurat_obj)
labels <- seurat_obj@meta.data[[opt$annotation_level]]
scores <- seurat_obj@meta.data[[paste0(opt$annotation_level, ".score")]]

# Handle missing scores
if (is.null(scores)) {
    scores <- rep(1.0, length(labels))
    log_msg("Warning: No scores available, using 1.0 for all cells")
}

# Create output DataFrame
results <- data.frame(
    cell_id = cell_ids,
    label = labels,
    score = scores,
    tool = "azimuth",
    stringsAsFactors = FALSE
)

log_msg(sprintf("Saving %d annotations to %s", nrow(results), opt$output))
write.csv(results, opt$output, row.names = FALSE, quote = FALSE)

log_msg("Azimuth annotation completed!")
cat(sprintf("\nUnique cell types: %d\n", length(unique(labels))))
cat(sprintf("Mean confidence: %.3f\n", mean(scores, na.rm = TRUE)))
```

---

## COMPONENT 2: Azimuth Nextflow Module

**File:** `modules/local/auto_annot_azimuth.nf`

```nextflow
process AUTO_ANNOT_AZIMUTH {
    tag "$meta.id"
    label 'process_high'

    container "satijalab/azimuth:0.5.0"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*_azimuth.csv"), emit: annotations
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    auto_annot_azimuth.R \\
        --input ${rds} \\
        --output ${prefix}_azimuth.csv \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
        Azimuth: \$(Rscript -e "cat(as.character(packageVersion('Azimuth')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "cell_id,label,score,tool" > ${prefix}_azimuth.csv
    echo "cell1,CD4 T cells,0.92,azimuth" >> ${prefix}_azimuth.csv
    touch versions.yml
    """
}
```

---

## COMPONENT 3: Annotation Merge Script

**File:** `bin/auto_annot_merge.py`

```python
#!/usr/bin/env python3

"""
Annotation Merge - Integrator for Multi-Tool Cell-Type Annotations
Combines annotations from multiple tools into a single H5AD object
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import scanpy as sc

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge multiple annotation tools into H5AD"
    )
    parser.add_argument("--input", required=True, help="Input H5AD file")
    parser.add_argument("--output", required=True, help="Output H5AD file")
    parser.add_argument("--annotations", nargs='+', required=True,
                       help="List of annotation CSV files")
    parser.add_argument("--backed", action="store_true",
                       help="Use backed mode for large datasets")
    return parser.parse_args()

def load_annotation_csv(csv_path):
    """Load annotation CSV and validate format."""
    df = pd.read_csv(csv_path)
    required_cols = ['cell_id', 'label', 'score', 'tool']
    
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"CSV must have columns: {required_cols}")
    
    tool_name = df['tool'].iloc[0]
    print(f"  Loaded {tool_name}: {len(df)} cells, {df['label'].nunique()} types")
    return df

def merge_annotations(adata, annotation_files):
    """Merge all annotation CSVs into adata.obs."""
    print(f"\nMerging {len(annotation_files)} annotation files...")
    
    for csv_path in annotation_files:
        print(f"\nProcessing: {csv_path}")
        df = load_annotation_csv(csv_path)
        
        tool_name = df['tool'].iloc[0]
        
        # Create column names
        label_col = f'celltype_{tool_name}'
        score_col = f'celltype_{tool_name}_score'
        
        # Set index to cell_id for merging
        df_indexed = df.set_index('cell_id')
        
        # Add to adata.obs
        adata.obs[label_col] = df_indexed.loc[adata.obs_names, 'label']
        adata.obs[score_col] = df_indexed.loc[adata.obs_names, 'score']
        
        print(f"  Added columns: {label_col}, {score_col}")
    
    return adata

def main():
    args = parse_args()
    
    print("="*60)
    print("Annotation Merge - Multi-Tool Integrator")
    print("="*60)
    
    # Load H5AD
    print(f"\nLoading H5AD: {args.input}")
    if args.backed:
        adata = sc.read_h5ad(args.input, backed='r')
        adata = adata.to_memory()  # Load to memory for writing
    else:
        adata = sc.read_h5ad(args.input)
    
    print(f"  Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Merge annotations
    adata = merge_annotations(adata, args.annotations)
    
    # Save merged H5AD
    print(f"\nSaving merged annotations to: {args.output}")
    adata.write_h5ad(args.output)
    
    print(f"\n✓ Merged H5AD saved with {len(args.annotations)} annotation tools")
    print(f"  New columns: {[c for c in adata.obs.columns if c.startswith('celltype_')]}")
    print("="*60)

if __name__ == "__main__":
    main()
```

---

## COMPONENT 4: Annotation Merge Module

**File:** `modules/local/auto_annot_merge.nf`

```nextflow
process AUTO_ANNOT_MERGE {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/scanpy_anndata:latest"

    input:
    tuple val(meta), path(h5ad), path(csvs)

    output:
    tuple val(meta), path("*_annotated.h5ad"), emit: h5ad
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def csv_files = csvs.collect{ it }.join(' ')
    """
    auto_annot_merge.py \\
        --input ${h5ad} \\
        --output ${prefix}_annotated.h5ad \\
        --annotations ${csv_files} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.h5ad
    touch versions.yml
    """
}
```

---

## COMPONENT 5: Configuration Updates

**Add to `nextflow.config`:**

```groovy
    // Annotation suite
    annotation {
        enabled           = true
        tools             = ['celltypist']  // Options: 'celltypist', 'azimuth'
        celltypist_model  = 'Immune_All_Low.pkl'
        celltypist_majority_voting = true
        azimuth_reference = 'pbmc'
    }
```

**Add to `conf/modules.config`:**

```groovy
    withName: 'H5AD_TO_RDS' {
        ext.args = '--use-raw'
        publishDir = [
            path: { "${params.outdir}/annotation/conversion" },
            mode: params.publish_dir_mode,
            pattern: '*.rds',
            enabled: params.annotation.tools.any { it in ['azimuth', 'singler'] }
        ]
    }

    withName: 'AUTO_ANNOT_CELLTYPIST' {
        ext.args = [
            "--model ${params.annotation.celltypist_model}",
            params.annotation.celltypist_majority_voting ? "--majority-voting" : "",
            "--download-model"
        ].join(' ').trim()
        publishDir = [
            path: { "${params.outdir}/annotation/celltypist" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'AUTO_ANNOT_AZIMUTH' {
        ext.args = "--reference ${params.annotation.azimuth_reference}"
        publishDir = [
            path: { "${params.outdir}/annotation/azimuth" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'AUTO_ANNOT_MERGE' {
        publishDir = [
            path: { "${params.outdir}/annotation" },
            mode: params.publish_dir_mode,
            pattern: '*_annotated.h5ad'
        ]
    }
```

---

## COMPONENT 6: Workflow Updates

**Add to `workflows/scannex.nf` (after NORMALIZE_INTEGRATE):**

```groovy
include { H5AD_TO_RDS          } from '../modules/local/h5ad_to_rds'
include { AUTO_ANNOT_CELLTYPIST} from '../modules/local/auto_annot_celltypist'
include { AUTO_ANNOT_AZIMUTH   } from '../modules/local/auto_annot_azimuth'
include { AUTO_ANNOT_MERGE     } from '../modules/local/auto_annot_merge'

    //
    // STEP 6: Multi-tool cell-type annotation (optional)
    //
    def integration_output = NORMALIZE_INTEGRATE.out.h5ad
    def annotation_csvs = Channel.empty()
    
    if (params.annotation.enabled && params.annotation.tools.size() > 0) {
        
        // Check if R-based tools are requested
        def r_tools = params.annotation.tools.findAll { it in ['azimuth', 'singler'] }
        def rds_channel = Channel.empty()
        
        if (r_tools.size() > 0) {
            // Convert H5AD to RDS for R-based tools
            H5AD_TO_RDS(integration_output)
            rds_channel = H5AD_TO_RDS.out.rds
        }
        
        // Run CellTypist (Python)
        if ('celltypist' in params.annotation.tools) {
            AUTO_ANNOT_CELLTYPIST(integration_output)
            annotation_csvs = annotation_csvs.mix(AUTO_ANNOT_CELLTYPIST.out.annotations)
        }
        
        // Run Azimuth (R)
        if ('azimuth' in params.annotation.tools) {
            AUTO_ANNOT_AZIMUTH(rds_channel)
            annotation_csvs = annotation_csvs.mix(AUTO_ANNOT_AZIMUTH.out.annotations)
        }
        
        // Merge all annotations
        annotation_csvs
            .map { meta, csv -> csv }
            .collect()
            .set { all_csvs }
        
        integration_output
            .combine(all_csvs)
            .set { merge_input }
        
        AUTO_ANNOT_MERGE(merge_input)
        integration_output = AUTO_ANNOT_MERGE.out.h5ad
    }
```

---

## COMPONENT 7: Dashboard Updates

**Update `dashboard/global.R` (add function):**

```r
#' Detect annotation columns automatically
#' @return Vector of column names starting with 'celltype_'
detect_annotation_columns <- function(metadata) {
  anno_cols <- grep("^celltype_", colnames(metadata), value = TRUE)
  # Exclude score columns
  anno_cols <- anno_cols[!grepl("_score$", anno_cols)]
  return(anno_cols)
}
```

**Update `dashboard/server.R` (in Color by dropdown logic):**

```r
# Dynamically detect annotation columns
annotation_cols <- detect_annotation_columns(data$metadata)

# Update color_by choices to include annotations
color_choices <- c(
  "Sample ID" = "sample_id",
  "Batch" = "batch",
  "Condition" = "condition"
)

# Add detected annotation columns
if (length(annotation_cols) > 0) {
  anno_choices <- setNames(
    annotation_cols,
    sapply(annotation_cols, function(x) {
      # Format: "CellType (celltypist)" from "celltype_celltypist"
      tool_name <- gsub("^celltype_", "", x)
      paste0("CellType (", tool_name, ")")
    })
  )
  color_choices <- c(color_choices, anno_choices)
}

updateSelectInput(session, "color_by", choices = color_choices)
```

---

## TESTING COMMANDS

```bash
# Test CellTypist only
nextflow run main.nf -profile test,docker \\
  --annotation.enabled true \\
  --annotation.tools celltypist \\
  --annotation.celltypist_model Immune_All_Low.pkl

# Test both tools
nextflow run main.nf -profile test,docker \\
  --annotation.enabled true \\
  --annotation.tools celltypist,azimuth \\
  --annotation.azimuth_reference pbmc

# Disable annotation
nextflow run main.nf -profile test,docker \\
  --annotation.enabled false
```

---

## OUTPUT STRUCTURE

```
results/
└── annotation/
    ├── conversion/
    │   └── sample.rds               # Seurat object (if R tools used)
    ├── celltypist/
    │   ├── sample_celltypist.csv    # Standardized annotations
    │   └── versions.yml
    ├── azimuth/
    │   ├── sample_azimuth.csv       # Standardized annotations
    │   └── versions.yml
    └── sample_annotated.h5ad        # ⭐ FINAL OUTPUT with all annotations
```

The final H5AD will contain columns like:
- `celltype_celltypist`: Labels from CellTypist
- `celltype_celltypist_score`: Confidence scores
- `celltype_azimuth`: Labels from Azimuth
- `celltype_azimuth_score`: Confidence scores

Dashboard will automatically detect and display these in the "Color by" dropdown!

---

**Status:** Complete modular annotation suite architecture ready for implementation.
