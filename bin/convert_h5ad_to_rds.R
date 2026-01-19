#!/usr/bin/env Rscript

#==============================================================================
# H5AD to RDS Conversion Script
#==============================================================================
# Converts AnnData H5AD files to Seurat v5 RDS objects for R-based annotation
# tools (e.g., Azimuth, SingleR)
#
# Requirements:
#   - Seurat >= 5.0.0
#   - SeuratDisk >= 0.0.0.9020
#   - reticulate
#   - anndata (Python)
#
# Usage:
#   Rscript convert_h5ad_to_rds.R --input integrated.h5ad --output seurat.rds
#==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(SeuratDisk)
    library(reticulate)
})

#==============================================================================
# Command Line Arguments
#==============================================================================

option_list <- list(
    make_option(c("-i", "--input"), 
                type = "character", 
                default = NULL,
                help = "Input H5AD file path", 
                metavar = "FILE"),
    
    make_option(c("-o", "--output"), 
                type = "character", 
                default = "seurat_object.rds",
                help = "Output RDS file path [default: %default]", 
                metavar = "FILE"),
    
    make_option(c("--assay-name"), 
                type = "character", 
                default = "RNA",
                help = "Name for the Seurat assay [default: %default]", 
                metavar = "STRING"),
    
    make_option(c("--use-raw"), 
                action = "store_true", 
                default = FALSE,
                help = "Use .raw.X from H5AD (raw counts) instead of .X"),
    
    make_option(c("--min-cells"), 
                type = "integer", 
                default = 0,
                help = "Minimum cells for feature filtering [default: %default]"),
    
    make_option(c("--min-features"), 
                type = "integer", 
                default = 0,
                help = "Minimum features for cell filtering [default: %default]"),
    
    make_option(c("--verbose"), 
                action = "store_true", 
                default = FALSE,
                help = "Print detailed progress messages")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "\nConvert H5AD (AnnData) to Seurat RDS object")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input H5AD file is required (--input)", call. = FALSE)
}

if (!file.exists(opt$input)) {
    stop(sprintf("Input file does not exist: %s", opt$input), call. = FALSE)
}

#==============================================================================
# Helper Functions
#==============================================================================

log_message <- function(msg) {
    if (opt$verbose) {
        cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
    }
}

#==============================================================================
# Main Conversion Pipeline
#==============================================================================

log_message("Starting H5AD to RDS conversion")
log_message(sprintf("Input: %s", opt$input))
log_message(sprintf("Output: %s", opt$output))

# Step 1: Load H5AD using reticulate + anndata
log_message("Loading H5AD file using Python anndata...")

tryCatch({
    # Import anndata
    ad <- import("anndata")
    
    # Read H5AD
    adata <- ad$read_h5ad(opt$input)
    
    log_message(sprintf("Loaded AnnData: %d cells x %d genes", 
                       nrow(adata$obs), 
                       nrow(adata$var)))
    
}, error = function(e) {
    stop(sprintf("Failed to load H5AD file: %s", e$message), call. = FALSE)
})

# Step 2: Extract count matrix
log_message("Extracting count matrix...")

if (opt$`use-raw` && !is.null(adata$raw)) {
    log_message("Using .raw.X (raw counts)")
    counts <- t(as.matrix(adata$raw$X))
    var_names <- adata$raw$var_names$to_list()
} else {
    log_message("Using .X (processed matrix)")
    counts <- t(as.matrix(adata$X))
    var_names <- adata$var_names$to_list()
}

obs_names <- adata$obs_names$to_list()

# Set row/column names
rownames(counts) <- var_names
colnames(counts) <- obs_names

log_message(sprintf("Count matrix dimensions: %d genes x %d cells", 
                   nrow(counts), ncol(counts)))

# Step 3: Extract metadata
log_message("Extracting cell metadata...")

metadata <- as.data.frame(adata$obs)
rownames(metadata) <- obs_names

log_message(sprintf("Metadata columns: %s", 
                   paste(colnames(metadata), collapse = ", ")))

# Step 4: Extract dimensional reductions (if available)
log_message("Checking for dimensional reductions...")

reductions <- list()

if ("X_pca" %in% names(adata$obsm)) {
    log_message("Found PCA coordinates")
    pca_coords <- as.matrix(adata$obsm[["X_pca"]])
    rownames(pca_coords) <- obs_names
    colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))
    reductions$pca <- pca_coords
}

if ("X_pca_harmony" %in% names(adata$obsm)) {
    log_message("Found Harmony-corrected PCA coordinates")
    harmony_coords <- as.matrix(adata$obsm[["X_pca_harmony"]])
    rownames(harmony_coords) <- obs_names
    colnames(harmony_coords) <- paste0("harmony_", 1:ncol(harmony_coords))
    reductions$harmony <- harmony_coords
}

if ("X_umap" %in% names(adata$obsm)) {
    log_message("Found UMAP coordinates")
    umap_coords <- as.matrix(adata$obsm[["X_umap"]])
    rownames(umap_coords) <- obs_names
    colnames(umap_coords) <- paste0("UMAP_", 1:ncol(umap_coords))
    reductions$umap <- umap_coords
}

# Step 5: Create Seurat object
log_message("Creating Seurat object...")

tryCatch({
    # Create Seurat object with counts
    seurat_obj <- CreateSeuratObject(
        counts = counts,
        meta.data = metadata,
        assay = opt$`assay-name`,
        min.cells = opt$`min-cells`,
        min.features = opt$`min-features`
    )
    
    log_message(sprintf("Created Seurat object: %d cells x %d genes", 
                       ncol(seurat_obj), nrow(seurat_obj)))
    
    # Add dimensional reductions
    if (length(reductions) > 0) {
        for (reduction_name in names(reductions)) {
            log_message(sprintf("Adding %s reduction", reduction_name))
            
            # Filter coordinates to match cells in Seurat object
            reduction_coords <- reductions[[reduction_name]]
            reduction_coords <- reduction_coords[colnames(seurat_obj), , drop = FALSE]
            
            # Create DimReduc object
            reduction_obj <- CreateDimReducObject(
                embeddings = reduction_coords,
                key = paste0(toupper(substring(reduction_name, 1, 1)), 
                           substring(reduction_name, 2), "_"),
                assay = opt$`assay-name`
            )
            
            seurat_obj[[reduction_name]] <- reduction_obj
        }
    }
    
    # Step 6: Add normalized data if available
    if (opt$`use-raw` == FALSE) {
        log_message("Adding normalized data to 'data' slot")
        seurat_obj <- SetAssayData(
            seurat_obj,
            slot = "data",
            new.data = counts  # Already normalized in H5AD .X
        )
    }
    
}, error = function(e) {
    stop(sprintf("Failed to create Seurat object: %s", e$message), call. = FALSE)
})

# Step 7: Save RDS
log_message(sprintf("Saving Seurat object to %s", opt$output))

tryCatch({
    saveRDS(seurat_obj, file = opt$output)
    log_message("Conversion completed successfully!")
    
}, error = function(e) {
    stop(sprintf("Failed to save RDS file: %s", e$message), call. = FALSE)
})

# Step 8: Print summary
cat("\n")
cat("==============================================\n")
cat("Conversion Summary\n")
cat("==============================================\n")
cat(sprintf("Input H5AD:        %s\n", opt$input))
cat(sprintf("Output RDS:        %s\n", opt$output))
cat(sprintf("Cells:             %d\n", ncol(seurat_obj)))
cat(sprintf("Genes:             %d\n", nrow(seurat_obj)))
cat(sprintf("Metadata columns:  %d\n", ncol(seurat_obj@meta.data)))
cat(sprintf("Reductions:        %s\n", 
           ifelse(length(reductions) > 0, 
                  paste(names(reductions), collapse = ", "), 
                  "None")))
cat("==============================================\n")

# Exit successfully
quit(status = 0)
