#!/usr/bin/env Rscript

#' Convert Seurat RDS to H5Seurat format
#'
#' This script loads a Seurat object from an RDS file and converts it to
#' the h5seurat format using SeuratDisk. This intermediate format preserves
#' all Seurat metadata and can be reliably converted to AnnData/H5AD in Python.
#'
#' @param input_rds Path to input RDS file containing Seurat object
#' @param output_h5seurat Path to output h5seurat file
#'
#' Usage: Rscript convert_rds_to_h5seurat.R <input.rds> <output.h5seurat>

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("ERROR: Incorrect number of arguments\n")
  cat("Usage: Rscript convert_rds_to_h5seurat.R <input.rds> <output.h5seurat>\n")
  quit(status = 1)
}

input_rds <- args[1]
output_h5seurat <- args[2]

# Validate input file exists
if (!file.exists(input_rds)) {
  cat(sprintf("ERROR: Input file does not exist: %s\n", input_rds))
  quit(status = 1)
}

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

# Print version information
cat("=== RDS to H5Seurat Conversion ===\n")
cat(sprintf("Seurat version: %s\n", packageVersion("Seurat")))
cat(sprintf("SeuratDisk version: %s\n", packageVersion("SeuratDisk")))
cat(sprintf("Input: %s\n", input_rds))
cat(sprintf("Output: %s\n", output_h5seurat))
cat("\n")

# Load Seurat object
cat("Loading Seurat object from RDS...\n")
tryCatch(
  {
    seurat_obj <- readRDS(input_rds)
  },
  error = function(e) {
    cat(sprintf("ERROR: Failed to load RDS file: %s\n", e$message))
    quit(status = 1)
  }
)

# Validate it's a Seurat object
if (!inherits(seurat_obj, "Seurat")) {
  cat(sprintf("ERROR: Object in RDS is not a Seurat object (class: %s)\n", 
              class(seurat_obj)[1]))
  quit(status = 1)
}

# Print basic information
cat(sprintf("Loaded Seurat object: %d cells x %d features\n",
            ncol(seurat_obj), nrow(seurat_obj)))

# Check for Seurat v5 assays and warn if present
if (packageVersion("Seurat") >= "5.0.0") {
  cat("Detected Seurat v5 - checking assay structure...\n")
  
  # Convert v5 assays to v3 for compatibility
  if (any(sapply(seurat_obj@assays, function(x) inherits(x, "Assay5")))) {
    cat("Converting Seurat v5 assays to v3 format for compatibility...\n")
    for (assay_name in names(seurat_obj@assays)) {
      if (inherits(seurat_obj@assays[[assay_name]], "Assay5")) {
        seurat_obj[[assay_name]] <- as(seurat_obj[[assay_name]], "Assay")
      }
    }
  }
}

# List available assays
cat(sprintf("Available assays: %s\n", paste(names(seurat_obj@assays), collapse = ", ")))

# List available reductions
if (length(names(seurat_obj@reductions)) > 0) {
  cat(sprintf("Available reductions: %s\n", 
              paste(names(seurat_obj@reductions), collapse = ", ")))
} else {
  cat("No dimensionality reductions found\n")
}

# Convert to H5Seurat format
cat("\nConverting to H5Seurat format...\n")
tryCatch(
  {
    SaveH5Seurat(
      seurat_obj, 
      filename = output_h5seurat,
      overwrite = TRUE,
      verbose = TRUE
    )
  },
  error = function(e) {
    cat(sprintf("ERROR: Failed to save H5Seurat file: %s\n", e$message))
    quit(status = 1)
  }
)

# Validate output file was created
if (!file.exists(output_h5seurat)) {
  cat("ERROR: H5Seurat file was not created successfully\n")
  quit(status = 1)
}

cat(sprintf("\n✓ Successfully converted RDS to H5Seurat: %s\n", output_h5seurat))
cat(sprintf("  File size: %.2f MB\n", file.size(output_h5seurat) / 1024^2))
