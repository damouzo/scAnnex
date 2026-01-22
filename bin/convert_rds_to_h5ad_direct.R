#!/usr/bin/env Rscript

#' Convert Seurat RDS to CSV format for Python import
#'
#' This script extracts the essential data from a Seurat object and exports
#' it as CSV files that can be easily read by Python/scanpy. This approach
#' avoids the need for SeuratDisk which is not available in conda.
#'
#' @param input_rds Path to input RDS file containing Seurat object
#' @param output_dir Path to output directory for CSV files
#'
#' Usage: Rscript convert_rds_to_h5ad_direct.R <input.rds> <output_dir>

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("ERROR: Incorrect number of arguments\n")
  cat("Usage: Rscript convert_rds_to_h5ad_direct.R <input.rds> <output_dir>\n")
  quit(status = 1)
}

input_rds <- args[1]
output_dir <- args[2]

# Validate input file exists
if (!file.exists(input_rds)) {
  cat(sprintf("ERROR: Input file does not exist: %s\n", input_rds))
  quit(status = 1)
}

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Print version information
cat("=== RDS to CSV Conversion ===\n")
cat(sprintf("Seurat version: %s\n", packageVersion("Seurat")))
cat(sprintf("Input: %s\n", input_rds))
cat(sprintf("Output directory: %s\n", output_dir))
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

# Get default assay
default_assay <- DefaultAssay(seurat_obj)
cat(sprintf("Using default assay: %s\n", default_assay))

# Extract count matrix (sparse format)
cat("Extracting count matrix...\n")

# Handle Seurat v5 vs v4 API differences
# v5 uses 'layer' parameter, v4 uses 'slot' parameter
seurat_version <- packageVersion("Seurat")
cat(sprintf("Seurat version: %s\n", seurat_version))

if (seurat_version >= "5.0.0") {
  # Seurat v5: use layer parameter
  cat("Using Seurat v5 API (layer parameter)...\n")
  counts <- tryCatch(
    {
      # Try to get counts layer
      GetAssayData(seurat_obj, layer = "counts", assay = default_assay)
    },
    error = function(e) {
      cat("Warning: 'counts' layer not found, trying 'data' layer...\n")
      GetAssayData(seurat_obj, layer = "data", assay = default_assay)
    }
  )
} else {
  # Seurat v4: use slot parameter
  cat("Using Seurat v4 API (slot parameter)...\n")
  counts <- GetAssayData(seurat_obj, slot = "counts", assay = default_assay)
}

# Convert sparse matrix to coordinate format and save
cat("Saving count matrix in Matrix Market format...\n")
writeMM(counts, file.path(output_dir, "matrix.mtx"))

# Save gene names
cat("Saving gene names...\n")
genes <- rownames(counts)
write.table(
  data.frame(gene = genes),
  file = file.path(output_dir, "genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Save cell barcodes
cat("Saving cell barcodes...\n")
barcodes <- colnames(counts)
write.table(
  data.frame(barcode = barcodes),
  file = file.path(output_dir, "barcodes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Save cell metadata
cat("Saving cell metadata...\n")
metadata <- seurat_obj@meta.data
write.csv(
  metadata,
  file = file.path(output_dir, "metadata.csv"),
  row.names = TRUE,
  quote = FALSE
)

# Save gene metadata (if available)
cat("Saving gene metadata...\n")
gene_meta <- seurat_obj[[default_assay]][[]]
if (nrow(gene_meta) > 0) {
  write.csv(
    gene_meta,
    file = file.path(output_dir, "gene_metadata.csv"),
    row.names = TRUE,
    quote = FALSE
  )
} else {
  cat("No gene metadata available\n")
}

# Print file sizes
cat("\n=== Output Files ===\n")
for (fname in c("matrix.mtx", "genes.tsv", "barcodes.tsv", "metadata.csv")) {
  fpath <- file.path(output_dir, fname)
  if (file.exists(fpath)) {
    fsize <- file.size(fpath) / 1024^2
    cat(sprintf("  %s: %.2f MB\n", fname, fsize))
  }
}

cat("\n✓ Successfully exported Seurat object to CSV format\n")
