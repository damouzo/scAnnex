#!/usr/bin/env Rscript

#' Create demo RDS file from PBMC MTX data for scAnnex
#' This script creates a Seurat object from 10x MTX data

cat("Generating RDS demo file for scAnnex...\n\n")

# Check for required packages
required_packages <- c("Seurat")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("ERROR: Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Install with: install.packages('Seurat')\n")
  quit(status = 1)
}

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
})

cat("Package Version:\n")
cat(sprintf("  Seurat: %s\n\n", packageVersion("Seurat")))

# Set paths
script_dir <- dirname(normalizePath(commandArgs()[4]))
mtx_dir <- file.path(script_dir, "10xMTX", "filtered_feature_bc_matrix")
rds_dir <- file.path(script_dir, "RDS")

# Create output directory
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# Check MTX directory exists
if (!dir.exists(mtx_dir)) {
  cat(sprintf("ERROR: MTX directory not found: %s\n", mtx_dir))
  quit(status = 1)
}

# Load MTX data
cat(sprintf("Loading MTX data from: %s\n", mtx_dir))

tryCatch(
  {
    pbmc_data <- Read10X(data.dir = mtx_dir)
    cat(sprintf("✓ Loaded: %d genes × %d cells\n", nrow(pbmc_data), ncol(pbmc_data)))
  },
  error = function(e) {
    cat(sprintf("ERROR loading MTX: %s\n", e$message))
    quit(status = 1)
  }
)

# Create Seurat object
cat("\nCreating Seurat object...\n")

tryCatch(
  {
    seurat_obj <- CreateSeuratObject(
      counts = pbmc_data,
      project = "scAnnex_demo",
      min.cells = 3,
      min.features = 200
    )
    cat(sprintf("✓ Created Seurat object: %d cells × %d features\n",
                ncol(seurat_obj), nrow(seurat_obj)))
  },
  error = function(e) {
    cat(sprintf("ERROR creating Seurat object: %s\n", e$message))
    quit(status = 1)
  }
)

# Add metadata
cat("\nAdding metadata...\n")
seurat_obj$sample_id <- "PBMC_1k"
seurat_obj$batch <- "batch1"
seurat_obj$condition <- "control"
seurat_obj$orig.ident <- "PBMC_1k_demo"

cat("✓ Added metadata columns: sample_id, batch, condition\n")

# Calculate QC metrics
seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
cat(sprintf("✓ Calculated percent.mt (mean: %.2f%%)\n", mean(seurat_obj$percent.mt)))

# Basic normalization
cat("\nRunning basic normalization...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
cat("✓ Normalization complete\n")

# Save RDS
cat("\nSaving RDS file...\n")
rds_output <- file.path(rds_dir, "pbmc_1k.rds")

tryCatch(
  {
    saveRDS(seurat_obj, file = rds_output)
    file_size <- file.size(rds_output) / 1024^2
    cat(sprintf("✓ Created: %s\n", rds_output))
    cat(sprintf("✓ Size: %.2f MB\n", file_size))
    cat(sprintf("✓ Dimensions: %d cells × %d genes\n", 
                ncol(seurat_obj), nrow(seurat_obj)))
  },
  error = function(e) {
    cat(sprintf("ERROR saving RDS: %s\n", e$message))
    quit(status = 1)
  }
)

cat("\n✓ RDS demo file created successfully!\n")
