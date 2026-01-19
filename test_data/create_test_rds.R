#!/usr/bin/env Rscript

#' Create test RDS file from PBMC MTX data
#'
#' This script:
#' 1. Loads the PBMC 1k MTX dataset
#' 2. Creates a Seurat object
#' 3. Adds basic metadata
#' 4. Saves as RDS file
#'
#' Requirements: Seurat, SeuratDisk (for validation)

cat("====================================================================\n")
cat("Creating scAnnex Test RDS File\n")
cat("====================================================================\n\n")

# Check for required packages
required_packages <- c("Seurat", "SeuratDisk")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("ERROR: Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Install with: install.packages(c('Seurat', 'remotes'))\n")
  cat("              remotes::install_github('mojaveazure/seurat-disk')\n")
  quit(status = 1)
}

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

cat("Package Versions:\n")
cat(sprintf("  Seurat: %s\n", packageVersion("Seurat")))
cat(sprintf("  SeuratDisk: %s\n\n", packageVersion("SeuratDisk")))

# Set paths
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
mtx_dir <- file.path(script_dir, "mtx", "filtered_feature_bc_matrix")
rds_dir <- file.path(script_dir, "rds")

# Create output directory
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# Check MTX directory exists
if (!dir.exists(mtx_dir)) {
  cat(sprintf("ERROR: MTX directory not found: %s\n", mtx_dir))
  quit(status = 1)
}

# Load MTX data
cat(sprintf("1. Loading MTX data from: %s\n", mtx_dir))

tryCatch(
  {
    pbmc_data <- Read10X(data.dir = mtx_dir)
    cat(sprintf("   ✓ Loaded: %d genes × %d cells\n", nrow(pbmc_data), ncol(pbmc_data)))
  },
  error = function(e) {
    cat(sprintf("ERROR loading MTX: %s\n", e$message))
    quit(status = 1)
  }
)

# Create Seurat object
cat("\n2. Creating Seurat object...\n")

tryCatch(
  {
    seurat_obj <- CreateSeuratObject(
      counts = pbmc_data,
      project = "scAnnex_test",
      min.cells = 3,
      min.features = 200
    )
    cat(sprintf("   ✓ Created Seurat object: %d cells × %d features\n",
                ncol(seurat_obj), nrow(seurat_obj)))
  },
  error = function(e) {
    cat(sprintf("ERROR creating Seurat object: %s\n", e$message))
    quit(status = 1)
  }
)

# Add metadata
cat("\n3. Adding metadata...\n")
seurat_obj$sample_id <- "test_seurat"
seurat_obj$batch <- "batch2"
seurat_obj$condition <- "treated"
seurat_obj$orig.ident <- "PBMC_test"

cat("   ✓ Added metadata columns: sample_id, batch, condition, orig.ident\n")

# Add some basic QC metrics
seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
cat(sprintf("   ✓ Calculated percent.mt (mean: %.2f%%)\n", mean(seurat_obj$percent.mt)))

# Basic normalization (optional but realistic)
cat("\n4. Running basic normalization...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
cat("   ✓ Normalization complete\n")

# Save RDS
cat("\n5. Saving RDS file...\n")
rds_output <- file.path(rds_dir, "pbmc_seurat.rds")

tryCatch(
  {
    saveRDS(seurat_obj, file = rds_output)
    file_size <- file.size(rds_output) / 1024^2
    cat(sprintf("   ✓ Created: %s\n", rds_output))
    cat(sprintf("   ✓ Size: %.2f MB\n", file_size))
    cat(sprintf("   ✓ Dimensions: %d cells × %d genes\n", 
                ncol(seurat_obj), nrow(seurat_obj)))
  },
  error = function(e) {
    cat(sprintf("ERROR saving RDS: %s\n", e$message))
    quit(status = 1)
  }
)

# Verify we can read it back
cat("\n6. Verifying RDS can be reloaded...\n")
tryCatch(
  {
    seurat_test <- readRDS(rds_output)
    cat(sprintf("   ✓ Successfully reloaded: %d cells × %d genes\n",
                ncol(seurat_test), nrow(seurat_test)))
    cat(sprintf("   ✓ Metadata columns: %s\n", 
                paste(colnames(seurat_test@meta.data), collapse = ", ")))
    cat(sprintf("   ✓ Assays: %s\n", 
                paste(names(seurat_test@assays), collapse = ", ")))
  },
  error = function(e) {
    cat(sprintf("ERROR reloading RDS: %s\n", e$message))
    quit(status = 1)
  }
)

cat("\n====================================================================\n")
cat("✓ Test RDS file creation complete!\n")
cat("====================================================================\n\n")
cat("Created files:\n")
cat(sprintf("  - %s\n", rds_output))
cat("\nThis file can now be used to test the RDS → H5AD conversion pipeline.\n")
