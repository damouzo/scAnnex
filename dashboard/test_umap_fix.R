#!/usr/bin/env Rscript
# Test script to verify UMAP loading works with the fix

# Setup Python path
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix != "") {
  python_path <- file.path(conda_prefix, "bin", "python3")
  message(sprintf("Using Conda Python: %s", python_path))
} else {
  python_path <- "/usr/bin/python3"
}

Sys.setenv(RETICULATE_PYTHON = python_path)

# Load libraries
library(reticulate)
use_python(python_path, required = TRUE)

# Import Python modules
ad <- import("anndata")

# Test file path
h5ad_path <- "/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad"

message("\n========================================")
message("Testing UMAP loading fix")
message("========================================\n")

# Test 1: Load with backed mode disabled (should work)
message("Test 1: Loading WITHOUT backed mode...")
tryCatch({
  adata <- ad$read_h5ad(h5ad_path, backed = NULL)
  obsm_keys <- py_to_r(adata$obsm$keys())
  message(sprintf("  ✓ .obsm keys available: %s", paste(obsm_keys, collapse=", ")))
  
  if ("X_umap" %in% obsm_keys) {
    umap_matrix <- py_to_r(adata$obsm["X_umap"])
    message(sprintf("  ✓ UMAP matrix shape: %d cells × %d dimensions", nrow(umap_matrix), ncol(umap_matrix)))
    message(sprintf("  ✓ UMAP range: [%.2f, %.2f] to [%.2f, %.2f]", 
                    min(umap_matrix[,1]), max(umap_matrix[,1]),
                    min(umap_matrix[,2]), max(umap_matrix[,2])))
  } else {
    message("  ✗ X_umap not found!")
  }
}, error = function(e) {
  message(sprintf("  ✗ ERROR: %s", e$message))
})

# Test 2: Load with backed mode enabled (problematic)
message("\nTest 2: Loading WITH backed mode...")
tryCatch({
  adata_backed <- ad$read_h5ad(h5ad_path, backed = "r")
  message("  ✓ File loaded in backed mode")
  
  # Try to access obsm
  tryCatch({
    obsm_keys <- py_to_r(adata_backed$obsm$keys())
    message(sprintf("  ✓ .obsm keys available: %s", paste(obsm_keys, collapse=", ")))
    
    if ("X_umap" %in% obsm_keys) {
      umap_matrix <- py_to_r(adata_backed$obsm["X_umap"])
      message(sprintf("  ✓ UMAP matrix accessible: %d cells × %d dimensions", 
                      nrow(umap_matrix), ncol(umap_matrix)))
    }
  }, error = function(e) {
    message(sprintf("  ✗ Cannot access .obsm in backed mode: %s", e$message))
    message("  → This is why we need to load separately or disable backed mode!")
  })
  
}, error = function(e) {
  message(sprintf("  ✗ ERROR: %s", e$message))
})

# Test 3: File size check (auto-disable backed mode logic)
message("\nTest 3: File size check for auto-disable backed mode...")
file_info <- file.info(h5ad_path)
file_size_mb <- file_info$size / (1024^2)
message(sprintf("  File size: %.1f MB", file_size_mb))

if (file_size_mb < 500) {
  message("  ✓ File is < 500 MB → auto-disabling backed mode is appropriate")
} else {
  message("  → File is large → backed mode recommended (with separate UMAP load)")
}

message("\n========================================")
message("Testing complete!")
message("========================================\n")
