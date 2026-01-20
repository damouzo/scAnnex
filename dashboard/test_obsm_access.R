#!/usr/bin/env Rscript
# Test different approaches to access .obsm in R

# Setup Python path
conda_prefix <- Sys.getenv("CONDA_PREFIX")
python_path <- file.path(conda_prefix, "bin", "python3")
Sys.setenv(RETICULATE_PYTHON = python_path)

library(reticulate)
use_python(python_path, required = TRUE)

ad <- import("anndata")
h5ad_path <- "/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad"

message("\nTest: Different ways to check for X_umap\n")

# Load data without backed mode
adata <- ad$read_h5ad(h5ad_path)

# Approach 1: Try direct access
message("Approach 1: Direct adata$obsm check...")
tryCatch({
  has_umap <- py_to_r(py_eval("'X_umap' in adata.obsm", convert = TRUE))
  message(sprintf("  Has X_umap: %s", has_umap))
}, error = function(e) {
  message(sprintf("  ERROR: %s", e$message))
})

# Approach 2: Use Python to get keys as list
message("\nApproach 2: Convert keys to list first...")
tryCatch({
  obsm_keys <- py_to_r(py_eval("list(adata.obsm.keys())", convert = TRUE))
  message(sprintf("  Keys: %s", paste(obsm_keys, collapse=", ")))
  message(sprintf("  'X_umap' in keys: %s", "X_umap" %in% obsm_keys))
}, error = function(e) {
  message(sprintf("  ERROR: %s", e$message))
})

# Approach 3: Try using py_has_attr
message("\nApproach 3: Get UMAP directly...")
tryCatch({
  umap_data <- py_to_r(adata$obsm["X_umap"])
  message(sprintf("  ✓ UMAP shape: %d × %d", nrow(umap_data), ncol(umap_data)))
  message(sprintf("  ✓ Range: [%.2f, %.2f] to [%.2f, %.2f]",
                  min(umap_data[,1]), max(umap_data[,1]),
                  min(umap_data[,2]), max(umap_data[,2])))
}, error = function(e) {
  message(sprintf("  ERROR: %s", e$message))
})

message("\nConclusion: Use py_eval to check, or just try-catch direct access\n")
