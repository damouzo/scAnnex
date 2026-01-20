#!/usr/bin/env Rscript
# Final test of the complete load_h5ad_data function

# Setup Python path
conda_prefix <- Sys.getenv("CONDA_PREFIX")
python_path <- file.path(conda_prefix, "bin", "python3")
Sys.setenv(RETICULATE_PYTHON = python_path)

library(reticulate)
use_python(python_path, required = TRUE)

ad <- import("anndata")

# Source the global.R to get our function
source("global.R", local = TRUE)

h5ad_path <- "/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad"

message("\n========================================")
message("FINAL TEST: load_h5ad_data() function")
message("========================================\n")

# Test 1: Load with backed mode = FALSE (default for small files)
message("Test 1: Load with backed=FALSE (default for small files)...")
data_obj <- load_h5ad_data(h5ad_path, backed = FALSE)

message("\nResults:")
message(sprintf("  Cells: %d", data_obj$n_cells))
message(sprintf("  Genes: %d", data_obj$n_genes))
message(sprintf("  Backed mode: %s", data_obj$backed))
message(sprintf("  UMAP available: %s", !is.null(data_obj$umap_coords)))

if (!is.null(data_obj$umap_coords)) {
  message(sprintf("  UMAP dimensions: %d cells × %d coords", 
                  nrow(data_obj$umap_coords), 
                  ncol(data_obj$umap_coords) - 1))  # -1 for cell_id column
  message(sprintf("  UMAP_1 range: [%.2f, %.2f]", 
                  min(data_obj$umap_coords$UMAP_1),
                  max(data_obj$umap_coords$UMAP_1)))
  message(sprintf("  UMAP_2 range: [%.2f, %.2f]",
                  min(data_obj$umap_coords$UMAP_2),
                  max(data_obj$umap_coords$UMAP_2)))
  message("  ✓ SUCCESS: UMAP data loaded correctly!")
} else {
  message("  ✗ FAILED: No UMAP data!")
}

message("\n========================================")
message("Test complete - Dashboard should work!")
message("========================================\n")
