# Test loading data with R/reticulate
library(reticulate)

# Set Python
python_path <- file.path(Sys.getenv("CONDA_PREFIX", "/home/damo/miniforge3/envs/scannex-dashboard"), "bin", "python3")
use_python(python_path, required = TRUE)

# Import Python modules
ad <- import("anndata")
np <- import("numpy")

# Load H5AD
h5ad_path <- "../results/standard/PBMC_1k_processed.h5ad"
cat("Loading H5AD:", h5ad_path, "\n")

adata <- ad$read_h5ad(h5ad_path)

cat("\n=== Basic Info ===\n")
cat("Shape:", adata$n_obs, "×", adata$n_vars, "\n")
cat("obs columns:", paste(names(adata$obs), collapse=", "), "\n")
cat("obsm keys:", paste(names(adata$obsm), collapse=", "), "\n")

cat("\n=== UMAP Check ===\n")
if ("X_umap" %in% names(adata$obsm)) {
  umap_data <- py_to_r(adata$obsm["X_umap"])
  cat("✓ X_umap shape:", dim(umap_data), "\n")
  cat("  First 3 rows:\n")
  print(head(umap_data, 3))
} else {
  cat("✗ X_umap NOT FOUND\n")
}

cat("\n=== Gene Names Check ===\n")
gene_names <- py_to_r(adata$var_names)
cat("Total genes:", length(gene_names), "\n")
cat("First 5 genes:", paste(gene_names[1:5], collapse=", "), "\n")
cat("Gene 'CD3D' present:", "CD3D" %in% gene_names, "\n")

cat("\n=== Cell IDs Check ===\n")
cell_ids <- py_to_r(adata$obs_names)
cat("Total cells:", length(cell_ids), "\n")
cat("First 3 cells:", paste(cell_ids[1:3], collapse=", "), "\n")

cat("\n✓ All checks passed!\n")
