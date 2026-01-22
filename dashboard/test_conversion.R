library(reticulate)

# Set Python
python_path <- "/home/damo/miniforge3/envs/scannex-dashboard/bin/python3"
Sys.setenv(RETICULATE_PYTHON = python_path)
use_python(python_path, required = TRUE)

# Import modules
ad <- import("anndata")
np <- import("numpy")

# Load data
h5ad_path <- "../results/standard/PBMC_1k_processed.h5ad"
cat("Loading:", h5ad_path, "\n")

adata <- ad$read_h5ad(h5ad_path)

# Test the exact conversion method from global.R
cat("\n=== Testing metadata conversion ===\n")
py$adata <- adata

py_run_string("
obs_dict = {}
for col in adata.obs.columns:
    obs_dict[col] = adata.obs[col].tolist()

obs_dict['cell_id'] = adata.obs_names.tolist()
")

obs_dict <- py$obs_dict
metadata <- as.data.frame(obs_dict, stringsAsFactors = FALSE, check.names = FALSE)

cat("✓ Metadata shape:", nrow(metadata), "×", ncol(metadata), "\n")
cat("  Columns:", paste(names(metadata), collapse=", "), "\n")
cat("  First cell_id:", metadata$cell_id[1], "\n")

# Test UMAP conversion
cat("\n=== Testing UMAP conversion ===\n")
umap_matrix <- py_to_r(adata$obsm["X_umap"])
cat("✓ UMAP shape:", dim(umap_matrix), "\n")

umap_coords <- data.frame(
  cell_id = metadata$cell_id,
  UMAP_1 = umap_matrix[, 1],
  UMAP_2 = umap_matrix[, 2]
)

cat("✓ UMAP coords shape:", nrow(umap_coords), "×", ncol(umap_coords), "\n")
cat("  First row:\n")
print(head(umap_coords, 1))

# Test gene names conversion
cat("\n=== Testing gene names ===\n")
py_run_string("
var_dict = {}
for col in adata.var.columns:
    var_dict[col] = adata.var[col].tolist()

var_dict['gene_id'] = adata.var_names.tolist()
")

var_dict <- py$var_dict
var_info <- as.data.frame(var_dict, stringsAsFactors = FALSE, check.names = FALSE)
rownames(var_info) <- var_info$gene_id

cat("✓ Var info shape:", nrow(var_info), "×", ncol(var_info), "\n")
cat("  Gene 'CD3D' present:", "CD3D" %in% rownames(var_info), "\n")

cat("\n✓ All conversions successful!\n")
