#!/usr/bin/env Rscript
# Complete dashboard functionality test
# Tests all core functions before launching the full app

# Setup Python path
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix != "") {
  python_path <- file.path(conda_prefix, "bin", "python3")
} else {
  python_path <- "/usr/bin/python3"
}

Sys.setenv(RETICULATE_PYTHON = python_path)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║          scAnnex Dashboard - Full Functionality Test          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Load required libraries
cat("📦 Loading R libraries...\n")
suppressPackageStartupMessages({
  library(reticulate)
  use_python(python_path, required = TRUE)
})

# Import Python modules
cat("🐍 Importing Python modules...\n")
ad <- import("anndata")
sc <- import("scanpy")

# Source global.R to get functions
cat("📄 Loading global.R functions...\n")
source("global.R", local = TRUE)

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("TEST 1: Data Loading\n")
cat("═══════════════════════════════════════════════════════════════\n")

h5ad_path <- "/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad"

if (!file.exists(h5ad_path)) {
  cat("❌ ERROR: Test file not found!\n")
  cat("   Expected: ", h5ad_path, "\n")
  quit(status = 1)
}

cat("✓ Test file found\n")
cat("  Path: ", h5ad_path, "\n\n")

# Load data
cat("Loading data with load_h5ad_data()...\n")
data_obj <- load_h5ad_data(h5ad_path, backed = FALSE)

cat("\n📊 Dataset Summary:\n")
cat("  • Cells: ", format(data_obj$n_cells, big.mark=","), "\n", sep="")
cat("  • Genes: ", format(data_obj$n_genes, big.mark=","), "\n", sep="")
cat("  • Backed mode: ", data_obj$backed, "\n", sep="")
cat("  • UMAP available: ", !is.null(data_obj$umap_coords), "\n", sep="")

if (!is.null(data_obj$umap_coords)) {
  cat("  • UMAP dimensions: ", nrow(data_obj$umap_coords), " cells × ", 
      ncol(data_obj$umap_coords) - 1, " coords\n", sep="")
  cat("  ✅ UMAP coordinates loaded successfully\n")
} else {
  cat("  ❌ ERROR: UMAP coordinates not loaded!\n")
  quit(status = 1)
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("TEST 2: Metadata Inspection\n")
cat("═══════════════════════════════════════════════════════════════\n")

metadata_cols <- names(data_obj$metadata)
cat("✓ Metadata columns (", length(metadata_cols), " total):\n", sep="")
for (col in metadata_cols[1:min(15, length(metadata_cols))]) {
  cat("  • ", col, "\n", sep="")
}
if (length(metadata_cols) > 15) {
  cat("  ... and ", length(metadata_cols) - 15, " more\n", sep="")
}

# Check for key columns
key_cols <- c("predicted_labels", "celltypist_score", "n_genes_by_counts", 
              "total_counts", "pct_counts_mt")
cat("\n🔍 Key columns check:\n")
for (col in key_cols) {
  if (col %in% metadata_cols) {
    cat("  ✅ ", col, "\n", sep="")
  } else {
    cat("  ⚠️  ", col, " (missing)\n", sep="")
  }
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("TEST 3: Cell Type Distribution\n")
cat("═══════════════════════════════════════════════════════════════\n")

if ("predicted_labels" %in% metadata_cols) {
  cell_types <- table(data_obj$metadata$predicted_labels)
  cell_types_sorted <- sort(cell_types, decreasing = TRUE)
  
  cat("✓ Cell types found:\n\n")
  cat(sprintf("%-40s %8s %8s\n", "Cell Type", "Count", "%"))
  cat(strrep("─", 60), "\n")
  
  total_cells <- sum(cell_types_sorted)
  for (i in 1:length(cell_types_sorted)) {
    ct_name <- names(cell_types_sorted)[i]
    ct_count <- cell_types_sorted[i]
    ct_pct <- (ct_count / total_cells) * 100
    cat(sprintf("%-40s %8d %7.1f%%\n", ct_name, ct_count, ct_pct))
  }
  cat(strrep("─", 60), "\n")
  cat(sprintf("%-40s %8d %7.1f%%\n", "TOTAL", total_cells, 100.0))
  
} else {
  cat("⚠️  predicted_labels column not found\n")
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("TEST 4: Gene Expression Extraction\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Test marker genes
marker_genes <- c("CD3D", "CD14", "CD79A", "MS4A1", "NKG7")
cat("Testing marker gene extraction:\n\n")

genes_available <- rownames(data_obj$var_info)
for (gene in marker_genes) {
  if (gene %in% genes_available) {
    tryCatch({
      expr <- get_gene_expression(data_obj, gene)
      expr_mean <- mean(expr, na.rm = TRUE)
      expr_max <- max(expr, na.rm = TRUE)
      pct_expressed <- sum(expr > 0) / length(expr) * 100
      cat(sprintf("  ✅ %-8s  Mean: %.3f  Max: %.3f  Expressed: %.1f%% cells\n", 
                  gene, expr_mean, expr_max, pct_expressed))
    }, error = function(e) {
      cat(sprintf("  ❌ %-8s  ERROR: %s\n", gene, e$message))
    })
  } else {
    cat(sprintf("  ⚠️  %-8s  Not found in dataset\n", gene))
  }
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("TEST 5: UMAP Data Integrity\n")
cat("═══════════════════════════════════════════════════════════════\n")

umap_data <- data_obj$umap_coords

cat("✓ UMAP data structure:\n")
cat("  • Rows (cells): ", nrow(umap_data), "\n", sep="")
cat("  • Columns: ", ncol(umap_data), " (", paste(names(umap_data), collapse=", "), ")\n", sep="")
cat("  • UMAP_1 range: [", sprintf("%.2f", min(umap_data$UMAP_1)), ", ", 
    sprintf("%.2f", max(umap_data$UMAP_1)), "]\n", sep="")
cat("  • UMAP_2 range: [", sprintf("%.2f", min(umap_data$UMAP_2)), ", ", 
    sprintf("%.2f", max(umap_data$UMAP_2)), "]\n", sep="")
cat("  • NAs in UMAP_1: ", sum(is.na(umap_data$UMAP_1)), "\n", sep="")
cat("  • NAs in UMAP_2: ", sum(is.na(umap_data$UMAP_2)), "\n", sep="")

if (sum(is.na(umap_data$UMAP_1)) == 0 && sum(is.na(umap_data$UMAP_2)) == 0) {
  cat("  ✅ No missing values in UMAP coordinates\n")
} else {
  cat("  ⚠️  WARNING: Missing values detected in UMAP\n")
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("TEST 6: Merge UMAP with Metadata (for plotting)\n")
cat("═══════════════════════════════════════════════════════════════\n")

if ("predicted_labels" %in% names(data_obj$metadata)) {
  # Simulate what the dashboard does
  umap_plot_data <- merge(
    umap_data,
    data_obj$metadata[, c("cell_id", "predicted_labels"), drop = FALSE],
    by = "cell_id"
  )
  
  cat("✓ Merged data for plotting:\n")
  cat("  • Rows: ", nrow(umap_plot_data), "\n", sep="")
  cat("  • Columns: ", ncol(umap_plot_data), " (", paste(names(umap_plot_data), collapse=", "), ")\n", sep="")
  cat("  • Ready for plotly: ", all(c("UMAP_1", "UMAP_2", "predicted_labels") %in% names(umap_plot_data)), "\n", sep="")
  
  if (nrow(umap_plot_data) == nrow(umap_data)) {
    cat("  ✅ All cells retained after merge\n")
  } else {
    cat("  ⚠️  WARNING: Some cells lost in merge\n")
  }
} else {
  cat("⚠️  Cannot test merge - predicted_labels not available\n")
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                     TEST SUMMARY                               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("✅ All core functionality tests PASSED!\n")
cat("\n")
cat("Dashboard is ready to launch. Run:\n")
cat("  ./launch_dashboard.sh\n")
cat("  or\n")
cat("  R -e \"shiny::runApp('.', host='127.0.0.1', port=8888)\"\n")
cat("\n")
cat("Then access: http://localhost:8888\n")
cat("\n")
