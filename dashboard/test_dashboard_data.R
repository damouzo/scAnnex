#!/usr/bin/env Rscript

# Test script to diagnose dashboard data loading
cat("=== Dashboard Data Diagnostic ===\n\n")

# Set environment
Sys.setenv(SCANNEX_DATA_PATH = "/home/damo/scAnnex/results")

# Source global.R to load functions
setwd("/home/damo/scAnnex/dashboard")
source("global.R")

cat("\n1. Testing H5AD file loading...\n")
h5ad_path <- "/home/damo/scAnnex/results/auto/PBMC_1k_annotated.h5ad"
cat(sprintf("   Loading: %s\n", h5ad_path))

result <- load_h5ad_data(h5ad_path)

if (result$success) {
  cat("   ✓ Data loaded successfully\n")
  cat(sprintf("   - Cells: %d\n", nrow(result$metadata)))
  cat(sprintf("   - Genes: %d\n", result$n_genes))
  cat(sprintf("   - Metadata columns: %d\n", ncol(result$metadata)))
  cat(sprintf("   - UMAP coords: %d cells\n", nrow(result$umap_coords)))
  
  cat("\n2. Checking UMAP coordinates structure:\n")
  cat(sprintf("   - Column names: %s\n", paste(names(result$umap_coords), collapse=", ")))
  cat(sprintf("   - First few rows:\n"))
  print(head(result$umap_coords, 3))
  
  cat("\n3. Checking metadata structure:\n")
  cat(sprintf("   - Column names: %s\n", paste(names(result$metadata), collapse=", ")))
  cat(sprintf("   - First few rows:\n"))
  print(head(result$metadata[, 1:5], 3))
  
  cat("\n4. Testing merge of UMAP with metadata:\n")
  test_color_by <- "leiden_0.5"
  umap_merged <- merge(
    result$umap_coords,
    result$metadata[, c("cell_id", test_color_by), drop = FALSE],
    by = "cell_id"
  )
  cat(sprintf("   - Merged rows: %d\n", nrow(umap_merged)))
  cat(sprintf("   - Merged columns: %s\n", paste(names(umap_merged), collapse=", ")))
  cat(sprintf("   - First few rows:\n"))
  print(head(umap_merged, 3))
  
  cat("\n5. Testing gene expression retrieval:\n")
  gene_result <- get_gene_expression(h5ad_path, "CD3D")
  if (gene_result$success) {
    cat(sprintf("   ✓ CD3D expression loaded: %d cells\n", length(gene_result$expression)))
    cat(sprintf("   - Expression range: [%.2f, %.2f]\n", 
                min(gene_result$expression), max(gene_result$expression)))
    cat(sprintf("   - Non-zero cells: %d (%.1f%%)\n", 
                sum(gene_result$expression > 0),
                100 * sum(gene_result$expression > 0) / length(gene_result$expression)))
  } else {
    cat(sprintf("   ✗ Failed to load CD3D: %s\n", gene_result$message))
  }
  
} else {
  cat(sprintf("   ✗ Failed to load data: %s\n", result$message))
}

cat("\n=== Diagnostic Complete ===\n")
