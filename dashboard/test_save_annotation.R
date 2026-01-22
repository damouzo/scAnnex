#!/usr/bin/env Rscript

# Test script to verify annotation saving works
cat("=== Testing Annotation Save Function ===\n\n")

# Set environment
Sys.setenv(SCANNEX_DATA_PATH = "/home/damo/scAnnex/results")

# Source global.R to load functions
setwd("/home/damo/scAnnex/dashboard")
source("global.R")

cat("1. Loading test H5AD file...\n")
h5ad_path <- "/home/damo/scAnnex/results/auto/PBMC_1k_annotated.h5ad"
data_obj <- load_h5ad_data(h5ad_path)

cat("\n2. Creating test annotation...\n")
# Create simple test annotation
test_labels <- rep("Group_A", nrow(data_obj$metadata))
names(test_labels) <- data_obj$metadata$cell_id
test_labels[1:100] <- "Group_B"  # Mark first 100 cells as Group_B

cat(sprintf("   - Total cells: %d\n", length(test_labels)))
cat(sprintf("   - Group_A: %d\n", sum(test_labels == "Group_A")))
cat(sprintf("   - Group_B: %d\n", sum(test_labels == "Group_B")))

cat("\n3. Testing save_annotation_to_h5ad function...\n")
tryCatch({
  output_path <- save_annotation_to_h5ad(
    h5ad_path = h5ad_path,
    annotation_name = "test_annotation",
    labels = test_labels,
    output_path = "/tmp/test_annotated.h5ad",
    create_copy = TRUE
  )
  
  cat("   ✓ Save successful!\n")
  cat(sprintf("   - Output file: %s\n", output_path))
  cat(sprintf("   - File exists: %s\n", file.exists(output_path)))
  
  cat("\n4. Verifying saved annotation...\n")
  adata_verify <- ad$read_h5ad(output_path)
  cat(sprintf("   - 'test_annotation' in obs: %s\n", 
              "test_annotation" %in% names(adata_verify$obs)))
  
  if ("test_annotation" %in% names(adata_verify$obs)) {
    saved_labels <- py_to_r(adata_verify$obs$test_annotation)
    cat(sprintf("   - Saved labels length: %d\n", length(saved_labels)))
    cat(sprintf("   - Unique labels: %s\n", 
                paste(unique(as.character(saved_labels)), collapse=", ")))
  }
  
  cat("\n✓ All tests passed!\n")
  
}, error = function(e) {
  cat(sprintf("   ✗ Error: %s\n", e$message))
  cat("\nFull error:\n")
  print(e)
})

cat("\n=== Test Complete ===\n")
