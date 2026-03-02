#!/usr/bin/env Rscript
# Test script to validate DGE results can be loaded

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# Set DGE results directory
dge_dir <- "results_dge_SUCCESS/dge/dge_results"

cat("Testing DGE results loading\n")
cat("================================\n\n")

# Check directory exists
if (!dir.exists(dge_dir)) {
  stop(sprintf("Directory not found: %s", dge_dir))
}

cat(sprintf("✓ DGE directory found: %s\n\n", dge_dir))

# Find all *_results.csv files (exclude all_contrasts_*)
result_files <- list.files(
  dge_dir,
  pattern = "^[^all].*_results\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0) {
  stop("No DGE results files found (looking for *_results.csv)")
}

cat(sprintf("✓ Found %d contrast result files:\n", length(result_files)))

# Extract contrast names
contrast_names <- gsub("_results\\.csv$", "", basename(result_files))

for (i in seq_along(result_files)) {
  cat(sprintf("  %d. %s\n", i, contrast_names[i]))
}

cat("\n")

# Load first contrast as test
cat(sprintf("Loading contrast: %s\n", contrast_names[1]))
cat("--------------------------------\n")

dge_df <- read.csv(result_files[1], stringsAsFactors = FALSE)

cat(sprintf("✓ Loaded %d genes\n", nrow(dge_df)))
cat(sprintf("✓ Columns: %s\n\n", paste(names(dge_df), collapse = ", ")))

# Show summary statistics
pval_threshold <- 0.05
logfc_threshold <- 0.25

dge_df$significant <- with(dge_df, 
  abs(log2_fc) >= logfc_threshold & 
  pvalue_adj < pval_threshold
)

n_sig <- sum(dge_df$significant)
n_up <- sum(dge_df$significant & dge_df$log2_fc > 0)
n_down <- sum(dge_df$significant & dge_df$log2_fc < 0)

cat("Summary Statistics:\n")
cat("--------------------------------\n")
cat(sprintf("Total genes tested: %d\n", nrow(dge_df)))
cat(sprintf("Significant genes (p < %.3f, |log2FC| >= %.2f): %d\n", 
            pval_threshold, logfc_threshold, n_sig))
cat(sprintf("  Upregulated: %d\n", n_up))
cat(sprintf("  Downregulated: %d\n\n", n_down))

# Create volcano plot
cat("Creating test volcano plot...\n")

dge_df$direction <- ifelse(
  !dge_df$significant, "Not Significant",
  ifelse(dge_df$log2_fc > 0, "Upregulated", "Downregulated")
)

p <- ggplot(dge_df, aes(x = log2_fc, y = -log10(pvalue_adj))) +
  geom_point(aes(color = direction), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#d62728",
      "Downregulated" = "#1f77b4",
      "Not Significant" = "gray70"
    )
  ) +
  geom_hline(yintercept = -log10(pval_threshold), 
             linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "gray30") +
  labs(
    title = sprintf("Test Volcano Plot: %s", contrast_names[1]),
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Add labels for top 10 significant genes
top_genes <- dge_df %>%
  filter(significant) %>%
  arrange(pvalue_adj) %>%
  head(10)

if (nrow(top_genes) > 0) {
  p <- p + 
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3
    )
}

# Save plot
output_file <- "test_volcano_plot.png"
ggsave(output_file, plot = p, width = 10, height = 8, dpi = 150)

cat(sprintf("✓ Volcano plot saved: %s\n\n", output_file))

# Show top 10 significant genes
if (nrow(top_genes) > 0) {
  cat("Top 10 Significant Genes:\n")
  cat("--------------------------------\n")
  
  top_genes_display <- top_genes %>%
    select(gene, log2_fc, pvalue, pvalue_adj, mean_expr_group1, mean_expr_group2) %>%
    mutate(
      log2_fc = round(log2_fc, 3),
      pvalue = format(pvalue, scientific = TRUE, digits = 3),
      pvalue_adj = format(pvalue_adj, scientific = TRUE, digits = 3),
      mean_expr_group1 = round(mean_expr_group1, 3),
      mean_expr_group2 = round(mean_expr_group2, 3)
    )
  
  print(top_genes_display, row.names = FALSE)
} else {
  cat("No significant genes found\n")
}

cat("\n✓ DGE loading test completed successfully!\n")
