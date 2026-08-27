#!/usr/bin/env Rscript

# Pseudobulk differential expression with DESeq2, per cell type.
#
# The unit of biological replication is the sample x cell type aggregate
# (never individual cells), avoiding pseudo-replication. A cell type is only
# tested when BOTH conditions of the contrast have >= 2 biological replicates
# (samples); otherwise it is catalogued as SKIP_INSUFFICIENT_REPLICATES.

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
  library(jsonlite)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--counts"), dest = "counts", type = "character", help = "Genes x aggregates counts TSV (integer)"),
  make_option(c("--metadata"), dest = "metadata", type = "character", help = "Aggregate metadata CSV"),
  make_option(c("--outdir"), dest = "outdir", type = "character", default = "pseudobulk_dge", help = "Output directory"),
  make_option(c("--condition-col"), dest = "condition_col", type = "character", default = "condition", help = "Condition column"),
  make_option(c("--groupby"), dest = "groupby", type = "character", default = "cell_type", help = "Cell type column"),
  make_option(c("--sample-col"), dest = "sample_col", type = "character", default = "sample_id", help = "Sample column"),
  make_option(c("--control-group"), dest = "control_group", type = "character", help = "Reference/control condition"),
  make_option(c("--padj-cutoff"), dest = "padj_cutoff", type = "double", default = 0.05, help = "Adjusted p-value cutoff"),
  make_option(c("--shrink"), dest = "shrink", type = "character", default = "apeglm", help = "LFC shrinkage type (apeglm/ashr/normal)"),
  make_option(c("--min-replicates"), dest = "min_replicates", type = "integer", default = 2, help = "Min samples per contrast group")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

slugify <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

agg_id <- function(sample_id, cell_type) paste(sample_id, cell_type, sep = "::")

meta <- fread(opt$metadata, data.table = FALSE)
meta[[opt$groupby]] <- as.character(meta[[opt$groupby]])
meta[[opt$condition_col]] <- as.character(meta[[opt$condition_col]])
meta[[opt$sample_col]] <- as.character(meta[[opt$sample_col]])
meta$agg_id <- agg_id(meta[[opt$sample_col]], meta[[opt$groupby]])

if (is.null(opt$control_group) || !nzchar(opt$control_group)) {
  stop("--control-group (reference condition) is required")
}
control <- opt$control_group

counts <- fread(opt$counts, data.table = FALSE)
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
counts <- as.matrix(counts)
mode(counts) <- "integer"

disease <- setdiff(unique(meta[[opt$condition_col]]), control)
if (length(disease) == 0) {
  message("No non-control conditions present; aborting")
  status <- list(success = FALSE, message = "No non-control conditions present", cell_types = list())
  writeLines(toJSON(status, pretty = TRUE, auto_unbox = TRUE), file.path(opt$outdir, "pseudobulk_status.json"))
  file.create(file.path(opt$outdir, "pseudobulk_summary.csv"))
  quit(status = 0)
}

make_status <- function(success, msg, entries) {
  logs <- lapply(entries, function(e) {
    list(
      cell_type = e$cell_type,
      contrast = e$contrast,
      status = e$status,
      message = e$message,
      n_control_samples = e$n_control_samples,
      n_treatment_samples = e$n_treatment_samples
    )
  })
  list(success = success, message = msg, contrasts = logs)
}

shrink_res <- function(dds, cond) {
  coef_name <- paste0("condition_", cond, "_vs_", control)
  shrink_type <- tolower(opt$shrink)
  if (!shrink_type %in% c("apeglm", "ashr", "normal")) shrink_type <- "normal"
  if (shrink_type == "normal") {
    return(DESeq2::lfcShrink(dds, coef = coef_name, type = "normal"))
  }
  res <- tryCatch(
    DESeq2::lfcShrink(dds, coef = coef_name, type = shrink_type),
    error = function(e) {
      message(sprintf("Shrink '%s' failed for %s, falling back to normal: %s", shrink_type, cond, conditionMessage(e)))
      DESeq2::lfcShrink(dds, coef = coef_name, type = "normal")
    }
  )
  res
}

plot_volcano <- function(res_df, contrast_id, file_id, outdir) {
  df <- data.frame(
    log2fc = res_df$log2fc,
    padj = res_df$padj,
    stringsAsFactors = FALSE
  )
  df$neg_log10_padj <- -log10(df$padj + 1e-300)
  df$sig <- ifelse(df$padj < opt$padj_cutoff & df$log2fc > 1, "Up",
    ifelse(df$padj < opt$padj_cutoff & df$log2fc < -1, "Down", "Not Sig"))
  p <- ggplot(df, aes(x = log2fc, y = neg_log10_padj)) +
    geom_point(aes(color = sig), size = 0.6, alpha = 0.7) +
    scale_color_manual(values = c(Up = "red", Down = "blue", `Not Sig` = "grey60")) +
    geom_hline(yintercept = -log10(opt$padj_cutoff), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    labs(title = contrast_id, x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    theme_minimal(base_size = 12)
  ggsave(file.path(outdir, paste0(file_id, "_volcano.pdf")), plot = p, width = 8, height = 7)
  ggsave(file.path(outdir, paste0(file_id, "_volcano.png")), plot = p, width = 8, height = 7, dpi = 150)
}

entries <- list()
results_log <- list()
n_run <- 0L
skipped <- 0L

# Defensive guard: if two distinct contrast ids slugify to the same file id,
# the module (and downstream GSEA) would silently overwrite one another.
# Fail loudly instead of losing data.
seen_file_ids <- list()
for (ct in sort(unique(meta[[opt$groupby]]))) {
  idx_ct <- meta[[opt$groupby]] == ct
  meta_ct <- meta[idx_ct, , drop = FALSE]
  ct_present <- setdiff(unique(meta_ct[[opt$condition_col]]), control)
  for (cond in sort(ct_present)) {
    contrast_id <- paste0(cond, "_vs_", control, "__", ct)
    file_id <- slugify(contrast_id)
    if (file_id %in% names(seen_file_ids)) {
      stop(sprintf(
        "Pseudo-bulk contrast file-id collision: '%s' and '%s' both slugify to '%s'. ",
        seen_file_ids[[file_id]], contrast_id, file_id
      ))
    }
    seen_file_ids[[file_id]] <- contrast_id
  }
}

for (ct in sort(unique(meta[[opt$groupby]]))) {
  idx_ct <- meta[[opt$groupby]] == ct
  meta_ct <- meta[idx_ct, , drop = FALSE]
  ct_present <- setdiff(unique(meta_ct[[opt$condition_col]]), control)

  for (cond in sort(ct_present)) {
    contrast_id <- paste0(cond, "_vs_", control, "__", ct)
    file_id <- slugify(contrast_id)
    sub <- meta_ct[meta_ct[[opt$condition_col]] %in% c(control, cond), , drop = FALSE]

    n_control <- length(unique(sub[[opt$sample_col]][sub[[opt$condition_col]] == control]))
    n_treat <- length(unique(sub[[opt$sample_col]][sub[[opt$condition_col]] == cond]))

    ok <- n_control >= opt$min_replicates && n_treat >= opt$min_replicates
    if (!ok) {
      skipped <- skipped + 1L
      results_log[[length(results_log) + 1]] <- list(
        cell_type = ct, contrast = contrast_id, status = "SKIP_INSUFFICIENT_REPLICATES",
        message = "Fewer than 2 samples in at least one contrast group (biological replicates must be >= 2 per group)",
        n_control_samples = n_control, n_treatment_samples = n_treat
      )
      next
    }

    cm <- counts[, sub$agg_id, drop = FALSE]
    col_data <- data.frame(condition = factor(sub[[opt$condition_col]], levels = c(control, cond)),
                           row.names = sub$agg_id)
    rownames(col_data) <- sub$agg_id
    colnames(cm) <- sub$agg_id

    res <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = cm,
        colData = col_data,
        design = as.formula(paste("~", opt$condition_col))
      )
      dds <- DESeq(dds, quiet = TRUE)
      shrink_res(dds, cond)
    }, error = function(e) {
      results_log[[length(results_log) + 1]] <<- list(
        cell_type = ct, contrast = contrast_id, status = "ERROR",
        message = paste0("DESeq2 failed: ", conditionMessage(e)),
        n_control_samples = n_control, n_treatment_samples = n_treat
      )
      NULL
    })

    if (is.null(res)) next

    res_df <- data.frame(
      gene = rownames(res),
      log2fc = res$log2FoldChange,
      pval = res$pvalue,
      padj = res$padj,
      stringsAsFactors = FALSE
    )
    res_df <- res_df[!is.na(res_df$log2fc) & is.finite(res_df$log2fc), ]
    fwrite(res_df, file.path(opt$outdir, paste0(file_id, "_results.csv")))

    plot_volcano(res_df, contrast_id, file_id, opt$outdir)
    n_run <- n_run + 1L
    results_log[[length(results_log) + 1]] <- list(
      cell_type = ct, contrast = contrast_id, status = "OK",
      message = sprintf("DESeq2 completed on %d genes (%d control / %d treatment samples)",
                        nrow(res_df), n_control, n_treat),
      n_control_samples = n_control, n_treatment_samples = n_treat
    )
  }
}

logs_df <- do.call(rbind, lapply(results_log, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
fwrite(logs_df, file.path(opt$outdir, "pseudobulk_summary.csv"))

status <- make_status(skipped == 0 && n_run > 0, sprintf("Pseudobulk DGE done: %d contrasts run, %d skipped", n_run, skipped), results_log)
writeLines(toJSON(status, pretty = TRUE, auto_unbox = TRUE), file.path(opt$outdir, "pseudobulk_status.json"))

message(sprintf("Done: %d run, %d skipped", n_run, skipped))