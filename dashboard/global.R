# scAnnex Dashboard - Global Settings and Functions
# Loaded once when the app starts

# ==============================================================================
# Python Setup - MUST BE FIRST (before loading reticulate)
# ==============================================================================

# Configure Python path BEFORE loading reticulate library
# This ensures we use the conda environment's Python with scanpy/anndata
conda_prefix <- Sys.getenv("CONDA_PREFIX")

if (conda_prefix != "") {
  # We're in a conda environment, use its Python
  python_path <- file.path(conda_prefix, "bin", "python3")
  message(sprintf("Using Conda Python from CONDA_PREFIX: %s", python_path))
} else {
  # Try to find the scannex-dashboard conda environment
  possible_paths <- c(
    "/home/damo/miniforge3/envs/scannex-dashboard/bin/python3",
    "~/miniforge3/envs/scannex-dashboard/bin/python3",
    "/opt/conda/envs/scannex-dashboard/bin/python3",
    "/usr/local/miniconda3/envs/scannex-dashboard/bin/python3"
  )
  
  # Expand ~ to home directory
  possible_paths <- sapply(possible_paths, path.expand)
  
  # Find first existing path
  python_path <- NULL
  for (path in possible_paths) {
    if (file.exists(path)) {
      python_path <- path
      message(sprintf("Using Conda Python from auto-detected path: %s", python_path))
      break
    }
  }
  
  if (is.null(python_path)) {
    python_path <- "/usr/bin/python3"
    warning("Could not find scannex-dashboard conda environment. Using system Python (may not work).")
  }
}

# CRITICAL: Set this BEFORE loading reticulate
Sys.setenv(RETICULATE_PYTHON = python_path)

# ==============================================================================
# Load R Libraries
# ==============================================================================

# Suppress all package startup messages for clean output
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(plotly)
  library(DT)
  library(ggplot2)
  library(reticulate)  # Loaded AFTER setting RETICULATE_PYTHON
  library(viridis)
  library(data.table)
  library(jsonlite)
})

# ==============================================================================
# Finalize Python Configuration
# ==============================================================================

use_python(python_path, required = TRUE)

# Import required Python modules (suppress warnings)
suppressWarnings({
  ad <- import("anndata")
  sc <- import("scanpy")
  np <- import("numpy")
  pd <- import("pandas")
})

# ==============================================================================
# Data Loading Functions
# ==============================================================================

#' Load H5AD file with backed mode for large datasets
#' 
#' @param h5ad_path Path to .h5ad file
#' @param backed Logical, use backed mode? (default: TRUE for >50k cells)
#' @return List containing metadata and UMAP coordinates
load_h5ad_data <- function(h5ad_path, backed = TRUE) {
  
  message(sprintf("Loading H5AD file: %s", h5ad_path))
  
  # For small datasets (<100k cells), disable backed mode to access .obsm
  # Check file size first to make smart decision
  file_size_mb <- file.info(h5ad_path)$size / (1024^2)
  
  if (backed && file_size_mb < 500) {
    message(sprintf("  File size: %.1f MB - disabling backed mode for full .obsm access", file_size_mb))
    backed <- FALSE
  }
  
  # Read H5AD with backed mode for large datasets
  if (backed) {
    message("  Using backed mode (read-only, low memory)")
    adata <- ad$read_h5ad(h5ad_path, backed = "r")
  } else {
    message("  Loading into memory (full access)")
    adata <- ad$read_h5ad(h5ad_path)
  }
  
  # Extract metadata as data.frame
  # Use py_run_string to convert in Python and avoid reticulate py_to_r() issues
  message("  Converting metadata to R data.frame...")
  
  tryCatch({
    # Assign adata to Python namespace so py_run_string can access it
    py$adata <- adata
    
    # Use py_run_string to convert data in Python before passing to R
    # This avoids pandas DataFrame conversion issues in reticulate
    py_run_string("
# Convert obs to simple Python dict of lists
obs_dict = {}
for col in adata.obs.columns:
    obs_dict[col] = adata.obs[col].tolist()

obs_dict['cell_id'] = adata.obs_names.tolist()
")
    
    # Access the dict from R (py$ namespace)
    obs_dict <- py$obs_dict
    
    # Create data.frame from the dict
    metadata <- as.data.frame(obs_dict, stringsAsFactors = FALSE, check.names = FALSE)
    
    message(sprintf("  ✓ Metadata loaded: %d cells, %d columns", nrow(metadata), ncol(metadata)))
    
  }, error = function(e) {
    stop(sprintf("Failed to load metadata: %s", e$message))
  })
  
  # Extract UMAP coordinates if available
  umap_coords <- NULL
  
  tryCatch({
    message("  Checking for UMAP coordinates in .obsm['X_umap']...")
    
    # For backed mode, load UMAP separately since .obsm is not accessible
    if (backed) {
      message("  Backed mode detected - loading UMAP from separate read...")
      adata_temp <- ad$read_h5ad(h5ad_path, backed = NULL)
      umap_matrix <- py_to_r(adata_temp$obsm["X_umap"])
    } else {
      # Direct access for in-memory
      umap_matrix <- py_to_r(adata$obsm["X_umap"])
    }
    
    # Successfully got UMAP data
    # Use cell_id column from metadata, not rownames (which are just 1,2,3...)
    umap_coords <- data.frame(
      cell_id = metadata$cell_id,
      UMAP_1 = umap_matrix[, 1],
      UMAP_2 = umap_matrix[, 2]
    )
    
    message(sprintf("  ✓ UMAP coordinates loaded: %d cells", nrow(umap_coords)))
    message(sprintf("  ✓ UMAP range: [%.2f, %.2f] to [%.2f, %.2f]",
                    min(umap_coords$UMAP_1), max(umap_coords$UMAP_1),
                    min(umap_coords$UMAP_2), max(umap_coords$UMAP_2)))
    
  }, error = function(e) {
    message(sprintf("  WARNING: Could not load UMAP coordinates: %s", e$message))
    message("  Plots requiring UMAP will not be available.")
  })
  
  # Get variable genes info
  # Use same py_run_string approach for var DataFrame
  tryCatch({
    py_run_string("
# Convert var to simple Python dict of lists
var_dict = {}
for col in adata.var.columns:
    var_dict[col] = adata.var[col].tolist()

var_dict['gene_id'] = adata.var_names.tolist()
")
    
    # Access the dict from R
    var_dict <- py$var_dict
    var_info <- as.data.frame(var_dict, stringsAsFactors = FALSE, check.names = FALSE)
    
    # Set gene names as rownames for easy lookup
    rownames(var_info) <- var_info$gene_id
    
  }, error = function(e) {
    message(sprintf("  WARNING: Could not load var info: %s", e$message))
    var_info <- NULL
  })
  
  # Get dataset dimensions
  n_cells <- adata$n_obs
  n_genes <- adata$n_vars
  
  message(sprintf("  Loaded: %d cells × %d genes", n_cells, n_genes))
  
  return(list(
    adata = adata,
    metadata = metadata,
    umap_coords = umap_coords,
    var_info = var_info,
    n_cells = n_cells,
    n_genes = n_genes,
    backed = backed
  ))
}

#' Extract gene expression for a specific gene
#' 
#' @param data_obj Data object from load_h5ad_data()
#' @param gene_name Gene name to extract
#' @return Numeric vector of expression values
get_gene_expression <- function(data_obj, gene_name) {
  
  adata <- data_obj$adata
  
  # Check if gene exists
  if (!(gene_name %in% rownames(data_obj$var_info))) {
    stop(sprintf("Gene '%s' not found in dataset", gene_name))
  }
  
  # Get gene index
  gene_idx <- which(rownames(data_obj$var_info) == gene_name) - 1  # Python 0-indexed
  
  # Extract expression (handles both backed and in-memory)
  if (data_obj$backed) {
    expr <- py_to_r(adata$X[, as.integer(gene_idx)]$toarray()$flatten())
  } else {
    expr <- py_to_r(adata$X[, as.integer(gene_idx)])
  }
  
  names(expr) <- data_obj$metadata$cell_id
  return(expr)
}

#' Load QC report JSON
#' 
#' @param qc_dir Path to QC results directory
#' @return List containing QC metrics and thresholds
load_qc_report <- function(qc_dir) {
  
  # Try multiple possible locations for qc_report.json
  possible_paths <- c(
    file.path(qc_dir, "qc_report.json"),             # Direct path
    file.path(qc_dir, "qc_results", "qc_report.json"), # In qc_results subdirectory (CURRENT)
    file.path(qc_dir, "results", "qc_report.json"),  # In results subdirectory (OLD)
    file.path(dirname(qc_dir), "qc", "qc_report.json")  # Parent qc dir
  )
  
  report_path <- NULL
  for (path in possible_paths) {
    if (file.exists(path)) {
      report_path <- path
      break
    }
  }
  
  if (is.null(report_path)) {
    warning(sprintf("QC report not found in: %s", qc_dir))
    warning(sprintf("  Tried: %s", paste(possible_paths, collapse=", ")))
    return(NULL)
  }
  
  message(sprintf("Loading QC report from: %s", report_path))
  qc_report <- fromJSON(report_path)
  
  return(qc_report)
}

#' Get available QC plots
#' 
#' @param qc_dir Path to QC results directory
#' @return Character vector of plot file paths
get_qc_plots <- function(qc_dir) {
  
  if (!dir.exists(qc_dir)) {
    return(character(0))
  }
  
  # Try multiple possible locations for plots
  possible_dirs <- c(
    qc_dir,                                  # Direct path
    file.path(qc_dir, "qc_results"),         # In qc_results subdirectory (CURRENT)
    file.path(qc_dir, "plots")               # In plots subdirectory
  )
  
  plot_files <- character(0)
  
  for (plot_dir in possible_dirs) {
    if (dir.exists(plot_dir)) {
      files <- list.files(
        plot_dir, 
        pattern = "\\.(png|jpg|jpeg)$", 
        full.names = TRUE
      )
      plot_files <- c(plot_files, files)
    }
  }
  
  # Remove duplicates and sort
  plot_files <- unique(plot_files)
  plot_files <- sort(plot_files)
  
  if (length(plot_files) > 0) {
    message(sprintf("Found %d QC plot(s)", length(plot_files)))
  }
  
  return(plot_files)
}

# ==============================================================================
# Plotting Functions
# ==============================================================================

#' Create interactive UMAP plot with plotly
#' 
#' @param umap_data Data frame with UMAP_1, UMAP_2, and color variable
#' @param color_by Column name to color points by
#' @param title Plot title
#' @return Plotly object
plot_umap_interactive <- function(umap_data, color_by, title = "UMAP") {
  
  if (is.null(umap_data)) {
    return(NULL)
  }
  
  # Merge metadata if color_by not in umap_data
  # (handled in server logic)
  
  p <- plot_ly(
    data = umap_data,
    x = ~UMAP_1,
    y = ~UMAP_2,
    type = 'scattergl',  # Use WebGL for performance
    mode = 'markers',
    marker = list(
      size = 3,
      opacity = 0.7
    ),
    text = ~paste("Cell:", cell_id),
    hoverinfo = 'text'
  )
  
  # Add color if specified
  if (!is.null(color_by) && color_by %in% names(umap_data)) {
    p <- p %>%
      add_trace(color = as.formula(paste0("~", color_by)))
  }
  
  p <- p %>%
    layout(
      title = title,
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2"),
      hovermode = 'closest'
    )
  
  return(p)
}

#' Create QC violin plot
#' 
#' @param metadata Data frame with QC metrics
#' @param metric_col Column name of metric to plot
#' @param threshold Numeric threshold to draw line (optional)
#' @return ggplot object
plot_qc_violin <- function(metadata, metric_col, threshold = NULL) {
  
  p <- ggplot(metadata, aes_string(x = "1", y = metric_col)) +
    geom_violin(fill = "lightblue", alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
    theme_minimal() +
    labs(
      title = sprintf("Distribution of %s", metric_col),
      x = "",
      y = metric_col
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (!is.null(threshold)) {
    p <- p + geom_hline(
      yintercept = threshold, 
      color = "red", 
      linetype = "dashed"
    )
  }
  
  return(p)
}

#' Calculate gene set score using ranking-based method (AUCell-inspired)
#' 
#' @param data_obj Data object from load_h5ad_data()
#' @param gene_list Character vector of gene names
#' @return Named numeric vector of normalized scores per cell (0-1 scale)
calculate_gene_set_score <- function(data_obj, gene_list) {
  
  adata <- data_obj$adata
  
  # Clean and validate gene list
  gene_list <- unique(trimws(gene_list))
  gene_list <- gene_list[gene_list != ""]
  
  if (length(gene_list) == 0) {
    stop("Gene list is empty")
  }
  
  # Find which genes exist in the dataset
  available_genes <- rownames(data_obj$var_info)
  found_genes <- gene_list[gene_list %in% available_genes]
  missing_genes <- gene_list[!(gene_list %in% available_genes)]
  
  if (length(found_genes) == 0) {
    stop(sprintf("None of the genes found in dataset. Missing: %s", 
                 paste(missing_genes, collapse = ", ")))
  }
  
  if (length(missing_genes) > 0) {
    message(sprintf("  Warning: %d gene(s) not found: %s", 
                    length(missing_genes),
                    paste(missing_genes, collapse = ", ")))
  }
  
  message(sprintf("  Calculating ranking-based score using %d gene(s)", 
                  length(found_genes)))
  
  # Extract expression for all genes in the set
  gene_expr_list <- list()
  
  for (gene in found_genes) {
    gene_idx <- which(rownames(data_obj$var_info) == gene) - 1  # Python 0-indexed
    
    # Extract expression (handles both backed and in-memory)
    if (data_obj$backed) {
      expr <- py_to_r(adata$X[, as.integer(gene_idx)]$toarray()$flatten())
    } else {
      expr <- py_to_r(adata$X[, as.integer(gene_idx)])
      if (is.matrix(expr)) {
        expr <- as.vector(expr)
      }
    }
    
    gene_expr_list[[gene]] <- expr
  }
  
  # Combine into matrix (genes x cells)
  expr_matrix <- do.call(rbind, gene_expr_list)
  rownames(expr_matrix) <- found_genes
  
  # Rank-based scoring (simplified AUCell approach)
  # For each cell: calculate mean percentile rank of signature genes
  # This gives a score from 0 (low expression) to 1 (high expression)
  
  n_cells <- ncol(expr_matrix)
  scores <- numeric(n_cells)
  
  for (i in 1:n_cells) {
    # Get expression values for this cell's signature genes
    cell_sig_expr <- expr_matrix[, i]
    
    # Calculate percentile rank within the signature genes
    # (comparing to other genes in the signature for this cell)
    # Higher expression = higher percentile
    ranks <- rank(cell_sig_expr, ties.method = "average")
    mean_rank <- mean(ranks)
    
    # Normalize to 0-1 scale
    n_genes_sig <- length(cell_sig_expr)
    normalized_score <- (mean_rank - 1) / (n_genes_sig - 1)
    
    scores[i] <- normalized_score
  }
  
  # Alternative: use mean expression but normalize to 0-1 scale
  # This is simpler and faster while still being meaningful
  mean_expr <- colMeans(expr_matrix)
  
  # Normalize to 0-1 range (min-max scaling)
  min_expr <- min(mean_expr)
  max_expr <- max(mean_expr)
  
  if (max_expr > min_expr) {
    normalized_scores <- (mean_expr - min_expr) / (max_expr - min_expr)
  } else {
    normalized_scores <- rep(0.5, length(mean_expr))  # All same value
  }
  
  names(normalized_scores) <- data_obj$metadata$cell_id
  
  message(sprintf("  Score range: %.3f - %.3f (mean: %.3f)", 
                  min(normalized_scores), max(normalized_scores), mean(normalized_scores)))
  
  return(normalized_scores)
}

#' Parse custom annotation rules from textarea input
#' 
#' @param rules_text Text with annotation rules (one per line)
#' Format: clustering_name,cluster_id,label
#' Example: leiden_0.5,4,HSC
#' @return Data frame with columns: clustering, cluster, label
parse_annotation_rules <- function(rules_text) {
  
  # Split by newlines and clean
  lines <- trimws(strsplit(rules_text, "\n")[[1]])
  lines <- lines[lines != ""]
  
  if (length(lines) == 0) {
    return(data.frame(
      clustering = character(),
      cluster = character(),
      label = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Parse each line
  rules_list <- list()
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # Split by comma
    parts <- trimws(strsplit(line, ",")[[1]])
    
    if (length(parts) != 3) {
      warning(sprintf("Skipping invalid line %d: '%s' (expected 3 comma-separated values)", i, line))
      next
    }
    
    rules_list[[length(rules_list) + 1]] <- data.frame(
      clustering = parts[1],
      cluster = parts[2],
      label = parts[3],
      stringsAsFactors = FALSE
    )
  }
  
  if (length(rules_list) == 0) {
    return(data.frame(
      clustering = character(),
      cluster = character(),
      label = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Combine into data frame
  rules_df <- do.call(rbind, rules_list)
  
  return(rules_df)
}

#' Apply custom annotation rules to create a new annotation column
#' 
#' @param metadata Data frame with cell metadata (must include cell_id)
#' @param rules_df Data frame with annotation rules (from parse_annotation_rules)
#' @param annotation_name Name for the new annotation column
#' @return Named character vector with cell_id as names and labels as values
apply_custom_annotation <- function(metadata, rules_df, annotation_name = "custom_annotation") {
  
  if (nrow(rules_df) == 0) {
    stop("No valid annotation rules provided")
  }
  
  # Initialize all cells as "Unknown"
  n_cells <- nrow(metadata)
  labels <- rep("Unknown", n_cells)
  names(labels) <- metadata$cell_id
  
  message(sprintf("Applying custom annotation '%s' with %d rule(s)...", annotation_name, nrow(rules_df)))
  
  # Apply rules in order (later rules override earlier ones)
  for (i in 1:nrow(rules_df)) {
    rule <- rules_df[i, ]
    
    clustering_col <- rule$clustering
    cluster_id <- rule$cluster
    new_label <- rule$label
    
    # Check if clustering column exists
    if (!(clustering_col %in% names(metadata))) {
      warning(sprintf("Clustering column '%s' not found in metadata. Skipping rule: %s,%s,%s",
                      clustering_col, clustering_col, cluster_id, new_label))
      next
    }
    
    # Find cells matching this cluster
    matching_cells <- metadata$cell_id[metadata[[clustering_col]] == cluster_id]
    
    if (length(matching_cells) == 0) {
      warning(sprintf("No cells found for cluster '%s' in clustering '%s'", 
                      cluster_id, clustering_col))
      next
    }
    
    # Apply label (overwrites previous labels)
    labels[matching_cells] <- new_label
    
    message(sprintf("  Rule %d: %s = %s → '%s' (%d cells)",
                    i, clustering_col, cluster_id, new_label, length(matching_cells)))
  }
  
  # Count labels
  label_counts <- table(labels)
  unknown_count <- sum(labels == "Unknown")
  
  message(sprintf("Custom annotation completed:"))
  message(sprintf("  Total cells: %d", n_cells))
  message(sprintf("  Annotated: %d (%.1f%%)", n_cells - unknown_count, 
                  100 * (n_cells - unknown_count) / n_cells))
  message(sprintf("  Unknown: %d (%.1f%%)", unknown_count, 100 * unknown_count / n_cells))
  message(sprintf("  Unique labels: %d", length(unique(labels))))
  
  return(labels)
}

#' Save custom annotation to H5AD file
#' 
#' @param h5ad_path Path to H5AD file
#' @param annotation_name Name of the annotation column
#' @param labels Named vector of labels (cell_id as names)
#' @param output_path Output path for H5AD (NULL = overwrite original)
#' @param create_copy Logical, create a copy instead of overwriting
save_annotation_to_h5ad <- function(h5ad_path, annotation_name, labels,
                                     output_path = NULL, create_copy = FALSE) {
  
  message(sprintf("Saving annotation '%s' to H5AD...", annotation_name))
  
  # Get cell IDs from labels (we already have them as names)
  cell_ids <- names(labels)
  labels_vec <- as.character(unname(labels))
  
  # Write labels to temporary CSV file
  temp_csv <- tempfile(fileext = ".csv")
  write.table(
    data.frame(cell_id = cell_ids, label = labels_vec),
    file = temp_csv,
    sep = ",",
    row.names = FALSE,
    quote = TRUE
  )
  
  # Determine output path
  if (is.null(output_path)) {
    if (create_copy) {
      # Create versioned copy
      base_name <- tools::file_path_sans_ext(h5ad_path)
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      output_path <- sprintf("%s_%s.h5ad", base_name, timestamp)
    } else {
      # Overwrite original
      output_path <- h5ad_path
    }
  }
  
  # Do everything in Python to avoid reticulate conversion issues
  py_run_string(sprintf("
import anndata as ad
import pandas as pd

# Allow writing nullable strings (needed for anndata >= 0.11)
ad.settings.allow_write_nullable_strings = True

# Load H5AD file
adata = ad.read_h5ad('%s')

# Read labels from CSV
labels_df = pd.read_csv('%s', index_col='cell_id')

# Ensure labels are in same order as adata.obs
labels_ordered = labels_df.loc[adata.obs_names, 'label']

# Add to adata.obs as Categorical
adata.obs['%s'] = pd.Categorical(labels_ordered.values)

# Save H5AD
adata.write_h5ad('%s')

# Return success
success = True
", h5ad_path, temp_csv, annotation_name, output_path))
  
  # Check if Python operation succeeded
  if (!py$success) {
    unlink(temp_csv)
    stop("Failed to save annotation in Python")
  }
  
  # Clean up temp file
  unlink(temp_csv)
  
  message(sprintf("  Writing to: %s", output_path))
  message(sprintf("  Successfully saved annotation '%s' with %d unique labels", 
                  annotation_name, length(unique(labels))))
  
  return(output_path)
}

# ==============================================================================
# Utility Functions
# ==============================================================================

#' Format large numbers with commas
#' @param x Numeric value
#' @return Character string with formatted number
format_number <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

#' Calculate percentage
#' @param part Numerator
#' @param total Denominator
#' @return Character string with percentage
format_percentage <- function(part, total) {
  sprintf("%.1f%%", (part / total) * 100)
}

# ==============================================================================
# Default Data Path (can be overridden via environment variable)
# ==============================================================================

# Set default data path
data_path_env <- Sys.getenv("SCANNEX_DATA_PATH")
DEFAULT_DATA_PATH <- if (data_path_env == "") {
  "/srv/shiny-server/data"
} else {
  data_path_env
}

# Auto-detect H5AD file and QC directory
DEFAULT_H5AD_FILE <- ""
DEFAULT_QC_DIR <- ""

if (dir.exists(DEFAULT_DATA_PATH)) {
  # Find H5AD file - try dashboard, annotated, or processed files
  patterns <- c(
    ".*dashboard.*\\.h5ad$",     # Dashboard-compatible files (priority)
    ".*annotated.*\\.h5ad$",      # Annotated files
    ".*processed.*\\.h5ad$"       # Processed files
  )
  
  h5ad_files <- character(0)
  for (pattern in patterns) {
    files <- list.files(
      DEFAULT_DATA_PATH,
      pattern = pattern,
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(files) > 0) {
      h5ad_files <- files
      break
    }
  }
  
  if (length(h5ad_files) > 0) {
    DEFAULT_H5AD_FILE <- h5ad_files[1]
    message(sprintf("  Auto-detected H5AD: %s", basename(DEFAULT_H5AD_FILE)))
  }
  
  # Find QC directory
  qc_dirs <- list.files(
    DEFAULT_DATA_PATH,
    pattern = "^qc$",
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = TRUE
  )
  
  if (length(qc_dirs) > 0) {
    DEFAULT_QC_DIR <- qc_dirs[1]
    message(sprintf("  Auto-detected QC dir: %s", DEFAULT_QC_DIR))
  }
}

message("=============================================================")
message("scAnnex Dashboard - Global environment initialized")
message(sprintf("Default data path: %s", DEFAULT_DATA_PATH))
message("=============================================================")
