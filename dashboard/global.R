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
  message(sprintf("Using Conda Python: %s", python_path))
} else {
  # Fallback to system Python (may not work without packages)
  python_path <- "/usr/bin/python3"
  warning("Not in Conda environment. Python packages may be missing.")
}

# CRITICAL: Set this BEFORE loading reticulate
Sys.setenv(RETICULATE_PYTHON = python_path)

# ==============================================================================
# Load R Libraries
# ==============================================================================

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

# ==============================================================================
# Finalize Python Configuration
# ==============================================================================

use_python(python_path, required = TRUE)

# Import required Python modules
ad <- import("anndata")
sc <- import("scanpy")
np <- import("numpy")
pd <- import("pandas")

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
  metadata <- py_to_r(adata$obs)
  metadata$cell_id <- rownames(metadata)
  
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
    umap_coords <- data.frame(
      cell_id = rownames(metadata),
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
  var_info <- py_to_r(adata$var)
  
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
  
  report_path <- file.path(qc_dir, "qc_report.json")
  
  if (!file.exists(report_path)) {
    warning(sprintf("QC report not found: %s", report_path))
    return(NULL)
  }
  
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
  
  plot_files <- list.files(
    qc_dir, 
    pattern = "\\.(png|jpg|jpeg)$", 
    full.names = TRUE
  )
  
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

message("=============================================================")
message("scAnnex Dashboard - Global environment initialized")
message(sprintf("Default data path: %s", DEFAULT_DATA_PATH))
message("=============================================================")
