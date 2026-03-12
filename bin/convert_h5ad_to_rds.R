#!/usr/bin/env Rscript

#==============================================================================
# H5AD to RDS Conversion Script
#==============================================================================
# Converts AnnData H5AD files to Seurat v5 RDS objects for R-based annotation
# tools (e.g., Azimuth, SingleR)
#
# Requirements:
#   - Seurat >= 5.0.0
#   - SeuratDisk >= 0.0.0.9020
#   - reticulate
#   - anndata (Python)
#
# Usage:
#   Rscript convert_h5ad_to_rds.R --input integrated.h5ad --output seurat.rds
#==============================================================================

suppressPackageStartupMessages({
    library(Seurat)
})

# Prevent reticulate from creating managed uv environments inside HPC scratch/cache
Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")

python_bin <- Sys.which("python3")
if (python_bin == "") {
    python_bin <- Sys.which("python")
}
if (python_bin != "") {
    Sys.setenv(RETICULATE_PYTHON = python_bin)
}

suppressPackageStartupMessages({
    library(reticulate)
})

#==============================================================================
# Command Line Arguments
#==============================================================================
print_help <- function() {
    cat("Convert H5AD (AnnData) to Seurat RDS object\n\n")
    cat("Usage:\n")
    cat("  Rscript convert_h5ad_to_rds.R --input integrated.h5ad --output seurat.rds [options]\n\n")
    cat("Options:\n")
    cat("  -i, --input FILE          Input H5AD file path (required)\n")
    cat("  -o, --output FILE         Output RDS file path [default: seurat_object.rds]\n")
    cat("      --assay-name STRING   Name for Seurat assay [default: RNA]\n")
    cat("      --use-raw             Use .raw.X from H5AD instead of .X\n")
    cat("      --min-cells INT       Minimum cells for feature filtering [default: 0]\n")
    cat("      --min-features INT    Minimum features for cell filtering [default: 0]\n")
    cat("      --verbose             Print detailed progress messages\n")
    cat("  -h, --help                Show this help\n")
}

parse_args_simple <- function() {
    argv <- commandArgs(trailingOnly = TRUE)

    opt <- list(
        input = NULL,
        output = "seurat_object.rds",
        `assay-name` = "RNA",
        `use-raw` = FALSE,
        `min-cells` = 0L,
        `min-features` = 0L,
        verbose = FALSE
    )

    i <- 1L
    while (i <= length(argv)) {
        arg <- argv[[i]]

        if (arg %in% c("-h", "--help")) {
            print_help()
            quit(status = 0)
        } else if (arg %in% c("-i", "--input")) {
            i <- i + 1L
            if (i > length(argv)) stop("Missing value for --input", call. = FALSE)
            opt$input <- argv[[i]]
        } else if (arg %in% c("-o", "--output")) {
            i <- i + 1L
            if (i > length(argv)) stop("Missing value for --output", call. = FALSE)
            opt$output <- argv[[i]]
        } else if (arg == "--assay-name") {
            i <- i + 1L
            if (i > length(argv)) stop("Missing value for --assay-name", call. = FALSE)
            opt$`assay-name` <- argv[[i]]
        } else if (arg == "--use-raw") {
            opt$`use-raw` <- TRUE
        } else if (arg == "--min-cells") {
            i <- i + 1L
            if (i > length(argv)) stop("Missing value for --min-cells", call. = FALSE)
            opt$`min-cells` <- as.integer(argv[[i]])
        } else if (arg == "--min-features") {
            i <- i + 1L
            if (i > length(argv)) stop("Missing value for --min-features", call. = FALSE)
            opt$`min-features` <- as.integer(argv[[i]])
        } else if (arg == "--verbose") {
            opt$verbose <- TRUE
        } else {
            stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
        }

        i <- i + 1L
    }

    opt
}

opt <- parse_args_simple()

# Validate required arguments
if (is.null(opt$input)) {
    print_help()
    stop("Input H5AD file is required (--input)", call. = FALSE)
}

if (!file.exists(opt$input)) {
    stop(sprintf("Input file does not exist: %s", opt$input), call. = FALSE)
}

#==============================================================================
# Helper Functions
#==============================================================================

log_message <- function(msg) {
    if (opt$verbose) {
        cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
    }
}

py_to_dense_matrix <- function(py_obj) {
    obj <- py_obj
    if (py_has_attr(obj, "toarray")) {
        obj <- obj$toarray()
    }
    mat <- py_to_r(obj)
    if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
    }
    mat
}

py_to_char_vector <- function(py_obj) {
    vals <- py_to_r(py_obj)
    as.character(unlist(vals, use.names = FALSE))
}



#==============================================================================
# Main Conversion Pipeline
#==============================================================================

log_message("Starting H5AD to RDS conversion")
log_message(sprintf("Input: %s", opt$input))
log_message(sprintf("Output: %s", opt$output))

# Step 1: Load H5AD using reticulate + anndata
log_message("Loading H5AD file using Python anndata...")

if (python_bin != "") {
    use_python(python_bin, required = TRUE)
}
if (!py_available(initialize = TRUE)) {
    stop("Python runtime not available in container", call. = FALSE)
}

ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(opt$input)
log_message(sprintf("Loaded AnnData: %d cells x %d genes",
    as.integer(py_to_r(adata$n_obs)), as.integer(py_to_r(adata$n_vars))))

# Step 2: Extract matrices from pipeline layers
# Schema (guaranteed after unify_input.py): layers['counts'] = raw counts, layers['log1p_norm'] = log-normalized
log_message("Extracting expression matrices...")

obs_names <- py_to_char_vector(adata$obs_names$to_list())
var_names <- py_to_char_vector(adata$var_names$to_list())

counts <- t(py_to_dense_matrix(adata$layers[["counts"]]))
rownames(counts) <- var_names
colnames(counts) <- obs_names

data_matrix <- t(py_to_dense_matrix(adata$layers[["log1p_norm"]]))
rownames(data_matrix) <- var_names
colnames(data_matrix) <- obs_names

log_message(sprintf("Matrix dimensions: %d genes x %d cells", length(var_names), length(obs_names)))

# Step 3: Extract metadata
log_message("Extracting cell metadata...")

meta_tmp <- tempfile(pattern = "adata_obs_", fileext = ".csv")
adata$obs$to_csv(meta_tmp, index = TRUE)
metadata <- read.csv(meta_tmp, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
unlink(meta_tmp)

if (nrow(metadata) != length(obs_names)) {
    stop(
        sprintf(
            "Metadata row count mismatch: %d rows in obs vs %d cells in matrix",
            nrow(metadata),
            length(obs_names)
        ),
        call. = FALSE
    )
}

metadata <- metadata[obs_names, , drop = FALSE]

log_message(sprintf("Metadata columns: %s", 
                   paste(colnames(metadata), collapse = ", ")))

# Step 4: Extract dimensional reductions (if available)
log_message("Checking for dimensional reductions...")

reductions <- list()
obsm_keys <- py_to_char_vector(adata$obsm_keys())

if ("X_pca" %in% obsm_keys) {
    log_message("Found PCA coordinates")
    pca_coords <- py_to_dense_matrix(adata$obsm[["X_pca"]])
    rownames(pca_coords) <- obs_names
    colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))
    reductions$pca <- pca_coords
}

if ("X_pca_harmony" %in% obsm_keys) {
    log_message("Found Harmony-corrected PCA coordinates")
    harmony_coords <- py_to_dense_matrix(adata$obsm[["X_pca_harmony"]])
    rownames(harmony_coords) <- obs_names
    colnames(harmony_coords) <- paste0("harmony_", 1:ncol(harmony_coords))
    reductions$harmony <- harmony_coords
}

if ("X_umap" %in% obsm_keys) {
    log_message("Found UMAP coordinates")
    umap_coords <- py_to_dense_matrix(adata$obsm[["X_umap"]])
    rownames(umap_coords) <- obs_names
    colnames(umap_coords) <- paste0("UMAP_", 1:ncol(umap_coords))
    reductions$umap <- umap_coords
}

# Step 5: Create Seurat object
log_message("Creating Seurat object...")

seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    assay = opt$`assay-name`,
    min.cells = opt$`min-cells`,
    min.features = opt$`min-features`
)

log_message(sprintf("Created Seurat object: %d cells x %d genes", ncol(seurat_obj), nrow(seurat_obj)))

for (reduction_name in names(reductions)) {
    log_message(sprintf("Adding %s reduction", reduction_name))
    reduction_coords <- reductions[[reduction_name]][colnames(seurat_obj), , drop = FALSE]
    seurat_obj[[reduction_name]] <- CreateDimReducObject(
        embeddings = reduction_coords,
        key = paste0(toupper(substring(reduction_name, 1, 1)), substring(reduction_name, 2), "_"),
        assay = opt$`assay-name`
    )
}

# Step 6: Add log-normalized data to the 'data' slot
log_message("Adding log-normalized data to Seurat 'data' layer")
seurat_obj <- SetAssayData(
    object = seurat_obj,
    assay = opt$`assay-name`,
    layer = "data",
    new.data = data_matrix[rownames(seurat_obj), colnames(seurat_obj)]
)

# Step 7: Save RDS
log_message(sprintf("Saving Seurat object to %s", opt$output))

saveRDS(seurat_obj, file = opt$output)
log_message("Conversion completed successfully!")

# Step 8: Print summary
cat("\n")
cat("==============================================\n")
cat("Conversion Summary\n")
cat("==============================================\n")
cat(sprintf("Input H5AD:        %s\n", opt$input))
cat(sprintf("Output RDS:        %s\n", opt$output))
cat(sprintf("Cells:             %d\n", ncol(seurat_obj)))
cat(sprintf("Genes:             %d\n", nrow(seurat_obj)))
cat(sprintf("Metadata columns:  %d\n", ncol(seurat_obj@meta.data)))
cat(sprintf("Reductions:        %s\n", 
           ifelse(length(reductions) > 0, 
                  paste(names(reductions), collapse = ", "), 
                  "None")))
cat("==============================================\n")

# Exit successfully
quit(status = 0)
