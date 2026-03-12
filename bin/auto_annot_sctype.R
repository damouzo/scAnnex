#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(jsonlite)
    library(Seurat)
})

print_help <- function() {
    cat("Usage: Rscript auto_annot_sctype.R --input <in.rds> --output <out.csv> --status <status.json> --markers-file <markers.csv> [--continue-on-error]\n")
}

parse_args_simple <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    opt <- list(
        input = NULL,
        output = NULL,
        status = NULL,
        `markers-file` = NULL,
        `continue-on-error` = FALSE
    )

    i <- 1L
    while (i <= length(argv)) {
        arg <- argv[[i]]
        if (arg %in% c("-h", "--help")) {
            print_help()
            quit(status = 0)
        } else if (arg %in% c("--input", "--output", "--status", "--markers-file")) {
            i <- i + 1L
            if (i > length(argv)) stop(sprintf("Missing value for %s", arg), call. = FALSE)
            key <- sub("^--", "", arg)
            opt[[key]] <- argv[[i]]
        } else if (arg == "--continue-on-error") {
            opt$`continue-on-error` <- TRUE
        } else {
            stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
        }
        i <- i + 1L
    }

    if (is.null(opt$input) || is.null(opt$output) || is.null(opt$status) || is.null(opt$`markers-file`)) {
        print_help()
        stop("Missing required arguments", call. = FALSE)
    }

    opt
}

normalize_marker_table <- function(markers_df) {
    colnames(markers_df) <- tolower(trimws(colnames(markers_df)))
    required <- c("cell_type", "gene")
    missing <- setdiff(required, colnames(markers_df))
    if (length(missing) > 0) {
        stop(sprintf("Markers file missing required columns: %s", paste(missing, collapse = ", ")))
    }

    markers_df$cell_type <- as.character(markers_df$cell_type)
    markers_df$gene <- as.character(markers_df$gene)
    if (!"weight" %in% colnames(markers_df)) {
        markers_df$weight <- 1.0
    }
    markers_df$weight <- suppressWarnings(as.numeric(markers_df$weight))
    markers_df$weight[is.na(markers_df$weight)] <- 1.0
    markers_df
}

run_sctype_like <- function(obj, markers_file) {
    if (!file.exists(markers_file)) {
        stop(sprintf("Markers file not found: %s", markers_file))
    }

    markers_raw <- read.csv(markers_file, stringsAsFactors = FALSE, check.names = FALSE)
    markers_tbl <- normalize_marker_table(markers_raw)

    if (nrow(markers_tbl) == 0) {
        stop("Markers file is empty")
    }

    if ("data" %in% Layers(obj[[DefaultAssay(obj)]])) {
        expr <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")
    } else {
        expr <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
    }

    genes <- rownames(expr)
    cells <- colnames(expr)
    labels <- rep("Unknown", length(cells))
    names(labels) <- cells
    scores <- rep(NA_real_, length(cells))
    names(scores) <- cells

    by_type <- split(markers_tbl, markers_tbl$cell_type)
    score_matrix <- matrix(NA_real_, nrow = length(by_type), ncol = length(cells))
    rownames(score_matrix) <- names(by_type)
    colnames(score_matrix) <- cells

    row_idx <- 1L
    for (ct in names(by_type)) {
        marker_set <- by_type[[ct]]
        marker_genes <- marker_set$gene
        marker_weights <- marker_set$weight
        keep <- marker_genes %in% genes

        if (!any(keep)) {
            row_idx <- row_idx + 1L
            next
        }

        marker_genes <- marker_genes[keep]
        marker_weights <- marker_weights[keep]
        sub_expr <- expr[marker_genes, , drop = FALSE]

        if (inherits(sub_expr, "dgCMatrix")) {
            weighted <- Matrix::Diagonal(x = marker_weights) %*% sub_expr
            ct_score <- Matrix::colMeans(weighted)
        } else {
            ct_score <- colMeans(sweep(as.matrix(sub_expr), 1, marker_weights, `*`))
        }

        score_matrix[row_idx, ] <- as.numeric(ct_score)
        row_idx <- row_idx + 1L
    }

    if (nrow(score_matrix) > 0) {
        max_idx <- apply(score_matrix, 2, function(x) {
            if (all(is.na(x))) return(NA_integer_)
            which.max(x)
        })

        for (i in seq_along(cells)) {
            idx <- max_idx[[i]]
            if (!is.na(idx)) {
                labels[[i]] <- rownames(score_matrix)[idx]
                scores[[i]] <- score_matrix[idx, i]
            }
        }
    }

    meta <- obj@meta.data
    cluster_col <- NULL
    if ("seurat_clusters" %in% colnames(meta)) {
        cluster_col <- "seurat_clusters"
    } else if ("leiden_0.5" %in% colnames(meta)) {
        cluster_col <- "leiden_0.5"
    } else if (ncol(meta) > 0) {
        cluster_col <- colnames(meta)[1]
    }

    cluster_vals <- if (!is.null(cluster_col)) as.character(meta[colnames(obj), cluster_col]) else rep(NA_character_, length(cells))

    data.frame(
        cell_id = cells,
        auto_annot_sctype = as.character(labels),
        auto_annot_sctype_score = as.numeric(scores),
        auto_annot_sctype_cluster = cluster_vals,
        stringsAsFactors = FALSE
    )
}

opt <- parse_args_simple()
status <- list(tool = "sctype", success = TRUE, errors = list())

output_df <- NULL
tryCatch({
    obj <- readRDS(opt$input)
    output_df <- run_sctype_like(obj, opt$`markers-file`)
}, error = function(e) {
    status$success <- FALSE
    status$errors[[1]] <- e$message

    fallback_cells <- character(0)
    obj_try <- try(readRDS(opt$input), silent = TRUE)
    if (!inherits(obj_try, "try-error")) {
        fallback_cells <- colnames(obj_try)
    }

    output_df <<- data.frame(
        cell_id = fallback_cells,
        auto_annot_sctype = NA_character_,
        auto_annot_sctype_score = NA_real_,
        auto_annot_sctype_cluster = NA_character_,
        stringsAsFactors = FALSE
    )

    if (!opt$`continue-on-error`) {
        write.csv(output_df, opt$output, row.names = FALSE, quote = TRUE)
        write_json(status, opt$status, pretty = TRUE, auto_unbox = TRUE)
        quit(status = 1)
    }
})

write.csv(output_df, opt$output, row.names = FALSE, quote = TRUE)
write_json(status, opt$status, pretty = TRUE, auto_unbox = TRUE)
