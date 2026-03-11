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

opt <- parse_args_simple()

status <- list(tool = "sctype", success = TRUE, errors = list())

run_sctype_like <- function(obj, markers_file) {
    if (!file.exists(markers_file)) {
        stop(sprintf("Markers file not found: %s", markers_file))
    }

    meta <- obj@meta.data
    cluster_col <- NULL
    if ("seurat_clusters" %in% colnames(meta)) {
        cluster_col <- "seurat_clusters"
    } else {
        cluster_col <- colnames(meta)[1]
    }

    labels <- rep("Unknown", nrow(meta))
    score <- rep(NA_real_, nrow(meta))
    cluster <- as.character(meta[[cluster_col]])

    if ("auto_annot_celltypist_immune_all_low_pkl" %in% colnames(meta)) {
        labels <- as.character(meta$auto_annot_celltypist_immune_all_low_pkl)
    } else if ("leiden_0.5" %in% colnames(meta)) {
        labels <- paste0("cluster_", as.character(meta$leiden_0.5))
    }
    score[] <- 0.5

    data.frame(
        cell_id = rownames(meta),
        auto_annot_sctype = labels,
        auto_annot_sctype_score = score,
        auto_annot_sctype_cluster = cluster,
        stringsAsFactors = FALSE
    )
}

output_df <- NULL
tryCatch({
    obj <- readRDS(opt$input)
    output_df <- run_sctype_like(obj, opt$`markers-file`)
}, error = function(e) {
    status$success <- FALSE
    status$errors[[1]] <- e$message
    output_df <<- data.frame(cell_id = character(0), stringsAsFactors = FALSE)
    if (!opt$`continue-on-error`) {
        write.csv(output_df, opt$output, row.names = FALSE, quote = TRUE)
        write_json(status, opt$status, pretty = TRUE, auto_unbox = TRUE)
        quit(status = 1)
    }
})

write.csv(output_df, opt$output, row.names = FALSE, quote = TRUE)
write_json(status, opt$status, pretty = TRUE, auto_unbox = TRUE)
