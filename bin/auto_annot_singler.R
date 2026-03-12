#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleR)
    library(celldex)
    library(SingleCellExperiment)
})

print_help <- function() {
    cat("Usage: Rscript auto_annot_singler.R --input <in.rds> --output <out.csv> --status <status.json> [--refs \"BlueprintEncodeData\"] [--prune] [--continue-on-error]\n")
}

parse_args_simple <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    opt <- list(
        input = NULL,
        output = NULL,
        status = NULL,
        refs = "BlueprintEncodeData",
        prune = FALSE,
        `continue-on-error` = FALSE
    )

    i <- 1L
    while (i <= length(argv)) {
        arg <- argv[[i]]
        if (arg %in% c("-h", "--help")) {
            print_help()
            quit(status = 0)
        } else if (arg %in% c("--input", "--output", "--status", "--refs")) {
            i <- i + 1L
            if (i > length(argv)) stop(sprintf("Missing value for %s", arg), call. = FALSE)
            key <- sub("^--", "", arg)
            opt[[key]] <- argv[[i]]
        } else if (arg == "--prune") {
            opt$prune <- TRUE
        } else if (arg == "--continue-on-error") {
            opt$`continue-on-error` <- TRUE
        } else {
            stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
        }
        i <- i + 1L
    }

    if (is.null(opt$input) || is.null(opt$output) || is.null(opt$status)) {
        print_help()
        stop("Missing required arguments", call. = FALSE)
    }

    opt
}

opt <- parse_args_simple()

slugify <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

load_ref <- function(ref_name) {
    if (ref_name == "BlueprintEncodeData") {
        return(celldex::BlueprintEncodeData())
    }
    if (ref_name == "HumanPrimaryCellAtlasData") {
        return(celldex::HumanPrimaryCellAtlasData())
    }
    if (ref_name == "MonacoImmuneData") {
        return(celldex::MonacoImmuneData())
    }
    stop(sprintf("Unsupported reference: %s", ref_name))
}

status <- list(tool = "singler", success = TRUE, refs = list(), errors = list())
results <- NULL
all_cells <- NULL

write_status <- function(status_obj, path) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
        jsonlite::write_json(status_obj, path, pretty = TRUE, auto_unbox = TRUE)
        return(invisible(NULL))
    }

    success_val <- ifelse(isTRUE(status_obj$success), "true", "false")
    msg <- ""
    if (length(status_obj$errors) > 0) {
        msg <- as.character(status_obj$errors[[1]])
    }
    msg <- gsub("\\", "\\\\", msg, fixed = TRUE)
    msg <- gsub("\"", "\\\"", msg, fixed = TRUE)
    payload <- sprintf('{"tool":"singler","success":%s,"message":"%s"}', success_val, msg)
    writeLines(payload, con = path, useBytes = TRUE)
}

obj <- readRDS(opt$input)
all_cells <- colnames(obj)
sce <- as.SingleCellExperiment(obj)
refs <- trimws(strsplit(opt$refs, ",")[[1]])
refs <- refs[refs != ""]

for (ref_name in refs) {
    ref_slug <- slugify(ref_name)
    tryCatch({
        ref <- load_ref(ref_name)
        pred <- SingleR(
            test = sce,
            ref = ref,
            labels = ref$label.main,
            prune = opt$prune
        )

        label_col <- paste0("auto_annot_singler_", ref_slug)
        score_col <- paste0(label_col, "_score")
        delta_col <- paste0(label_col, "_delta_next")

        pred_df <- as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE)

        score_vals <- rep(NA_real_, nrow(pred))
        if (!is.null(pred$scores)) {
            score_vals <- apply(pred$scores, 1, max)
        }

        labels_vals <- rep(NA_character_, nrow(pred_df))
        if ("labels" %in% colnames(pred_df)) {
            labels_vals <- as.character(pred_df[["labels"]])
        }

        delta_vals <- rep(NA_real_, nrow(pred_df))
        if ("delta.next" %in% colnames(pred_df)) {
            delta_vals <- as.numeric(pred_df[["delta.next"]])
        }

        block <- data.frame(
            cell_id = rownames(pred),
            stringsAsFactors = FALSE
        )
        block[[label_col]] <- labels_vals
        block[[score_col]] <- as.numeric(score_vals)
        block[[delta_col]] <- delta_vals

        if (is.null(results)) {
            results <- block
        } else {
            results <- merge(results, block, by = "cell_id", all = TRUE)
        }

        status$refs[[length(status$refs) + 1]] <- list(ref = ref_name, success = TRUE)
    }, error = function(e) {
        status$refs[[length(status$refs) + 1]] <- list(ref = ref_name, success = FALSE)
        status$errors[[length(status$errors) + 1]] <- paste0(ref_name, ": ", e$message)
        status$success <- FALSE

        label_col <- paste0("auto_annot_singler_", ref_slug)
        score_col <- paste0(label_col, "_score")
        delta_col <- paste0(label_col, "_delta_next")
        placeholder <- data.frame(cell_id = all_cells, stringsAsFactors = FALSE)
        placeholder[[label_col]] <- NA_character_
        placeholder[[score_col]] <- NA_real_
        placeholder[[delta_col]] <- NA_real_

        if (is.null(results)) {
            results <<- placeholder
        } else {
            results <<- merge(results, placeholder, by = "cell_id", all = TRUE)
        }

        if (!opt$`continue-on-error`) {
            stop(e)
        }
    })
}

if (is.null(results)) {
    results <- data.frame(cell_id = all_cells, stringsAsFactors = FALSE)
    status$success <- FALSE
}

write.csv(results, opt$output, row.names = FALSE, quote = TRUE)
write_status(status, opt$status)

if (!status$success && !opt$`continue-on-error`) {
    quit(status = 1)
}
