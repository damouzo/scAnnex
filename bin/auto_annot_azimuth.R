#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(jsonlite)
    library(Seurat)
    library(Azimuth)
})

print_help <- function() {
    cat("Usage: Rscript auto_annot_azimuth.R --input <in.rds> --output <out.csv> --status <status.json> [--refs \"Human - PBMC\"] [--continue-on-error]\n")
}

parse_args_simple <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    opt <- list(
        input = NULL,
        output = NULL,
        status = NULL,
        refs = "Human - PBMC",
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

map_ref_name <- function(ref_name) {
    key <- trimws(tolower(ref_name))
    mapping <- list(
        "human - pbmc" = "pbmcref",
        "human - bone marrow" = "bonemarrowref"
    )
    if (!is.null(mapping[[key]])) {
        return(mapping[[key]])
    }
    return(ref_name)
}

slugify <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

status <- list(tool = "azimuth", success = TRUE, refs = list(), errors = list())
results <- NULL

obj <- readRDS(opt$input)
refs <- trimws(strsplit(opt$refs, ",")[[1]])
refs <- refs[refs != ""]

for (ref in refs) {
    ref_id <- map_ref_name(ref)
    ref_slug <- slugify(ref)

    tryCatch({
        annotated <- RunAzimuth(obj, reference = ref_id)
        meta <- annotated@meta.data

        label_col <- paste0("auto_annot_azimuth_", ref_slug)
        score_col <- paste0(label_col, "_score")
        l1_col <- paste0(label_col, "_l1")
        l2_col <- paste0(label_col, "_l2")

        block <- data.frame(
            cell_id = rownames(meta),
            stringsAsFactors = FALSE
        )

        if ("predicted.celltype.l2" %in% colnames(meta)) {
            block[[label_col]] <- as.character(meta$predicted.celltype.l2)
        } else if ("predicted.celltype.l1" %in% colnames(meta)) {
            block[[label_col]] <- as.character(meta$predicted.celltype.l1)
        } else {
            block[[label_col]] <- NA_character_
        }

        if ("prediction.score.max" %in% colnames(meta)) {
            block[[score_col]] <- as.numeric(meta$prediction.score.max)
        } else {
            block[[score_col]] <- NA_real_
        }

        if ("predicted.celltype.l1" %in% colnames(meta)) {
            block[[l1_col]] <- as.character(meta$predicted.celltype.l1)
        } else {
            block[[l1_col]] <- NA_character_
        }

        if ("predicted.celltype.l2" %in% colnames(meta)) {
            block[[l2_col]] <- as.character(meta$predicted.celltype.l2)
        } else {
            block[[l2_col]] <- NA_character_
        }

        if (is.null(results)) {
            results <- block
        } else {
            results <- merge(results, block, by = "cell_id", all = TRUE)
        }

        status$refs[[length(status$refs) + 1]] <- list(ref = ref, mapped_ref = ref_id, success = TRUE)
    }, error = function(e) {
        status$refs[[length(status$refs) + 1]] <- list(ref = ref, mapped_ref = ref_id, success = FALSE)
        status$errors[[length(status$errors) + 1]] <- paste0(ref, ": ", e$message)
        if (!opt$`continue-on-error`) {
            stop(e)
        }
    })
}

if (is.null(results)) {
    results <- data.frame(cell_id = character(0), stringsAsFactors = FALSE)
    status$success <- FALSE
}

write.csv(results, opt$output, row.names = FALSE, quote = TRUE)
write_json(status, opt$status, pretty = TRUE, auto_unbox = TRUE)

if (!status$success && !opt$`continue-on-error`) {
    quit(status = 1)
}
