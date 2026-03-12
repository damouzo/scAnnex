#!/usr/bin/env Rscript

# In read-only containers, SeuratData needs a writable location to install
# reference data packages. Set this up before loading any libraries.
local_lib <- file.path(Sys.getenv("HOME", tempdir()), ".scannex_R_libs")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
    library(jsonlite)
    library(Seurat)
    library(Azimuth)
    library(SeuratData)
})


print_help <- function() {
    cat("Usage: Rscript auto_annot_azimuth.R --input <in.rds> --output <out.csv> --status <status.json> [--refs \"pbmcref\"] [--continue-on-error]\n")
}


parse_args_simple <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    opt <- list(
        input = NULL,
        output = NULL,
        status = NULL,
        refs = "pbmcref",
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

normalize_ref_code <- function(ref_name) {
    key <- trimws(tolower(ref_name))
    key <- gsub("\\s+", "", key)

    aliases <- list(
        "pbmcref" = "pbmcref",
        "human-pbmc" = "pbmcref",
        "human_pbmc" = "pbmcref",
        "humanpbmc" = "pbmcref",
        "bonemarrowref" = "bonemarrowref",
        "human-bonemarrow" = "bonemarrowref",
        "human_bonemarrow" = "bonemarrowref",
        "humanbonemarrow" = "bonemarrowref"
    )

    if (!is.null(aliases[[key]])) {
        return(aliases[[key]])
    }
    return(trimws(ref_name))
}

slugify <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

pick_col <- function(meta, candidates = character(0), regex = NULL, exclude_regex = NULL) {
    cols <- colnames(meta)

    for (name in candidates) {
        if (name %in% cols) {
            return(list(values = meta[[name]], column = name))
        }
    }

    if (!is.null(regex)) {
        hits <- grep(regex, cols, value = TRUE, ignore.case = TRUE)
        if (!is.null(exclude_regex) && length(hits) > 0) {
            hits <- hits[!grepl(exclude_regex, hits, ignore.case = TRUE)]
        }
        if (length(hits) > 0) {
            return(list(values = meta[[hits[[1]]]], column = hits[[1]]))
        }
    }

    return(list(values = rep(NA, nrow(meta)), column = NA_character_))
}

write_status <- function(status_obj, path) {
    write_json(status_obj, path, pretty = TRUE, auto_unbox = TRUE)
}

opt <- parse_args_simple()

status <- list(tool = "azimuth", success = TRUE, refs = list(), errors = list())
results <- NULL

obj <- readRDS(opt$input)
refs <- trimws(strsplit(opt$refs, ",")[[1]])
refs <- refs[refs != ""]
if (length(refs) == 0) {
    refs <- c("pbmcref")
    status$errors[[length(status$errors) + 1]] <- "Empty --refs detected; using default reference 'pbmcref'"
}

if (!("RNA" %in% Assays(obj))) {
    stop("Seurat object must contain RNA assay for Azimuth", call. = FALSE)
}

query_gene_count <- nrow(obj[["RNA"]])

for (ref in refs) {
    ref_id <- normalize_ref_code(ref)
    ref_slug <- slugify(ref_id)

    tryCatch({
        # Ensure the reference data package is installed to a writable location.
        # SeuratData::InstallData() is safe to call even if already installed.
        pkg_name <- paste0(ref_id, ".SeuratData")
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
            message("Installing Azimuth reference: ", ref_id)
            SeuratData::InstallData(ref_id)
        }

        # Strip reductions, graphs and extra assays before mapping.
        # Integrated objects carry Harmony/UMAP reductions that cause
        # dimension mismatches during Azimuth's internal projection step.
        obj_query <- DietSeurat(obj, assays = "RNA", dimreducs = character(0), graphs = character(0))

        annotated <- RunAzimuth(obj_query, reference = ref_id)
        meta <- annotated@meta.data

        label_col <- paste0("auto_annot_azimuth_", ref_slug)
        score_col <- paste0(label_col, "_score")
        l1_col <- paste0(label_col, "_l1")
        l2_col <- paste0(label_col, "_l2")

        l2_detect <- pick_col(
            meta,
            candidates = c("predicted.celltype.l2", "predicted.celltype_l2"),
            regex = "predicted.*celltype.*l2"
        )
        l1_detect <- pick_col(
            meta,
            candidates = c("predicted.celltype.l1", "predicted.celltype_l1"),
            regex = "predicted.*celltype.*l1"
        )
        label_detect <- pick_col(
            meta,
            candidates = c("predicted.celltype", "predicted_celltype"),
            regex = "predicted.*celltype",
            exclude_regex = "score|prob"
        )
        score_detect <- pick_col(
            meta,
            candidates = c(
                "predicted.celltype.l2.score",
                "predicted.celltype_l2_score",
                "predicted.celltype.score",
                "predicted_celltype_score",
                "prediction.score.max",
                "predicted.score.max",
                "mapping.score"
            ),
            regex = "(predicted|prediction|mapping).*(score|max)"
        )

        l2_vals <- l2_detect$values
        l1_vals <- l1_detect$values
        label_vals <- label_detect$values
        score_vals <- suppressWarnings(as.numeric(score_detect$values))

        final_label <- as.character(l2_vals)
        missing_label <- is.na(final_label) | final_label == ""
        final_label[missing_label] <- as.character(l1_vals[missing_label])
        missing_label <- is.na(final_label) | final_label == ""
        final_label[missing_label] <- as.character(label_vals[missing_label])

        block <- data.frame(
            cell_id = rownames(meta),
            stringsAsFactors = FALSE
        )
        block[[label_col]] <- final_label
        block[[score_col]] <- score_vals
        block[[l1_col]] <- as.character(l1_vals)
        block[[l2_col]] <- as.character(l2_vals)

        no_labels <- all(is.na(block[[label_col]]) | block[[label_col]] == "")
        no_scores <- all(is.na(block[[score_col]]))
        if (no_labels && no_scores) {
            status$refs[[length(status$refs) + 1]] <- list(
                ref = ref,
                mapped_ref = ref_id,
                success = FALSE,
                query_genes = query_gene_count,
                overlap_genes = NA,
                detected_label_col = label_detect$column,
                detected_label_l1_col = l1_detect$column,
                detected_label_l2_col = l2_detect$column,
                detected_score_col = score_detect$column,
                metadata_n_cols = ncol(meta)
            )
            status$errors[[length(status$errors) + 1]] <- paste0(
                ref,
                ": Azimuth returned no label/score columns. Available metadata columns: ",
                paste(colnames(meta), collapse = ",")
            )
            status$success <- FALSE

            if (is.null(results)) {
                results <- block
            } else {
                results <- merge(results, block, by = "cell_id", all = TRUE)
            }

            if (!opt$`continue-on-error`) {
                stop("Azimuth reference produced empty annotations")
            }

            next
        }

        if (is.null(results)) {
            results <- block
        } else {
            results <- merge(results, block, by = "cell_id", all = TRUE)
        }

        status$refs[[length(status$refs) + 1]] <- list(
            ref = ref,
            mapped_ref = ref_id,
            success = TRUE,
            query_genes = query_gene_count,
            overlap_genes = NA,
            detected_label_col = label_detect$column,
            detected_label_l1_col = l1_detect$column,
            detected_label_l2_col = l2_detect$column,
            detected_score_col = score_detect$column,
            metadata_n_cols = ncol(meta)
        )
    }, error = function(e) {
        msg <- paste0(ref, ": ", conditionMessage(e))
        message("ERROR running Azimuth for ref '", ref, "': ", conditionMessage(e))

        # Use <<- to update variables in the enclosing (script) scope.
        # Plain <- inside an error handler creates local copies only.
        new_refs <- status$refs
        new_refs[[length(new_refs) + 1]] <- list(ref = ref, mapped_ref = ref_id, success = FALSE)
        status$refs <<- new_refs

        new_errors <- status$errors
        new_errors[[length(new_errors) + 1]] <- msg
        status$errors <<- new_errors

        status$success <<- FALSE

        placeholder <- data.frame(cell_id = colnames(obj), stringsAsFactors = FALSE)
        label_col <- paste0("auto_annot_azimuth_", ref_slug)
        score_col <- paste0(label_col, "_score")
        l1_col <- paste0(label_col, "_l1")
        l2_col <- paste0(label_col, "_l2")
        placeholder[[label_col]] <- NA_character_
        placeholder[[score_col]] <- NA_real_
        placeholder[[l1_col]] <- NA_character_
        placeholder[[l2_col]] <- NA_character_

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
    results <- data.frame(cell_id = colnames(obj), stringsAsFactors = FALSE)
    status$success <- FALSE
}

write.csv(results, opt$output, row.names = FALSE, quote = TRUE)
write_status(status, opt$status)

if (!status$success && !opt$`continue-on-error`) {
    quit(status = 1)
}
