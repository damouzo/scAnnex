#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggridges)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(pathview)
  library(ReactomePA)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(jsonlite)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    "--top-pathways" = 10L,
    "--pathview-top-n" = 10L,
    "--dashboard-max-pathways" = 50L
  )

  required <- c("--dge-file", "--organism", "--comparison", "--outdir")
  values <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      stop(sprintf("Missing value for argument: %s", key))
    }
    values[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  for (k in names(defaults)) {
    if (is.null(values[[k]])) {
      values[[k]] <- defaults[[k]]
    }
  }

  missing <- required[!required %in% names(values)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }

  list(
    dge_file = values[["--dge-file"]],
    organism = values[["--organism"]],
    comparison = values[["--comparison"]],
    outdir = values[["--outdir"]],
    top_pathways = as.integer(values[["--top-pathways"]]),
    pathview_top_n = as.integer(values[["--pathview-top-n"]]),
    dashboard_max_pathways = as.integer(values[["--dashboard-max-pathways"]])
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

safe_write_csv <- function(df, file_path) {
  if (is.null(df)) {
    data.table::fwrite(data.frame(), file_path)
  } else {
    data.table::fwrite(df, file_path)
  }
}

build_gene_rank <- function(dge) {
  if (!all(c("gene", "log2fc") %in% names(dge))) {
    stop("DGE file must include 'gene' and 'log2fc' columns")
  }

  rank_df <- dge %>%
    dplyr::select(gene, log2fc) %>%
    dplyr::mutate(gene = as.character(gene), log2fc = as.numeric(log2fc)) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2fc), is.finite(log2fc)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(log2fc = mean(log2fc), .groups = "drop")

  if (nrow(rank_df) < 20) {
    stop("Insufficient ranked genes after cleanup (<20)")
  }

  rank_df$gene
  gene_rank <- rank_df$log2fc
  names(gene_rank) <- rank_df$gene
  sort(gene_rank, decreasing = TRUE)
}

map_symbols_to_entrez <- function(symbols, orgdb) {
  mapping <- tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
  }, error = function(e) {
    data.frame()
  })

  if (nrow(mapping) == 0) {
    return(mapping)
  }

  mapping %>%
    dplyr::filter(!is.na(SYMBOL), !is.na(ENTREZID)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

safe_enrich_plot <- function(res_obj, file_base, top_n, gene_rank = NULL, orig_obj = NULL) {
  # orig_obj: gseaResult before setReadable; ridgeplot requires @geneList and
  # core_enrichment IDs to be consistent, which setReadable breaks.
  out <- list(dotplot = NULL, ridgeplot = NULL)
  if (is.null(res_obj)) {
    return(out)
  }

  res_df <- as.data.frame(res_obj)
  if (nrow(res_df) == 0) {
    return(out)
  }

  n_show <- min(top_n, nrow(res_df))

  dot <- enrichplot::dotplot(res_obj, showCategory = n_show) +
    ggplot2::ggtitle(sprintf("Dotplot (%d pathways)", n_show))
  ggplot2::ggsave(sprintf("%s_dotplot.pdf", file_base), plot = dot, width = 10, height = 8)
  out$dotplot <- sprintf("%s_dotplot.pdf", file_base)

  # Use original object for ridgeplot to avoid geneList/core_enrichment ID mismatch
  ridge_source <- if (!is.null(orig_obj)) orig_obj else res_obj
  ridge_df <- as.data.frame(ridge_source)
  ridge_ok <- "core_enrichment" %in% colnames(ridge_df) && nrow(ridge_df) > 0
  if (ridge_ok) {
    # enrichplot::ridgeplot checks `inherits(showCategory, "numeric")`;
    # integer values fail that check, so cast to numeric explicitly.
    n_ridge <- as.numeric(min(top_n, nrow(ridge_df)))

    tryCatch({
      ridge <- enrichplot::ridgeplot(
        ridge_source,
        showCategory = n_ridge,
        fill = "p.adjust",
        core_enrichment = TRUE,
        label_format = 30
      ) + ggplot2::ggtitle(sprintf("Ridgeplot (%d pathways)", as.integer(n_ridge)))

      ggplot2::ggsave(sprintf("%s_ridgeplot.pdf", file_base), plot = ridge, width = 11, height = 8)
      out$ridgeplot <- sprintf("%s_ridgeplot.pdf", file_base)
    }, error = function(e) {
      message(sprintf("Ridgeplot skipped: %s", conditionMessage(e)))
    })
  }

  if (!is.null(gene_rank) && nrow(res_df) > 0) {
    idx <- seq_len(n_show)
    tryCatch({
      gsea_p <- enrichplot::gseaplot2(
        res_obj,
        geneSetID = idx,
        title = sprintf("Running score (%d pathways)", length(idx)),
        color = grDevices::hcl.colors(length(idx), "Dark 3")
      )
      ggplot2::ggsave(sprintf("%s_gseaplot.pdf", file_base), plot = gsea_p, width = 12, height = 9)
      out$gseaplot <- sprintf("%s_gseaplot.pdf", file_base)
    }, error = function(e) {
      message(sprintf("Gseaplot skipped: %s", conditionMessage(e)))
    })
  }

  out
}

run_gsea_set <- function(rank_values, rank_entrez, orgdb, kegg_org, reactome_org, top_n) {
  results <- list()

  go_bp <- tryCatch({
    gseGO(
      geneList = rank_entrez,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1.0,
      eps = 1e-10,
      verbose = FALSE
    )
  }, error = function(e) NULL)
  results$go_bp <- go_bp

  kegg <- tryCatch({
    gseKEGG(
      geneList = rank_entrez,
      organism = kegg_org,
      pAdjustMethod = "BH",
      pvalueCutoff = 1.0,
      eps = 1e-10,
      verbose = FALSE
    )
  }, error = function(e) NULL)
  results$kegg <- kegg

  reactome <- tryCatch({
    gsePathway(
      geneList = rank_entrez,
      organism = reactome_org,
      pAdjustMethod = "BH",
      pvalueCutoff = 1.0,
      eps = 1e-10,
      verbose = FALSE
    )
  }, error = function(e) NULL)
  results$reactome <- reactome

  results
}

main <- function() {
  args <- parse_args()
  ensure_dir(args$outdir)
  ensure_dir(file.path(args$outdir, "pathview"))

  if (!(args$organism %in% c("human", "mouse"))) {
    stop("organism must be 'human' or 'mouse'")
  }

  if (args$organism == "human") {
    orgdb <- org.Hs.eg.db
    kegg_org <- "hsa"
    reactome_org <- "human"
  } else {
    orgdb <- org.Mm.eg.db
    kegg_org <- "mmu"
    reactome_org <- "mouse"
  }

  dge <- data.table::fread(args$dge_file, data.table = FALSE)
  if (!all(c("gene", "log2fc") %in% names(dge))) {
    stop("Input DGE CSV requires columns: gene, log2fc")
  }

  gene_rank_symbol <- build_gene_rank(dge)
  symbols <- names(gene_rank_symbol)
  mapping <- map_symbols_to_entrez(symbols, orgdb)

  if (nrow(mapping) == 0) {
    stop("No SYMBOL->ENTREZ mappings found. Cannot run GSEA")
  }

  map_idx <- match(mapping$SYMBOL, names(gene_rank_symbol))
  rank_entrez <- gene_rank_symbol[map_idx]
  names(rank_entrez) <- mapping$ENTREZID
  rank_entrez <- sort(rank_entrez, decreasing = TRUE)

  gsea_results <- run_gsea_set(
    rank_values = gene_rank_symbol,
    rank_entrez = rank_entrez,
    orgdb = orgdb,
    kegg_org = kegg_org,
    reactome_org = reactome_org,
    top_n = args$top_pathways
  )

  output_tables <- list(
    GSEA_GO_BP = NULL,
    GSEA_KEGG = NULL,
    GSEA_REACTOME = NULL
  )

  output_plots <- list()

  # Keep pre-setReadable objects for dashboard ridgeplot compatibility
  gsea_results_orig <- list(go_bp = NULL, kegg = NULL, reactome = NULL)

  if (!is.null(gsea_results$go_bp) && nrow(as.data.frame(gsea_results$go_bp)) > 0) {
    gsea_results_orig$go_bp <- gsea_results$go_bp
    go_obj <- clusterProfiler::setReadable(gsea_results$go_bp, OrgDb = orgdb, keyType = "ENTREZID")
    output_tables$GSEA_GO_BP <- as.data.frame(go_obj)
    p <- safe_enrich_plot(go_obj, file.path(args$outdir, "GSEA_GO_BP"), args$top_pathways, rank_entrez,
                          orig_obj = gsea_results$go_bp)
    output_plots$go_bp <- p
    gsea_results$go_bp <- go_obj
  }

  if (!is.null(gsea_results$kegg) && nrow(as.data.frame(gsea_results$kegg)) > 0) {
    gsea_results_orig$kegg <- gsea_results$kegg
    kegg_obj <- clusterProfiler::setReadable(gsea_results$kegg, OrgDb = orgdb, keyType = "ENTREZID")
    output_tables$GSEA_KEGG <- as.data.frame(kegg_obj)
    p <- safe_enrich_plot(kegg_obj, file.path(args$outdir, "GSEA_KEGG"), args$top_pathways, rank_entrez,
                          orig_obj = gsea_results$kegg)
    output_plots$kegg <- p
    gsea_results$kegg <- kegg_obj
  }

  if (!is.null(gsea_results$reactome) && nrow(as.data.frame(gsea_results$reactome)) > 0) {
    gsea_results_orig$reactome <- gsea_results$reactome
    reactome_obj <- gsea_results$reactome
    output_tables$GSEA_REACTOME <- as.data.frame(reactome_obj)
    # Reactome does not go through setReadable, orig_obj == res_obj
    p <- safe_enrich_plot(reactome_obj, file.path(args$outdir, "GSEA_REACTOME"), args$top_pathways, rank_entrez)
    output_plots$reactome <- p
    gsea_results$reactome <- reactome_obj
  }

  for (name in names(output_tables)) {
    safe_write_csv(output_tables[[name]], file.path(args$outdir, sprintf("%s.csv", name)))
  }

  if (!is.null(output_tables$GSEA_KEGG) && nrow(output_tables$GSEA_KEGG) > 0) {
    top_kegg <- head(output_tables$GSEA_KEGG, args$pathview_top_n)
    for (i in seq_len(nrow(top_kegg))) {
      pid <- as.character(top_kegg$ID[i])
      pname <- as.character(top_kegg$Description[i])
      cat(sprintf("Running Pathview: [%s] %s\n", pid, pname))
      pathview::pathview(
        gene.data = rank_entrez,
        pathway.id = pid,
        species = kegg_org,
        out.suffix = args$comparison,
        kegg.dir = file.path(args$outdir, "pathview")
      )
    }

    generated <- list.files(pattern = paste0("^", kegg_org, ".*", args$comparison, ".*\\.png$"), full.names = TRUE)
    if (length(generated) > 0) {
      file.copy(generated, file.path(args$outdir, "pathview", basename(generated)), overwrite = TRUE)
      file.remove(generated)
    }
  }

  dashboard_data <- list(
    comparison = args$comparison,
    organism = args$organism,
    top_pathways_default = args$top_pathways,
    max_pathways_dashboard = args$dashboard_max_pathways,
    tables = output_tables,
    gsea_results = gsea_results,
    gsea_results_orig = gsea_results_orig
  )
  saveRDS(dashboard_data, file.path(args$outdir, "gsea_dashboard_data.rds"))

  summary <- list(
    comparison = args$comparison,
    organism = args$organism,
    n_ranked_genes_symbol = length(gene_rank_symbol),
    n_ranked_genes_entrez = length(rank_entrez),
    mapping_rate = round(length(rank_entrez) / length(gene_rank_symbol), 4),
    n_go_bp = if (!is.null(output_tables$GSEA_GO_BP)) nrow(output_tables$GSEA_GO_BP) else 0,
    n_kegg = if (!is.null(output_tables$GSEA_KEGG)) nrow(output_tables$GSEA_KEGG) else 0,
    n_reactome = if (!is.null(output_tables$GSEA_REACTOME)) nrow(output_tables$GSEA_REACTOME) else 0,
    generated_at = as.character(Sys.time())
  )
  writeLines(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE), file.path(args$outdir, "gsea_summary.json"))
}

main()
