# scAnnex Dashboard - Server Logic
# Handles reactive data processing and plot generation

server <- function(input, output, session) {
  
  # ===========================================================================
  # REACTIVE VALUES
  # ===========================================================================
  
  rv <- reactiveValues(
    data_obj        = NULL,
    qc_report       = NULL,
    qc_reports_multi = list(),
    qc_plots        = list(),
    qc_cells_data   = NULL,
    data_loaded     = FALSE,
    h5ad_path       = DEFAULT_MERGED_H5AD,
    results_dir     = RESULTS_DIR,
    umap_color_choices   = character(0),
    umap_color_label_map = list()
  )

  sanitize_color_name <- function(x) {
    out <- x
    out <- gsub("^auto_annot_", "", out)
    out <- gsub("^celltypist_", "celltypist ", out)
    out <- gsub("^singler_", "singler ", out)
    out <- gsub("^azimuth_", "azimuth ", out)
    out <- gsub("^sctype", "sctype", out)
    out <- gsub("_score$", " score", out)
    out <- gsub("_delta_next$", " delta", out)
    out <- gsub("_pruned$", " pruned", out)
    out <- gsub("_l1$", " l1", out)
    out <- gsub("_l2$", " l2", out)
    out <- gsub("_", " ", out)
    out <- gsub("\\bpkl\\b", "PKL", out, ignore.case = TRUE)
    out <- trimws(out)
    ifelse(out == "", x, out)
  }

  # =========================================================================
  # AUTO-LOAD on startup (replaces manual Load Data button)
  # =========================================================================
  observeEvent(TRUE, {
    withProgress(message = "Initializing scAnnex dashboard...", value = 0, {

      # -- 1. QC reports ----------------------------------------------------
      tryCatch({
        incProgress(0.1, detail = "Loading QC reports")
        rv$qc_reports_multi <- load_all_qc_reports(RESULTS_DIR)
        if (length(rv$qc_reports_multi) > 0) {
          sample_choices <- names(rv$qc_reports_multi)
          rv$qc_report   <- rv$qc_reports_multi[[sample_choices[1]]]
          updateSelectInput(session, "qc_sample_select",
                            choices = sample_choices, selected = sample_choices[1])
        }
      }, error = function(e) message("QC reports: ", e$message))

      # -- 2. QC cell data for density plots --------------------------------
      tryCatch({
        incProgress(0.05, detail = "Loading per-cell QC data")
        rv$qc_cells_data <- load_qc_cells_data(RESULTS_DIR)
      }, error = function(e) message("QC cells data: ", e$message))

      # -- 3. H5AD (main bottleneck) ----------------------------------------
      if (nzchar(DEFAULT_MERGED_H5AD) && file.exists(DEFAULT_MERGED_H5AD)) {
        tryCatch({
          incProgress(0.4, detail = "Loading H5AD (this may take a moment)")
          rv$data_obj <- load_h5ad_data(DEFAULT_MERGED_H5AD, backed = FALSE)
          rv$h5ad_path <- DEFAULT_MERGED_H5AD

          rv$umap_color_choices <- setdiff(names(rv$data_obj$metadata), "cell_id")
          label_choices <- vapply(rv$umap_color_choices, sanitize_color_name, character(1))
          rv$umap_color_label_map <- as.list(label_choices)
          names(rv$umap_color_label_map) <- rv$umap_color_choices
          display_choices <- setNames(rv$umap_color_choices, label_choices)

          updateSelectInput(session, "umap_color_by",
                            choices  = display_choices,
                            selected = if ("batch" %in% rv$umap_color_choices) "batch"
                                       else rv$umap_color_choices[1])
          rv$data_loaded <- TRUE
        }, error = function(e) message("H5AD load: ", e$message))
      }

      # -- 4. DGE -----------------------------------------------------------
      if (nzchar(DEFAULT_DGE_DIR) && dir.exists(DEFAULT_DGE_DIR)) {
        tryCatch({
          incProgress(0.15, detail = "Loading DGE results")
          load_dge_results(DEFAULT_DGE_DIR, show_progress = FALSE)
        }, error = function(e) message("DGE: ", e$message))
      }

      # -- 5. GSEA ----------------------------------------------------------
      if (nzchar(DEFAULT_GSEA_DIR) && dir.exists(DEFAULT_GSEA_DIR)) {
        tryCatch({
          incProgress(0.25, detail = "Loading GSEA results")
          load_gsea_results(DEFAULT_GSEA_DIR, show_progress = FALSE)
        }, error = function(e) message("GSEA: ", e$message))
      }

      incProgress(0.05, detail = "Ready")
    })
  }, once = TRUE, ignoreNULL = FALSE, ignoreInit = FALSE)

  # is_integrated: TRUE when multi-sample QC reports are present
  output$is_integrated <- reactive({
    length(rv$qc_reports_multi) > 1
  })
  outputOptions(output, "is_integrated", suspendWhenHidden = FALSE)

  # Overview: results directory label
  output$overview_results_dir_label <- renderText({
    if (nzchar(rv$results_dir) && dir.exists(rv$results_dir)) {
      paste0("Results: ", rv$results_dir)
    } else {
      "Results directory not found — set SCANNEX_RESULTS_DIR environment variable"
    }
  })

  # Overview: KPI card values
  output$ov_n_samples <- renderText({
    n <- length(rv$qc_reports_multi)
    if (n == 0) "—" else as.character(n)
  })

  output$ov_cells_before <- renderText({
    reps <- rv$qc_reports_multi
    if (length(reps) == 0) return("—")
    total <- sum(vapply(reps, function(r) {
      as.integer(r$filtering_statistics$cells_initial %||% 0)
    }, integer(1)))
    format_number(total)
  })

  output$ov_cells_after <- renderText({
    reps <- rv$qc_reports_multi
    if (length(reps) == 0) return("—")
    total <- sum(vapply(reps, function(r) {
      as.integer(r$filtering_statistics$cells_final %||% 0)
    }, integer(1)))
    format_number(total)
  })

  output$ov_avg_retention <- renderText({
    reps <- rv$qc_reports_multi
    if (length(reps) == 0) return("—")
    vals <- vapply(reps, function(r) {
      as.numeric(r$filtering_statistics$cells_retained_pct %||% NA_real_)
    }, numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return("—")
    sprintf("%.1f%%", mean(vals))
  })

  output$ov_avg_genes <- renderText({
    reps <- rv$qc_reports_multi
    if (length(reps) == 0) return("—")
    vals <- vapply(reps, function(r) {
      as.numeric(r$qc_metrics_after$n_genes_by_counts$median %||% NA_real_)
    }, numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return("—")
    format_number(round(mean(vals)))
  })

  output$ov_avg_mt <- renderText({
    reps <- rv$qc_reports_multi
    if (length(reps) == 0) return("—")
    vals <- vapply(reps, function(r) {
      as.numeric(r$qc_metrics_after$pct_counts_mt$median %||% NA_real_)
    }, numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return("—")
    sprintf("%.1f%%", mean(vals))
  })

  # Overview: per-sample summary table
  output$overview_sample_table <- renderDT({
    reps <- rv$qc_reports_multi
    if (length(reps) == 0) {
      return(datatable(data.frame(Message = "QC reports not loaded"),
                       options = list(dom = "t"), rownames = FALSE))
    }
    rows <- lapply(names(reps), function(sid) {
      r  <- reps[[sid]]
      fs <- r$filtering_statistics
      ma <- r$qc_metrics_after
      data.frame(
        Sample            = sid,
        Cells_Before      = as.integer(fs$cells_initial    %||% NA),
        Cells_After       = as.integer(fs$cells_final      %||% NA),
        Retention_pct     = round(as.numeric(fs$cells_retained_pct %||% NA), 1),
        Genes_After       = as.integer(fs$genes_final      %||% NA),
        Median_Genes      = round(as.numeric(ma$n_genes_by_counts$median %||% NA), 0),
        Median_Counts     = round(as.numeric(ma$total_counts$median      %||% NA), 0),
        Median_MT_pct     = round(as.numeric(ma$pct_counts_mt$median     %||% NA), 2),
        stringsAsFactors  = FALSE, check.names = FALSE
      )
    })
    df <- do.call(rbind, rows)
    datatable(df,
      options = list(pageLength = 20, scrollX = TRUE, dom = "frtip"),
      rownames = FALSE
    ) %>%
    formatStyle("Cells_Before",
      background = styleColorBar(c(0, max(df$Cells_Before, na.rm = TRUE)), "#d4e6f1"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    ) %>%
    formatStyle("Cells_After",
      background = styleColorBar(c(0, max(df$Cells_After, na.rm = TRUE)), "#a9dfbf"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    ) %>%
    formatStyle("Retention_pct",
      background = styleColorBar(c(0, 100), "#b3d9f7"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    ) %>%
    formatStyle("Genes_After",
      background = styleColorBar(c(0, max(df$Genes_After, na.rm = TRUE)), "#d7bde2"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    ) %>%
    formatStyle("Median_Genes",
      background = styleColorBar(c(0, max(df$Median_Genes, na.rm = TRUE)), "#fad7a0"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    ) %>%
    formatStyle("Median_Counts",
      background = styleColorBar(c(0, max(df$Median_Counts, na.rm = TRUE)), "#a2d9ce"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    ) %>%
    formatStyle("Median_MT_pct",
      background = styleColorBar(c(0, max(df$Median_MT_pct, na.rm = TRUE)), "#f9e79f"),
      backgroundSize = "100% 90%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    )
  })

  # Shared color palette for density plots (derived from sample IDs in qc_cells_data)
  sa_palette <- reactive({
    req(rv$qc_cells_data)
    get_sample_palette(unique(rv$qc_cells_data$sample_id))
  })

  overview_density_plot <- function(metric_col, x_label) {
    req(rv$qc_cells_data)
    df  <- rv$qc_cells_data
    pal <- sa_palette()
    req(metric_col %in% names(df))
    ggplot(df, aes_string(x = metric_col, color = "sample_id", fill = "sample_id")) +
      geom_density(alpha = 0.12, linewidth = 0.7) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values  = pal) +
      labs(x = x_label, y = "Density", color = "Sample", fill = "Sample") +
      theme_bw(base_size = 12) +
      theme(
        legend.position  = "bottom",
        legend.text      = element_text(size = 11),
        legend.key.size  = unit(0.5, "lines"),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(nrow = 3), color = guide_legend(nrow = 3))
  }

  output$overview_density_mt <- renderPlot({
    overview_density_plot("pct_counts_mt", "Mitochondrial %")
  })
  output$overview_density_genes <- renderPlot({
    overview_density_plot("n_genes_by_counts", "Genes per cell")
  })
  output$overview_density_counts <- renderPlot({
    overview_density_plot("total_counts", "Total counts")
  })

  output$overview_density_status <- renderUI({
    if (is.null(rv$qc_cells_data)) {
      tags$p(class = "text-muted", style = "padding: 8px;",
             "Per-cell QC data not available (pipeline must have generated *_qc.h5ad files).")
    }
  })
  
  # ===========================================================================
  # DATA LOADING
  # ===========================================================================
  # (Data loading is now handled automatically by the startup observer above)

  
  # ===========================================================================
  # TAB 2: QC OVERVIEW
  # ===========================================================================

  qc_report_active <- reactive({
    # Multi-sample: use per-sample selector when more than one report is loaded
    if (length(rv$qc_reports_multi) > 0) {
      selected_sample <- input$qc_sample_select
      if (!is.null(selected_sample) && selected_sample %in% names(rv$qc_reports_multi)) {
        return(rv$qc_reports_multi[[selected_sample]])
      }
    }

    rv$qc_report
  })

  observeEvent(input$qc_sample_select, {
    req(input$qc_sample_select)
    results_dir <- rv$results_dir
    if (!nzchar(results_dir) || !dir.exists(results_dir)) return()

    if (input$qc_sample_select %in% names(rv$qc_reports_multi)) {
      rv$qc_report <- rv$qc_reports_multi[[input$qc_sample_select]]
    }

    sample_plots <- get_qc_plots_for_sample(results_dir, input$qc_sample_select)
    if (length(sample_plots) > 0) {
      rv$qc_plots <- sample_plots
    }
  }, ignoreInit = TRUE)
  
  # QC Summary value boxes
  output$qc_box_cells_before <- renderUI({
    report <- qc_report_active()
    req(report)
    n <- format_number(report$filtering_statistics$cells_initial)
    value_box(title = "Cells (Before QC)", value = n,
              showcase = icon("circle"), theme = "primary")
  })

  output$qc_box_cells_after <- renderUI({
    report <- qc_report_active()
    req(report)
    n <- format_number(report$filtering_statistics$cells_final)
    value_box(title = "Cells (After QC)", value = n,
              showcase = icon("check-circle"), theme = "success")
  })

  output$qc_box_genes_after <- renderUI({
    report <- qc_report_active()
    req(report)
    n <- format_number(report$filtering_statistics$genes_final)
    value_box(title = "Genes (After QC)", value = n,
              showcase = icon("dna"), theme = "info")
  })

  output$qc_box_retention <- renderUI({
    report <- qc_report_active()
    req(report)
    pct <- sprintf("%.1f%%", as.numeric(report$filtering_statistics$cells_retained_pct))
    value_box(title = "Cell Retention", value = pct,
              showcase = icon("percent"), theme = "warning")
  })

  output$qc_box_median_genes <- renderUI({
    report <- qc_report_active()
    req(report)
    ma <- report$qc_metrics_after
    val <- if (!is.null(ma)) format_number(round(ma$n_genes_by_counts$median, 0)) else "—"
    value_box(title = "Median Genes / Cell", value = val,
              showcase = icon("dna"), theme = "info")
  })

  output$qc_box_median_mt <- renderUI({
    report <- qc_report_active()
    req(report)
    ma <- report$qc_metrics_after
    val <- if (!is.null(ma)) sprintf("%.2f%%", ma$pct_counts_mt$median) else "—"
    value_box(title = "Median MT%", value = val,
              showcase = icon("biohazard"), theme = "secondary")
  })
  
  # QC Metrics Table
  output$qc_metrics_table <- renderDT({
    report <- qc_report_active()
    req(report)

    # Extract metrics
    metrics_before <- report$qc_metrics_before
    metrics_after <- report$qc_metrics_after
    
    # Build table
    metrics_df <- data.frame(
      Metric = c("Genes per cell", "UMI counts", "% Mitochondrial"),
      stringsAsFactors = FALSE
    )
    
    # Add Before columns
    metrics_df$Before_Median <- c(
      round(metrics_before$n_genes_by_counts$median, 0),
      round(metrics_before$total_counts$median, 0),
      round(metrics_before$pct_counts_mt$median, 2)
    )
    
    metrics_df$Before_Mean <- c(
      round(metrics_before$n_genes_by_counts$mean, 0),
      round(metrics_before$total_counts$mean, 0),
      round(metrics_before$pct_counts_mt$mean, 2)
    )
    
    # Add After columns if available
    if (!is.null(metrics_after)) {
      metrics_df$After_Median <- c(
        round(metrics_after$n_genes_by_counts$median, 0),
        round(metrics_after$total_counts$median, 0),
        round(metrics_after$pct_counts_mt$median, 2)
      )
      
      metrics_df$After_Mean <- c(
        round(metrics_after$n_genes_by_counts$mean, 0),
        round(metrics_after$total_counts$mean, 0),
        round(metrics_after$pct_counts_mt$mean, 2)
      )
    }
    
    datatable(
      metrics_df,
      options = list(
        dom = 't',
        pageLength = 10,
        ordering = FALSE
      ),
      rownames = FALSE
    )
  })
  
  # QC Thresholds Table
  output$qc_thresholds_table <- renderDT({
    report <- qc_report_active()
    req(report)
    
    thresholds <- report$thresholds_applied
    
    if (is.character(thresholds) && thresholds == "manual") {
      # Return simple message for manual thresholds
      threshold_df <- data.frame(
        Metric = "Manual thresholds",
        `Lower Bound` = "-",
        `Upper Bound` = "-",
        check.names = FALSE
      )
    } else {
      # Build table from thresholds
      metric_names <- c(
        "n_genes_by_counts" = "Genes per cell",
        "total_counts" = "UMI counts",
        "pct_counts_mt" = "% Mitochondrial",
        "pct_counts_ribo" = "% Ribosomal",
        "pct_counts_hb" = "% Hemoglobin"
      )
      
      threshold_list <- list()
      
      for (metric in names(thresholds)) {
        vals <- thresholds[[metric]]
        
        # Format lower bound
        lower <- if (!is.null(vals[1]) && !is.na(vals[1])) {
          sprintf("%.1f", vals[1])
        } else {
          "No limit"
        }
        
        # Format upper bound
        upper <- if (length(vals) > 1 && !is.null(vals[2]) && !is.na(vals[2])) {
          sprintf("%.1f", vals[2])
        } else {
          "No limit"
        }
        
        # Use friendly name or fallback to metric name
        metric_display <- if (metric %in% names(metric_names)) {
          metric_names[[metric]]
        } else {
          metric
        }
        
        threshold_list[[length(threshold_list) + 1]] <- data.frame(
          Metric = metric_display,
          `Lower Bound` = lower,
          `Upper Bound` = upper,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
      
      threshold_df <- do.call(rbind, threshold_list)
    }
    
    datatable(
      threshold_df,
      options = list(
        dom = 't',
        pageLength = 10,
        ordering = FALSE,
        columnDefs = list(
          list(className = 'dt-center', targets = c(1, 2))
        )
      ),
      rownames = FALSE
    )
  })
  
  # QC Plot Before
  output$qc_plot_before <- renderImage({
    req(rv$qc_plots)
    req(input$qc_plot_before_select)
    
    # Map selection to file pattern
    plot_type <- tolower(input$qc_plot_before_select)
    pattern <- sprintf("qc_before_%s.png", plot_type)
    
    # Find matching plot
    plot_file <- rv$qc_plots[grep(pattern, rv$qc_plots)]
    
    if (length(plot_file) == 0) {
      return(list(src = "", alt = "Plot not found"))
    }
    
    list(
      src = plot_file[1],
      contentType = "image/png",
      alt = sprintf("QC Before - %s", input$qc_plot_before_select)
    )
    
  }, deleteFile = FALSE)
  
  # QC Plot After
  output$qc_plot_after <- renderImage({
    req(rv$qc_plots)
    req(input$qc_plot_after_select)
    
    # Map selection to file pattern
    plot_type <- tolower(input$qc_plot_after_select)
    pattern <- sprintf("qc_after_%s.png", plot_type)
    
    # Find matching plot
    plot_file <- rv$qc_plots[grep(pattern, rv$qc_plots)]
    
    if (length(plot_file) == 0) {
      return(list(src = "", alt = "Plot not found"))
    }
    
    list(
      src = plot_file[1],
      contentType = "image/png",
      alt = sprintf("QC After - %s", input$qc_plot_after_select)
    )
    
  }, deleteFile = FALSE)
  
  # ===========================================================================
  # TAB 3: CLUSTERING & UMAP
  # ===========================================================================
  
  # UMAP Plot
  output$umap_plot <- renderPlotly({
    req(rv$data_loaded)
    req(rv$data_obj$umap_coords)
    req(input$umap_color_by)
    
    # Merge UMAP coords with metadata
    umap_data <- merge(
      rv$data_obj$umap_coords,
      rv$data_obj$metadata[, c("cell_id", input$umap_color_by), drop = FALSE],
      by = "cell_id"
    )

    if (!is.null(umap_data[[input$umap_color_by]])) {
      if (is.logical(umap_data[[input$umap_color_by]])) {
        umap_data[[input$umap_color_by]] <- as.character(umap_data[[input$umap_color_by]])
      }
      if (is.character(umap_data[[input$umap_color_by]])) {
        umap_data[[input$umap_color_by]] <- as.factor(umap_data[[input$umap_color_by]])
      }
      if (all(is.na(umap_data[[input$umap_color_by]]))) {
        validate(need(FALSE, sprintf("Column '%s' contains only NA values", input$umap_color_by)))
      }

      non_na_vals <- umap_data[[input$umap_color_by]][!is.na(umap_data[[input$umap_color_by]])]
      if (is.numeric(non_na_vals) && length(unique(non_na_vals)) <= 1) {
        validate(need(FALSE, sprintf("Column '%s' has no numeric range (constant value)", input$umap_color_by)))
      }
    }

    legend_title <- input$umap_color_by
    if (!is.null(rv$umap_color_label_map[[input$umap_color_by]])) {
      legend_title <- rv$umap_color_label_map[[input$umap_color_by]]
    }
    
    color_values <- umap_data[[input$umap_color_by]]
    is_numeric_color <- is.numeric(color_values)

    # Add hover label column (used in text aesthetic)
    umap_data$hover_cat <- as.character(umap_data[[input$umap_color_by]])

    if (is_numeric_color) {
      p <- plot_ly(
        data = umap_data,
        x = ~UMAP_1,
        y = ~UMAP_2,
        type = 'scattergl',
        mode = 'markers',
        marker = list(
          size = input$umap_point_size,
          opacity = input$umap_opacity,
          color = color_values,
          colorscale = 'Viridis',
          showscale = TRUE,
          colorbar = list(title = legend_title)
        ),
        text = ~paste0("Cell: ", cell_id, "<br>", legend_title, ": ", round(as.numeric(hover_cat), 3)),
        hoverinfo = 'text'
      )
    } else {
      cats <- levels(umap_data[[input$umap_color_by]])
      pal  <- setNames(SC_COLORS[((seq_along(cats) - 1L) %% length(SC_COLORS)) + 1L], cats)
      p <- plot_ly(
        data = umap_data,
        x = ~UMAP_1,
        y = ~UMAP_2,
        type = 'scattergl',
        mode = 'markers',
        marker = list(
          size = input$umap_point_size,
          opacity = input$umap_opacity
        ),
        color  = as.formula(paste0("~`", input$umap_color_by, "`")),
        colors = pal,
        text   = ~paste0("Cell: ", cell_id, "<br>", legend_title, ": ", hover_cat),
        hoverinfo = 'text'
      )
    }

    p <- p %>%
      layout(
        title  = sprintf("UMAP colored by %s", legend_title),
        xaxis  = list(title = "UMAP 1", constrain = "domain"),
        yaxis  = list(title = "UMAP 2", scaleanchor = "x", scaleratio = 1, constrain = "domain"),
        legend = list(title = list(text = legend_title),
                      x = 1.02, y = 0.5, xanchor = "left", yanchor = "middle",
                      orientation = "v", itemsizing = "constant"),
        margin    = list(r = 160, b = 40),
        hovermode = "closest"
      ) %>% config(responsive = TRUE)
    
    return(p)
  })
  
  # Metadata Table
  output$metadata_table <- renderDT({
    req(rv$data_loaded)
    
    # Reorder columns to put cell_id first
    metadata <- rv$data_obj$metadata
    
    # Get cell_id column and other columns
    if ("cell_id" %in% names(metadata)) {
      other_cols <- setdiff(names(metadata), "cell_id")
      metadata <- metadata[, c("cell_id", other_cols), drop = FALSE]
    }
    
    datatable(
      metadata,
      extensions = 'Buttons',
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        search = list(search = ''),
        dom = 'Bfrtip',
        buttons = list(
          list(
            extend = 'csv',
            text = 'Download CSV',
            filename = 'cell_metadata',
            exportOptions = list(
              modifier = list(page = "all", search = "applied")
            )
          ),
          list(
            extend = 'excel',
            text = 'Download Excel',
            filename = 'cell_metadata',
            exportOptions = list(
              modifier = list(page = "all", search = "applied")
            )
          )
        ),
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      filter = 'top',
      rownames = FALSE
    )
  })
  
  # ===========================================================================
  # TAB 4: GENE EXPRESSION
  # ===========================================================================
  
  # Gene expression UMAP - auto-detects single gene vs gene set
  output$gene_expression_umap <- renderPlotly({
    req(input$btn_plot_genes)
    
    isolate({
      req(input$gene_input)
      req(rv$data_loaded)
      req(rv$data_obj$umap_coords)
      
      # Parse input (split by newlines and clean)
      gene_list <- trimws(strsplit(input$gene_input, "\n")[[1]])
      gene_list <- gene_list[gene_list != ""]
      
      if (length(gene_list) == 0) {
        showNotification(
          "Please enter at least one gene name",
          type = "warning",
          duration = 5
        )
        return(NULL)
      }
      
      # Auto-detect: single gene vs gene set
      if (length(gene_list) == 1) {
        # ===== SINGLE GENE EXPRESSION =====
        gene_name <- gene_list[1]
        
        tryCatch({
          
          # Get gene expression
          expr <- get_gene_expression(rv$data_obj, gene_name)

          # Merge with UMAP coords
          umap_data <- rv$data_obj$umap_coords
          umap_data$expression <- expr[umap_data$cell_id]

          top_n <- as.integer(input$gene_top_n_cells %||% 0)

          if (top_n > 0) {
            # ---- Top-N highlight mode: red vs gray ----
            thresh <- sort(umap_data$expression, decreasing = TRUE)[min(top_n, nrow(umap_data))]
            umap_data$is_top <- umap_data$expression >= thresh
            # Sort so top cells render on top
            umap_data <- umap_data[order(umap_data$is_top), ]

            p <- plot_ly(
              data   = umap_data,
              x      = ~UMAP_1, y = ~UMAP_2,
              type   = "scattergl", mode = "markers",
              marker = list(size    = input$umap_point_size %||% 3,
                            opacity = input$umap_opacity    %||% 0.7,
                            color   = ifelse(umap_data$is_top, "red", "lightgray")),
              text     = ~paste0("Cell: ", cell_id, "<br>Expression: ", round(expression, 2),
                                 "<br>Top ", top_n, ": ", ifelse(is_top, "yes", "no")),
              hoverinfo = "text"
            ) %>% layout(
              title  = sprintf("Expression of %s (top %d cells highlighted)", gene_name, top_n),
              xaxis  = list(title = "UMAP 1", constrain = "domain"),
              yaxis  = list(title = "UMAP 2", scaleanchor = "x", scaleratio = 1, constrain = "domain"),
              margin = list(r = 60, b = 40),
              hovermode = "closest"
            ) %>% config(responsive = TRUE)
          } else {
            # ---- Normal gradient mode ----
            p <- plot_ly(
              data   = umap_data,
              x      = ~UMAP_1, y = ~UMAP_2,
              type   = "scattergl", mode = "markers",
              marker = list(
                size       = input$umap_point_size %||% 3,
                opacity    = input$umap_opacity    %||% 0.7,
                color      = ~expression,
                colorscale = "Viridis",
                showscale  = TRUE,
                colorbar   = list(title = "Expression")
              ),
              text     = ~paste("Cell:", cell_id, "<br>Expression:", round(expression, 2)),
              hoverinfo = "text"
            ) %>% layout(
              title  = sprintf("Expression of %s", gene_name),
              xaxis  = list(title = "UMAP 1", constrain = "domain"),
              yaxis  = list(title = "UMAP 2", scaleanchor = "x", scaleratio = 1, constrain = "domain"),
              margin = list(r = 80, b = 40),
              hovermode = "closest"
            ) %>% config(responsive = TRUE)
          }

          return(p)
          
        }, error = function(e) {
          showNotification(
            paste("Error plotting gene:", e$message),
            type = "error",
            duration = 5
          )
          return(NULL)
        })
        
      } else {
        # ===== GENE SET SCORING =====
        
        tryCatch({
          
          # Calculate gene set score
          score <- calculate_gene_set_score(rv$data_obj, gene_list)
          
          # Merge with UMAP coords
          umap_data <- rv$data_obj$umap_coords
          umap_data$score <- score[umap_data$cell_id]
          
          # Create plot
          p <- plot_ly(
            data = umap_data,
            x = ~UMAP_1,
            y = ~UMAP_2,
            type = 'scattergl',
            mode = 'markers',
            marker = list(
              size = 3,
              opacity = 0.7,
              color = ~score,
              colorscale = 'Viridis',
              showscale = TRUE,
              colorbar = list(title = "Gene Set Score")
            ),
            text = ~paste("Cell:", cell_id, "<br>Score:", round(score, 3)),
            hoverinfo = 'text'
          ) %>%
            layout(
              title  = sprintf("Gene Set Score (%d genes)", length(gene_list)),
              xaxis  = list(title = "UMAP 1", constrain = "domain"),
              yaxis  = list(title = "UMAP 2", scaleanchor = "x", scaleratio = 1, constrain = "domain"),
              margin = list(r = 80, b = 40),
              hovermode = 'closest'
            ) %>% config(responsive = TRUE)
          
          showNotification(
            sprintf("Gene set score calculated for %d genes", length(gene_list)),
            type = "message",
            duration = 3
          )
          
          return(p)
          
        }, error = function(e) {
          showNotification(
            paste("Error calculating gene set score:", e$message),
            type = "error",
            duration = 5
          )
          return(NULL)
        })
      }
    })
  })
  
  # ===========================================================================
  # TAB 5: DIFFERENTIAL EXPRESSION
  # ===========================================================================
  
  # Reactive values for DGE data
  rv_dge <- reactiveValues(
    dge_dir = NULL,
    contrasts = character(0),
    dge_results = list(),
    dge_loaded = FALSE
  )

  normalize_dge_df <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    # Harmonize column names from different DGE exporters
    if (!"log2_fc" %in% names(df) && "log2fc" %in% names(df)) {
      df$log2_fc <- df$log2fc
    }
    if (!"pvalue" %in% names(df) && "pval" %in% names(df)) {
      df$pvalue <- df$pval
    }
    if (!"pvalue_adj" %in% names(df) && "pval_adj" %in% names(df)) {
      df$pvalue_adj <- df$pval_adj
    }

    # Ensure numeric columns are numeric
    for (col_name in c("log2_fc", "pvalue", "pvalue_adj")) {
      if (col_name %in% names(df)) {
        df[[col_name]] <- suppressWarnings(as.numeric(df[[col_name]]))
      }
    }

    df
  }

  load_dge_results <- function(dge_dir, show_progress = TRUE) {
    if (!dir.exists(dge_dir)) {
      stop(sprintf("Directory not found: %s", dge_dir))
    }

    if (show_progress) {
      incProgress(0.2, detail = "Scanning for contrasts")
    }

    # Find all *_results.csv files (exclude all_contrasts_*)
    result_files <- list.files(
      dge_dir,
      pattern = "^[^all].*_results\\.csv$",
      full.names = TRUE
    )

    if (length(result_files) == 0) {
      stop("No DGE results files found (looking for *_results.csv)")
    }

    if (show_progress) {
      incProgress(0.3, detail = sprintf("Loading %d contrasts", length(result_files)))
    }

    contrast_names <- gsub("_results\\.csv$", "", basename(result_files))

    dge_data <- list()
    for (i in seq_along(result_files)) {
      contrast_name <- contrast_names[i]
      df <- read.csv(result_files[i], stringsAsFactors = FALSE)
      dge_data[[contrast_name]] <- normalize_dge_df(df)

      if (show_progress) {
        incProgress(0.4 / length(result_files))
      }
    }

    rv_dge$dge_dir <- dge_dir
    rv_dge$contrasts <- contrast_names
    rv_dge$dge_results <- dge_data
    rv_dge$dge_loaded <- TRUE

    updateSelectInput(
      session,
      "dge_contrast_select",
      choices = contrast_names,
      selected = contrast_names[1]
    )

    if (show_progress) {
      incProgress(0.1, detail = "Done!")
    }
  }
  
  # (DGE loading is now handled by the startup auto-load observer)

  # Reactive: processed DGE data for the current contrast (used by plot + hover)
  dge_plot_data <- reactive({
    req(rv_dge$dge_loaded)
    req(input$dge_contrast_select)
    dge_df <- rv_dge$dge_results[[input$dge_contrast_select]]
    if (is.null(dge_df) || nrow(dge_df) == 0) return(NULL)
    required_cols <- c("log2_fc", "pvalue_adj")
    if (!all(required_cols %in% names(dge_df))) return(NULL)
    dge_df %>%
      filter(!is.na(log2_fc), !is.na(pvalue_adj), pvalue_adj > 0) %>%
      mutate(
        neg_log10_padj = -log10(pvalue_adj),
        significant = abs(log2_fc) >= input$dge_logfc_threshold & pvalue_adj < input$dge_pval_threshold,
        direction = ifelse(
          abs(log2_fc) >= input$dge_logfc_threshold & pvalue_adj < input$dge_pval_threshold,
          ifelse(log2_fc > 0, "Upregulated", "Downregulated"),
          "Not Significant"
        )
      )
  })

  select_dge_label_genes <- function(dge_df, top_n) {
    if (is.null(dge_df) || nrow(dge_df) == 0 || is.null(top_n) || top_n <= 0) {
      return(dge_df[0, , drop = FALSE])
    }

    sig_df <- dge_df %>% filter(significant)
    if (nrow(sig_df) == 0) {
      return(sig_df)
    }

    top_n <- min(as.integer(top_n), nrow(sig_df))

    if (!"gene" %in% names(sig_df)) {
      return(
        sig_df %>%
          arrange(desc(abs(log2_fc)), pvalue_adj) %>%
          head(top_n)
      )
    }

    sig_df <- sig_df %>%
      mutate(
        .gene_key = ifelse(
          !is.na(gene) & nzchar(as.character(gene)),
          as.character(gene),
          paste0("row_", row_number())
        )
      )

    # Reserve label slots for both LFC tails so positive/negative extremes are visible.
    side_quota <- min(max(1L, floor(top_n * 0.35)), top_n)

    up_extreme <- sig_df %>%
      filter(log2_fc > 0) %>%
      arrange(desc(log2_fc), pvalue_adj) %>%
      head(side_quota)

    down_extreme <- sig_df %>%
      filter(log2_fc < 0) %>%
      arrange(log2_fc, pvalue_adj) %>%
      head(side_quota)

    extreme_genes <- bind_rows(up_extreme, down_extreme) %>%
      distinct(.gene_key, .keep_all = TRUE)

    max_lfc <- max(abs(sig_df$log2_fc), na.rm = TRUE)
    max_nlp <- max(sig_df$neg_log10_padj, na.rm = TRUE)

    ranked_genes <- sig_df %>%
      mutate(
        score_lfc = abs(log2_fc) / pmax(max_lfc, 1e-9),
        score_padj = neg_log10_padj / pmax(max_nlp, 1e-9),
        score = 0.7 * score_lfc + 0.3 * score_padj
      ) %>%
      arrange(desc(score), pvalue_adj)

    remaining_n <- max(top_n - nrow(extreme_genes), 0L)
    remaining_genes <- if (remaining_n > 0) {
      ranked_genes %>%
        filter(!(.gene_key %in% extreme_genes$.gene_key)) %>%
        head(remaining_n)
    } else {
      ranked_genes[0, , drop = FALSE]
    }

    bind_rows(extreme_genes, remaining_genes) %>%
      distinct(.gene_key, .keep_all = TRUE) %>%
      head(top_n)
  }

  # Volcano plot
  output$dge_volcano_plot <- renderPlot({
    dge_df <- dge_plot_data()
    if (is.null(dge_df)) {
      plot.new(); text(0.5, 0.5, "No data available for this contrast", cex = 1.5); return()
    }
    
    # Create volcano plot
    p <- ggplot(dge_df, aes(x = log2_fc, y = neg_log10_padj)) +
      geom_point(aes(color = direction), alpha = 0.6, size = 2) +
      scale_color_manual(
        values = c(
          "Upregulated" = "#d62728",
          "Downregulated" = "#1f77b4",
          "Not Significant" = "gray70"
        )
      ) +
      geom_hline(yintercept = -log10(input$dge_pval_threshold), 
                 linetype = "dashed", color = "gray30") +
      geom_vline(xintercept = c(-input$dge_logfc_threshold, input$dge_logfc_threshold), 
                 linetype = "dashed", color = "gray30") +
      labs(
        title = sprintf("Volcano Plot: %s", input$dge_contrast_select),
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        color = "Regulation"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        plot.title      = element_text(hjust = 0.5, face = "bold")
      )
    
    # Add gene labels with guaranteed LFC-tail coverage plus weighted ranking.
    if (input$dge_top_n_genes > 0) {
      top_genes <- select_dge_label_genes(dge_df, input$dge_top_n_genes)

      if (nrow(top_genes) > 0) {
        p <- p +
          ggrepel::geom_text_repel(
            data  = top_genes,
            aes(label = gene),
            size  = input$dge_gene_label_size,
            max.overlaps  = 25,
            box.padding   = 0.5,
            point.padding = 0.3
          )
      }
    }

    print(p)
  })

  # DGE hover info panel (gene details when hovering over volcano)
  output$dge_hover_info <- renderUI({
    dge_df <- dge_plot_data()
    if (is.null(dge_df)) return(tags$p(class = "text-muted", "Load DGE results first."))

    hover <- input$plot_hover
    if (is.null(hover)) {
      return(tags$div(
        style = "padding: 12px;",
        icon("arrow-pointer"), tags$span(class = "text-muted ms-1", "Hover over a gene to see details.")
      ))
    }

    if (is.null(hover$x) || is.null(hover$y) || !is.finite(hover$x) || !is.finite(hover$y)) {
      return(tags$p(class = "text-muted small", style = "padding: 8px;", "No gene at cursor."))
    }

    hover_gene <- nearPoints(
      dge_df,
      hover,
      xvar = "log2_fc",
      yvar = "neg_log10_padj",
      threshold = 6,
      maxpoints = 1,
      addDist = FALSE
    )

    if (nrow(hover_gene) == 0) {
      return(tags$p(class = "text-muted small", style = "padding: 8px;", "No gene at cursor."))
    }

    g <- hover_gene[1, , drop = FALSE]
    direction <- if (!is.na(g$significant) && isTRUE(g$significant)) {
      if (g$log2_fc > 0) "Upregulated" else "Downregulated"
    } else "Not Significant"
    dir_color <- switch(direction,
      Upregulated   = "#d62728",
      Downregulated = "#1f77b4",
      "#555"
    )
    tags$div(
      class = "hover-gene-info",
      tags$div(
        style = paste0("font-weight:700; font-size:1.05rem; color:", dir_color, "; margin-bottom:10px;"),
        if (!is.null(g$gene)) g$gene else "Gene"
      ),
      tags$table(
        class = "table table-sm table-bordered hover-gene-table",
        tags$tbody(
          tags$tr(tags$th("Metric"), tags$th("Value")),
          tags$tr(tags$td("Log2 FC"),      tags$td(round(g$log2_fc, 3))),
          tags$tr(tags$td("-log10 padj"),  tags$td(round(g$neg_log10_padj, 2))),
          tags$tr(tags$td("p-value"),      tags$td(if (!is.null(g$pvalue))     format(g$pvalue,     scientific = TRUE, digits = 3) else "—")),
          tags$tr(tags$td("Adj. p-value"), tags$td(format(g$pvalue_adj, scientific = TRUE, digits = 3))),
          if (!is.null(g$mean_expr_group1) && !is.na(g$mean_expr_group1))
            tags$tr(tags$td("Mean group 1"), tags$td(round(g$mean_expr_group1, 3))),
          if (!is.null(g$mean_expr_group2) && !is.na(g$mean_expr_group2))
            tags$tr(tags$td("Mean group 2"), tags$td(round(g$mean_expr_group2, 3))),
          tags$tr(tags$td("Direction"),
            tags$td(tags$span(style = paste0("color:", dir_color, "; font-weight:600;"), direction)))
        )
      )
    )
  })
  
  # Significant genes table
  output$dge_significant_genes_table <- renderDT({
    req(rv_dge$dge_loaded)
    req(input$dge_contrast_select)
    
    # Get selected contrast data
    dge_df <- rv_dge$dge_results[[input$dge_contrast_select]]
    
    if (is.null(dge_df) || nrow(dge_df) == 0) {
      return(data.frame(Message = "No data available"))
    }

    required_cols <- c("log2_fc", "pvalue_adj")
    missing_cols <- setdiff(required_cols, names(dge_df))
    if (length(missing_cols) > 0) {
      return(data.frame(Message = sprintf("Missing required DGE columns: %s", paste(missing_cols, collapse = ", "))))
    }

    dge_df <- dge_df %>%
      filter(!is.na(log2_fc), !is.na(pvalue_adj), pvalue_adj > 0)
    
    # Filter for significant genes
    sig_genes <- dge_df %>%
      filter(
        abs(log2_fc) >= input$dge_logfc_threshold,
        pvalue_adj < input$dge_pval_threshold
      ) %>%
      arrange(pvalue_adj)

    if (nrow(sig_genes) == 0) {
      return(datatable(
        data.frame(Message = "No significant genes found with current thresholds"),
        options = list(dom = 't', ordering = FALSE),
        rownames = FALSE
      ))
    }

    display_cols <- c("gene", "log2_fc", "pvalue", "pvalue_adj", "mean_expr_group1", "mean_expr_group2")
    display_cols <- intersect(display_cols, names(sig_genes))
    sig_genes <- sig_genes %>% select(all_of(display_cols))
    
    # Round numeric columns
    if ("log2_fc" %in% names(sig_genes)) {
      sig_genes$log2_fc <- round(sig_genes$log2_fc, 3)
    }
    if ("pvalue" %in% names(sig_genes)) {
      sig_genes$pvalue <- format(sig_genes$pvalue, scientific = TRUE, digits = 3)
    }
    if ("pvalue_adj" %in% names(sig_genes)) {
      sig_genes$pvalue_adj <- format(sig_genes$pvalue_adj, scientific = TRUE, digits = 3)
    }
    if ("mean_expr_group1" %in% names(sig_genes)) {
      sig_genes$mean_expr_group1 <- round(sig_genes$mean_expr_group1, 3)
    }
    if ("mean_expr_group2" %in% names(sig_genes)) {
      sig_genes$mean_expr_group2 <- round(sig_genes$mean_expr_group2, 3)
    }
    
    datatable(
      sig_genes,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(3, 'asc'))  # Sort by adjusted p-value
      ),
      rownames = FALSE,
      caption = sprintf(
        "Significant genes: %d (Log2FC >= %.2f, Adj. P-value < %.3f)",
        nrow(sig_genes),
        input$dge_logfc_threshold,
        input$dge_pval_threshold
      )
    )
  })
  
  # Download volcano plot
  output$btn_download_volcano <- downloadHandler(
    filename = function() {
      sprintf("volcano_%s.png", input$dge_contrast_select)
    },
    content = function(file) {
      # Get selected contrast data
      dge_df <- rv_dge$dge_results[[input$dge_contrast_select]]

      required_cols <- c("log2_fc", "pvalue_adj")
      if (is.null(dge_df) || nrow(dge_df) == 0 || any(!required_cols %in% names(dge_df))) {
        stop("Selected contrast does not have required columns for volcano plot")
      }

      dge_df <- dge_df %>%
        filter(!is.na(log2_fc), !is.na(pvalue_adj), pvalue_adj > 0)
      
      # Add significance column
      dge_df$significant <- with(dge_df, 
        abs(log2_fc) >= input$dge_logfc_threshold & 
        pvalue_adj < input$dge_pval_threshold
      )
      
      # Add direction column
      dge_df$direction <- ifelse(
        !dge_df$significant, "Not Significant",
        ifelse(dge_df$log2_fc > 0, "Upregulated", "Downregulated")
      )
      
      # Create plot
      p <- ggplot(dge_df, aes(x = log2_fc, y = -log10(pvalue_adj))) +
        geom_point(aes(color = direction), alpha = 0.6, size = 2) +
        scale_color_manual(
          values = c(
            "Upregulated" = "#d62728",
            "Downregulated" = "#1f77b4",
            "Not Significant" = "gray70"
          )
        ) +
        geom_hline(yintercept = -log10(input$dge_pval_threshold), 
                   linetype = "dashed", color = "gray30") +
        geom_vline(xintercept = c(-input$dge_logfc_threshold, input$dge_logfc_threshold), 
                   linetype = "dashed", color = "gray30") +
        labs(
          title = sprintf("Volcano Plot: %s", input$dge_contrast_select),
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted P-value",
          color = "Regulation"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "bottom",
          plot.title      = element_text(hjust = 0.5, face = "bold")
        )
      
      # Add gene labels with guaranteed LFC-tail coverage plus weighted ranking.
      if (input$dge_top_n_genes > 0) {
        top_genes <- select_dge_label_genes(dge_df, input$dge_top_n_genes)

        if (nrow(top_genes) > 0) {
          p <- p +
            ggrepel::geom_text_repel(
              data  = top_genes,
              aes(label = gene),
              size  = input$dge_gene_label_size,
              max.overlaps = 25
            )
        }
      }
      
      # Save to file
      ggsave(file, plot = p, width = 8, height = 8, dpi = 300)
    }
  )
  
  # Download significant genes CSV
  output$btn_download_sig_genes <- downloadHandler(
    filename = function() {
      sprintf("significant_genes_%s.csv", input$dge_contrast_select)
    },
    content = function(file) {
      dge_df <- rv_dge$dge_results[[input$dge_contrast_select]]

      required_cols <- c("log2_fc", "pvalue_adj")
      if (is.null(dge_df) || nrow(dge_df) == 0 || any(!required_cols %in% names(dge_df))) {
        stop("Selected contrast does not have required columns for filtering")
      }

      dge_df <- dge_df %>%
        filter(!is.na(log2_fc), !is.na(pvalue_adj), pvalue_adj > 0)
      
      sig_genes <- dge_df %>%
        filter(
          abs(log2_fc) >= input$dge_logfc_threshold,
          pvalue_adj < input$dge_pval_threshold
        ) %>%
        arrange(pvalue_adj)
      
      write.csv(sig_genes, file, row.names = FALSE)
    }
  )
  
  # Download all results CSV
  output$btn_download_all_results <- downloadHandler(
    filename = function() {
      sprintf("all_results_%s.csv", input$dge_contrast_select)
    },
    content = function(file) {
      dge_df <- rv_dge$dge_results[[input$dge_contrast_select]]
      write.csv(dge_df, file, row.names = FALSE)
    }
  )

  # ===========================================================================
  # TAB 6: GSEA
  # ===========================================================================

  rv_gsea <- reactiveValues(
    gsea_dir = NULL,
    contrasts = character(0),
    data = list(),
    max_pathways = 50,
    loaded = FALSE
  )

  get_selected_gsea <- reactive({
    req(rv_gsea$loaded)
    req(input$gsea_contrast_select)
    req(input$gsea_db_select)

    contrast_data <- rv_gsea$data[[input$gsea_contrast_select]]
    if (is.null(contrast_data)) {
      return(NULL)
    }

    obj <- contrast_data$gsea_results[[input$gsea_db_select]]
    if (is.null(obj)) {
      return(NULL)
    }

    # Filter by adj. p-value cutoff from slider
    padj_cutoff <- input$gsea_padj_cutoff
    if (!is.null(padj_cutoff) && padj_cutoff < 1.0) {
      obj <- tryCatch(
        clusterProfiler::filter(obj, p.adjust <= padj_cutoff),
        error = function(e) obj
      )
    }

    obj
  })

  load_gsea_results <- function(gsea_dir, show_progress = TRUE) {
    if (!dir.exists(gsea_dir)) {
      stop(sprintf("Directory not found: %s", gsea_dir))
    }

    contrast_dirs <- list.dirs(gsea_dir, recursive = FALSE, full.names = TRUE)
    contrast_dirs <- contrast_dirs[grepl("_gsea$", basename(contrast_dirs))]

    if (length(contrast_dirs) == 0) {
      stop("No *_gsea directories found")
    }

    gsea_data <- list()
    contrast_names <- character(0)
    detected_max <- 50

    for (i in seq_along(contrast_dirs)) {
      contrast_dir <- contrast_dirs[[i]]
      contrast_name <- sub("_gsea$", "", basename(contrast_dir))
      rds_file <- file.path(contrast_dir, "gsea_dashboard_data.rds")

      if (!file.exists(rds_file)) {
        next
      }

      obj <- readRDS(rds_file)
      gsea_data[[contrast_name]] <- obj
      contrast_names <- c(contrast_names, contrast_name)

      if (!is.null(obj$max_pathways_dashboard)) {
        detected_max <- max(detected_max, as.integer(obj$max_pathways_dashboard))
      }

      if (show_progress) {
        incProgress(0.75 / max(length(contrast_dirs), 1))
      }
    }

    if (length(contrast_names) == 0) {
      stop("No valid gsea_dashboard_data.rds found in *_gsea directories")
    }

    rv_gsea$gsea_dir <- gsea_dir
    rv_gsea$contrasts <- contrast_names
    rv_gsea$data <- gsea_data
    rv_gsea$loaded <- TRUE
    rv_gsea$max_pathways <- detected_max

    # Register gsea dir as static resource so pathview PNGs can be served
    addResourcePath("gsea_results", gsea_dir)

    updateSelectInput(
      session,
      "gsea_contrast_select",
      choices = contrast_names,
      selected = contrast_names[1]
    )

    # Dynamically update database selector to only show databases present in the RDS
    available_dbs <- intersect(
      names(gsea_data[[contrast_names[1]]]$gsea_results),
      c("go_bp", "go_mf", "go_cc", "kegg", "reactome")
    )
    db_label_map <- c(
      go_bp    = "GO Biological Process",
      go_mf    = "GO Molecular Function",
      go_cc    = "GO Cellular Component",
      kegg     = "KEGG",
      reactome = "Reactome"
    )
    if (length(available_dbs) > 0) {
      db_choices <- setNames(available_dbs, db_label_map[available_dbs])
      updateSelectInput(session, "gsea_db_select",
                        choices  = db_choices,
                        selected = available_dbs[1])
    }

    default_top <- gsea_data[[contrast_names[1]]]$top_pathways_default
    if (is.null(default_top) || is.na(default_top)) {
      default_top <- 10
    }

    updateNumericInput(
      session,
      "gsea_n_pathways",
      value = min(as.integer(default_top), detected_max),
      max   = detected_max
    )
  }

  # (GSEA loading is handled by the startup auto-load observer)

  # Dynamic GSEA navset panel (Dotplot / Ridgeplot / GSEA Plot / Table / Pathview / TreeDot)
  output$gsea_main_panels <- renderUI({
    if (!isTRUE(rv_gsea$loaded)) {
      return(card(
        card_header("GSEA"),
        tags$p(class = "text-muted p-3",
               "GSEA results are being loaded. If nothing appears, ensure the pipeline has run with GSEA enabled.")
      ))
    }

    is_kegg <- isTRUE(input$gsea_db_select == "kegg")

    base_panels <- list(
      nav_panel("Dot Plot",
                plotOutput("gsea_dotplot", height = "450px")),
      nav_panel("Ridgeplot",
                plotOutput("gsea_ridgeplot", height = "420px")),
      nav_panel("GSEA Plot",
                plotOutput("gsea_multiplot", height = "440px")),
      nav_panel("Results Table",
                DTOutput("gsea_table"))
    )

    if (is_kegg) {
      base_panels <- c(base_panels, list(
        nav_panel("Pathview",
                  tags$div(
                    style = "padding: 8px;",
                    selectInput("gsea_pathview_select",
                                "Select Pathway:",
                                choices  = c("Loading..." = ""),
                                selected = ""),
                    uiOutput("gsea_pathview")
                  ))
      ))
    }

    base_panels <- c(base_panels, list(
      nav_panel("TreeDot",
                tags$p(class = "text-muted p-3",
                       "TreeDot visualization requires a multi-contrast ORA analysis. Not yet available for this dataset."))
    ))

    do.call(navset_card_tab, base_panels)
  })

  output$gsea_dotplot <- renderPlot({
    if (!isTRUE(rv_gsea$loaded)) return(invisible(NULL))
    gsea_obj <- get_selected_gsea()
    if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
      plot.new()
      text(0.5, 0.5, "No enriched pathways for selection", cex = 1.2)
      return()
    }

    n_show <- min(input$gsea_n_pathways, nrow(as.data.frame(gsea_obj)))
    p <- enrichplot::dotplot(gsea_obj, showCategory = n_show) +
      ggplot2::ggtitle(sprintf("%s: Dotplot (%d pathways)", input$gsea_contrast_select, n_show))
    print(p)
  })

  output$gsea_ridgeplot <- renderPlot({
    if (!isTRUE(rv_gsea$loaded)) return(invisible(NULL))

    contrast_data <- rv_gsea$data[[input$gsea_contrast_select]]
    db_key        <- input$gsea_db_select

    # Use orig (non-setReadable) object to preserve gene IDs for ridgeplot
    orig_obj <- if (!is.null(contrast_data$gsea_results_orig))
      contrast_data$gsea_results_orig[[db_key]] else NULL
    if (is.null(orig_obj))
      orig_obj <- contrast_data$gsea_results[[db_key]]

    if (is.null(orig_obj)) {
      plot.new(); text(0.5, 0.5, "No results for selected database", cex = 1.2); return()
    }

    # Apply padj filter to orig object
    padj_cutoff   <- input$gsea_padj_cutoff
    orig_filtered <- tryCatch(
      clusterProfiler::filter(orig_obj, p.adjust <= padj_cutoff),
      error = function(e) orig_obj
    )

    if (nrow(as.data.frame(orig_filtered)) == 0) {
      plot.new(); text(0.5, 0.5, "No enriched pathways for selection", cex = 1.2); return()
    }

    ridge_df <- as.data.frame(orig_filtered)
    if (!("core_enrichment" %in% names(ridge_df))) {
      plot.new(); text(0.5, 0.5, "Ridgeplot not available (no core_enrichment column)", cex = 1.2); return()
    }

    # Use character vector of Descriptions (not IDs) to avoid showCategory matching bug
    n_show   <- min(input$gsea_n_pathways, nrow(ridge_df))
    top_descs <- head(ridge_df$Description[order(ridge_df$p.adjust)], n_show)

    tryCatch({
      p <- enrichplot::ridgeplot(orig_filtered, showCategory = top_descs) +
        ggplot2::ggtitle(sprintf("%s: Ridgeplot (%d pathways)", input$gsea_contrast_select, length(top_descs)))
      print(p)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Ridgeplot error:", conditionMessage(e)), cex = 0.9)
    })
  })

  output$gsea_multiplot <- renderPlot({
    if (!isTRUE(rv_gsea$loaded)) return(invisible(NULL))
    gsea_obj <- get_selected_gsea()
    if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
      plot.new()
      text(0.5, 0.5, "No enriched pathways for selection", cex = 1.2)
      return()
    }

    n_show <- min(input$gsea_n_pathways, nrow(as.data.frame(gsea_obj)))
    idx <- seq_len(n_show)
    p <- enrichplot::gseaplot2(
      gsea_obj,
      geneSetID = idx,
      title = sprintf("%s: Running score (%d pathways)", input$gsea_contrast_select, n_show),
      color = grDevices::hcl.colors(n_show, "Dark 3")
    )
    print(p)
  })

  # Update pathview pathway selector when contrast changes
  observe({
    req(rv_gsea$loaded)
    req(input$gsea_contrast_select)

    contrast_name <- input$gsea_contrast_select
    pathview_dir <- file.path(rv_gsea$gsea_dir, paste0(contrast_name, "_gsea"), "pathview")

    if (!dir.exists(pathview_dir)) {
      updateSelectInput(session, "gsea_pathview_select", choices = c("No images available" = ""))
      return()
    }

    png_files <- list.files(pathview_dir, pattern = paste0("\\.", contrast_name, "\\.png$"))

    if (length(png_files) == 0) {
      updateSelectInput(session, "gsea_pathview_select", choices = c("No images available" = ""))
      return()
    }

    pathway_ids <- sub(paste0("\\.", contrast_name, "\\.png$"), "", png_files)

    # Order by enrichment (absolute NES descending = most enriched first)
    kegg_table <- rv_gsea$data[[contrast_name]]$tables$GSEA_KEGG
    if (!is.null(kegg_table) && nrow(kegg_table) > 0 && "NES" %in% names(kegg_table)) {
      kegg_sorted <- kegg_table[order(abs(kegg_table$NES), decreasing = TRUE), ]
      ordered_ids <- kegg_sorted$ID[kegg_sorted$ID %in% pathway_ids]
      remaining_ids <- pathway_ids[!pathway_ids %in% ordered_ids]
      pathway_ids <- c(ordered_ids, remaining_ids)
      png_files <- paste0(pathway_ids, ".", contrast_name, ".png")

      # Build labels with description
      desc_map <- setNames(as.character(kegg_sorted$Description), as.character(kegg_sorted$ID))
      labels <- ifelse(
        pathway_ids %in% names(desc_map),
        paste0(pathway_ids, " - ", desc_map[pathway_ids]),
        pathway_ids
      )
    } else {
      labels <- pathway_ids
    }

    choices <- setNames(png_files, labels)
    updateSelectInput(session, "gsea_pathview_select", choices = choices, selected = png_files[1])
  })

  output$gsea_pathview <- renderUI({
    req(rv_gsea$loaded)
    req(input$gsea_contrast_select)
    req(input$gsea_pathview_select)
    req(nzchar(input$gsea_pathview_select))

    contrast_name <- input$gsea_contrast_select
    src_url <- paste0("gsea_results/", contrast_name, "_gsea/pathview/", input$gsea_pathview_select)

    tags$div(
      tags$img(
        src = src_url,
        style = "width:100%; max-width:900px; display:block; margin:0 auto; border:1px solid #ddd; border-radius:4px;"
      )
    )
  })

  output$gsea_table <- renderDT({
    gsea_obj <- get_selected_gsea()
    if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
      return(datatable(
        data.frame(Message = "No enriched pathways for selected contrast/database"),
        options = list(dom = "t", ordering = FALSE),
        rownames = FALSE
      ))
    }

    gsea_df <- as.data.frame(gsea_obj)
    gsea_df <- head(gsea_df, input$gsea_n_pathways)
    keep <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues")
    keep <- keep[keep %in% names(gsea_df)]
    gsea_df <- gsea_df[, keep, drop = FALSE]

    # Format numeric columns: 3 decimals for scores, scientific for p-values
    num_cols <- c("enrichmentScore", "NES", "qvalues")
    num_cols <- num_cols[num_cols %in% names(gsea_df)]
    for (col in num_cols) {
      gsea_df[[col]] <- round(as.numeric(gsea_df[[col]]), 3)
    }
    sci_cols <- c("pvalue", "p.adjust")
    sci_cols <- sci_cols[sci_cols %in% names(gsea_df)]
    for (col in sci_cols) {
      gsea_df[[col]] <- formatC(as.numeric(gsea_df[[col]]), digits = 3, format = "e")
    }

    datatable(
      gsea_df,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # ===========================================================================
  # TAB 7: ANNOTATION STATION
  # ===========================================================================
  
  # Reactive value to store custom annotation
  rv_annotation <- reactiveValues(
    labels = NULL,
    annotation_name = NULL,
    rules_df = NULL
  )
  
  # Plot custom annotation button
  observeEvent(input$btn_plot_annotation, {
    
    req(rv$data_loaded)
    req(input$annot_name)
    req(input$annot_rules)
    
    isolate({
      
      tryCatch({
        
        # Parse annotation rules
        rules_df <- parse_annotation_rules(input$annot_rules)
        
        if (nrow(rules_df) == 0) {
          showNotification(
            "No valid annotation rules found. Please check your input format.",
            type = "warning",
            duration = 5
          )
          return(NULL)
        }
        
        # Apply custom annotation
        labels <- apply_custom_annotation(
          metadata = rv$data_obj$metadata,
          rules_df = rules_df,
          annotation_name = input$annot_name
        )
        
        # Store in reactive values
        rv_annotation$labels <- labels
        rv_annotation$annotation_name <- input$annot_name
        rv_annotation$rules_df <- rules_df
        
        showNotification(
          sprintf("Custom annotation '%s' created successfully!", input$annot_name),
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error creating annotation:", e$message),
          type = "error",
          duration = 5
        )
        rv_annotation$labels <- NULL
      })
    })
  })
  
  # Render annotation UMAP
  output$annotation_umap <- renderPlotly({
    req(rv$data_loaded)
    req(rv$data_obj$umap_coords)
    
    # Check if we have custom annotation labels, otherwise show all as "Unknown"
    if (!is.null(rv_annotation$labels)) {
      # Use custom annotation labels
      umap_data <- rv$data_obj$umap_coords
      umap_data$annotation <- rv_annotation$labels[umap_data$cell_id]
      plot_title <- sprintf("Custom Annotation: %s", rv_annotation$annotation_name)
    } else {
      # Show all cells as "Unknown" by default
      umap_data <- rv$data_obj$umap_coords
      umap_data$annotation <- "Unknown"
      plot_title <- "Custom Annotation (all cells marked as Unknown)"
    }
    
    # Create plot
    p <- plot_ly(
      data = umap_data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = 'scattergl',
      mode = 'markers',
      marker = list(
        size = input$annot_point_size,
        opacity = input$annot_opacity
      ),
      color  = ~annotation,
      colors = SC_COLORS,
      text   = ~paste("Cell:", cell_id, "<br>Label:", annotation),
      hoverinfo = 'text'
    ) %>%
      layout(
        title  = plot_title,
        xaxis  = list(title = "UMAP 1", constrain = "domain"),
        yaxis  = list(title = "UMAP 2", scaleanchor = "x", scaleratio = 1, constrain = "domain"),
        legend = list(x = 1.02, y = 0.5, xanchor = "left", yanchor = "middle",
                      orientation = "v", itemsizing = "constant"),
        margin    = list(r = 160, b = 40),
        hovermode = 'closest'
      ) %>% config(responsive = TRUE)
    
    return(p)
  })
  
  # Render annotation statistics
  output$annotation_stats <- renderText({
    
    # Show initial message if no annotation has been created
    if (is.null(rv_annotation$labels)) {
      if (!rv$data_loaded) {
        return("Load data first to start creating custom annotations.")
      } else {
        n_total <- nrow(rv$data_obj$metadata)
        return(sprintf(
          "Ready to annotate\n\nTotal cells: %s\n\nEnter annotation rules above and click 'Plot' to create a custom annotation.",
          format_number(n_total)
        ))
      }
    }
    
    # Calculate statistics
    label_counts <- table(rv_annotation$labels)
    n_total <- length(rv_annotation$labels)
    n_unknown <- sum(rv_annotation$labels == "Unknown")
    n_annotated <- n_total - n_unknown
    
    # Calculate percentages for each label
    label_percentages <- (label_counts / n_total) * 100
    
    # Format label distribution with counts and percentages
    label_dist_text <- paste(
      sprintf("  %s: %s (%.1f%%)", 
              names(label_counts), 
              format_number(label_counts),
              label_percentages),
      collapse = "\n"
    )
    
    # Format output
    stats_text <- sprintf(
      "Annotation: %s\n\nTotal cells: %s\nAnnotated: %s (%.1f%%)\nUnknown: %s (%.1f%%)\nUnique labels: %d\n\nLabel distribution:\n%s",
      rv_annotation$annotation_name,
      format_number(n_total),
      format_number(n_annotated),
      100 * n_annotated / n_total,
      format_number(n_unknown),
      100 * n_unknown / n_total,
      length(unique(rv_annotation$labels)),
      label_dist_text
    )
    
    return(stats_text)
  })
  
  # Save annotation to H5AD button
  observeEvent(input$btn_save_annotation, {
    
    req(rv_annotation$labels)
    req(rv$data_obj)
    req(nzchar(rv$h5ad_path %||% ""))
    
    isolate({
      
      withProgress(message = 'Saving annotation...', value = 0, {
        
        tryCatch({
          
          incProgress(0.3, detail = "Preparing to save...")
          
          # Determine save mode
          create_copy <- (input$annot_save_mode == "create_copy")
          
          incProgress(0.3, detail = "Writing to H5AD...")
          
          # Save annotation
          output_path <- save_annotation_to_h5ad(
            h5ad_path       = rv$h5ad_path,
            annotation_name = rv_annotation$annotation_name,
            labels          = rv_annotation$labels,
            output_path     = NULL,
            create_copy     = create_copy
          )
          
          incProgress(0.4, detail = "Done!")
          
          # Show success message
          msg <- if (create_copy) {
            sprintf("Annotation saved to new file:\n%s", basename(output_path))
          } else {
            sprintf("Annotation saved to original file:\n%s", basename(output_path))
          }
          
          showNotification(
            msg,
            type = "message",
            duration = 10
          )
          
          # If we created a copy, ask if user wants to reload
          if (create_copy) {
            showNotification(
              "Tip: You can reload the new file from the Data Input tab to see the saved annotation.",
              type = "message",
              duration = 10
            )
          } else {
            # If we overwrote, reload the data
            showNotification(
              "Reloading data with new annotation...",
              type = "message",
              duration = 3
            )
            
            # Trigger data reload
            Sys.sleep(1)
            
            rv$data_obj <- load_h5ad_data(
              rv$h5ad_path,
              backed = FALSE
            )
            
            # Update color choices
            rv$umap_color_choices <- setdiff(
              names(rv$data_obj$metadata),
              c("cell_id")
            )
            
            updateSelectInput(
              session,
              "umap_color_by",
              choices = rv$umap_color_choices,
              selected = if(rv_annotation$annotation_name %in% rv$umap_color_choices) {
                rv_annotation$annotation_name
              } else if("batch" %in% rv$umap_color_choices) {
                "batch"
              } else {
                rv$umap_color_choices[1]
              }
            )
          }
          
        }, error = function(e) {
          showNotification(
            paste("Error saving annotation:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    })
  })
  
}
