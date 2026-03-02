# scAnnex Dashboard - Server Logic
# Handles reactive data processing and plot generation

server <- function(input, output, session) {
  
  # ===========================================================================
  # REACTIVE VALUES
  # ===========================================================================
  
  rv <- reactiveValues(
    data_obj = NULL,
    qc_report = NULL,
    qc_reports_multi = list(),
    qc_plots = list(),
    data_loaded = FALSE,
    umap_color_choices = character(0)
  )

  # Initialize sample selector choices
  observeEvent(TRUE, {
    if (length(DEFAULT_SAMPLE_IDS) > 0) {
      updateSelectInput(
        session,
        "sample_select",
        choices = DEFAULT_SAMPLE_IDS,
        selected = DEFAULT_SAMPLE_IDS[1]
      )
    }
  }, once = TRUE)

  # Auto-populate paths based on selected mode/sample
  observe({
    req(input$data_mode)

    if (identical(input$data_mode, "integrated")) {
      if (nzchar(DEFAULT_MERGED_H5AD)) {
        updateTextInput(session, "input_h5ad_path", value = DEFAULT_MERGED_H5AD)
      }
      if (nzchar(DEFAULT_QC_DIR)) {
        updateTextInput(session, "input_qc_dir", value = DEFAULT_QC_DIR)
      }
      if (nzchar(DEFAULT_DGE_DIR)) {
        updateTextInput(session, "input_dge_dir", value = DEFAULT_DGE_DIR)
      }
    }

    if (identical(input$data_mode, "single")) {
      req(input$sample_select)
      if (input$sample_select %in% names(DEFAULT_SAMPLE_H5AD_FILES)) {
        updateTextInput(
          session,
          "input_h5ad_path",
          value = DEFAULT_SAMPLE_H5AD_FILES[[input$sample_select]]
        )
      }
      if (nzchar(DEFAULT_QC_DIR)) {
        updateTextInput(session, "input_qc_dir", value = DEFAULT_QC_DIR)
      }
    }
  })
  
  # ===========================================================================
  # DATA LOADING
  # ===========================================================================
  
  # Load data when button clicked
  observeEvent(input$btn_load_data, {

    selected_h5ad_path <- input$input_h5ad_path
    selected_qc_dir <- input$input_qc_dir

    if (identical(input$data_mode, "single") &&
        !is.null(input$sample_select) &&
        input$sample_select %in% names(DEFAULT_SAMPLE_H5AD_FILES)) {
      selected_h5ad_path <- DEFAULT_SAMPLE_H5AD_FILES[[input$sample_select]]
    }

    req(selected_h5ad_path)
    
    withProgress(message = 'Loading data...', value = 0, {
      
      tryCatch({
        
        # Check file exists
        if (!file.exists(selected_h5ad_path)) {
          stop(sprintf("File not found: %s", selected_h5ad_path))
        }
        
        incProgress(0.2, detail = "Reading H5AD file")
        
        # Load H5AD data
        rv$data_obj <- load_h5ad_data(
          selected_h5ad_path,
          backed = input$input_backed_mode
        )
        
        incProgress(0.4, detail = "Loading QC report")
        
        # Load QC report if directory exists
        if (dir.exists(selected_qc_dir)) {
          rv$qc_plots <- get_qc_plots(selected_qc_dir)

          if (identical(input$data_mode, "integrated")) {
            normalized_qc_dir <- normalizePath(selected_qc_dir, mustWork = TRUE)
            results_dir <- if (basename(normalized_qc_dir) == "qc") {
              dirname(normalized_qc_dir)
            } else {
              normalized_qc_dir
            }

            rv$qc_reports_multi <- load_all_qc_reports(results_dir)

            if (length(rv$qc_reports_multi) > 0) {
              sample_choices <- names(rv$qc_reports_multi)
              updateSelectInput(
                session,
                "qc_sample_select",
                choices = sample_choices,
                selected = sample_choices[1]
              )
              rv$qc_report <- rv$qc_reports_multi[[sample_choices[1]]]

              sample_plots <- get_qc_plots_for_sample(results_dir, sample_choices[1])
              if (length(sample_plots) > 0) {
                rv$qc_plots <- sample_plots
              }
            } else {
              rv$qc_report <- load_qc_report(selected_qc_dir)
            }
          } else {
            rv$qc_reports_multi <- list()
            rv$qc_report <- load_qc_report(selected_qc_dir)
          }
        } else {
          rv$qc_report <- NULL
          rv$qc_reports_multi <- list()
          rv$qc_plots <- list()
        }
        
        incProgress(0.2, detail = "Preparing visualization")
        
        # Update color choices for UMAP
        rv$umap_color_choices <- setdiff(
          names(rv$data_obj$metadata),
          c("cell_id")
        )
        
        # Update selectInput choices
        updateSelectInput(
          session,
          "umap_color_by",
          choices = rv$umap_color_choices,
          selected = if("batch" %in% rv$umap_color_choices) "batch" else rv$umap_color_choices[1]
        )
        
        incProgress(0.2, detail = "Done!")
        
        rv$data_loaded <- TRUE
        
        showNotification(
          "Data loaded successfully!",
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error loading data:", e$message),
          type = "error",
          duration = 10
        )
        rv$data_loaded <- FALSE
      })
    })
  })
  
  # Data load status text
  output$data_load_status <- renderText({
    if (rv$data_loaded) {
      sprintf(
        "✓ Data loaded successfully\n\nMode: %s\nDataset: %s cells × %s genes\nBacked mode: %s\nQC report: %s",
        ifelse(identical(input$data_mode, "single"), "Single sample", "Integrated"),
        format_number(rv$data_obj$n_cells),
        format_number(rv$data_obj$n_genes),
        ifelse(rv$data_obj$backed, "Yes", "No"),
        ifelse(is.null(rv$qc_report), "Not available", "Loaded")
      )
    } else {
      "No data loaded. Click 'Load Data' to begin."
    }
  })
  
  # Sidebar dataset info
  output$sidebar_dataset_info <- renderText({
    if (rv$data_loaded) {
      sprintf(
        "%s cells\n%s genes",
        format_number(rv$data_obj$n_cells),
        format_number(rv$data_obj$n_genes)
      )
    } else {
      "No data loaded"
    }
  })
  
  # ===========================================================================
  # TAB 2: QC OVERVIEW
  # ===========================================================================

  qc_report_active <- reactive({
    if (identical(input$data_mode, "integrated") && length(rv$qc_reports_multi) > 0) {
      selected_sample <- input$qc_sample_select
      if (!is.null(selected_sample) && selected_sample %in% names(rv$qc_reports_multi)) {
        return(rv$qc_reports_multi[[selected_sample]])
      }
    }

    rv$qc_report
  })

  observeEvent(input$qc_sample_select, {
    req(identical(input$data_mode, "integrated"))
    req(input$input_qc_dir)
    req(input$qc_sample_select)

    if (!dir.exists(input$input_qc_dir)) {
      return()
    }

    normalized_qc_dir <- normalizePath(input$input_qc_dir, mustWork = TRUE)
    results_dir <- if (basename(normalized_qc_dir) == "qc") {
      dirname(normalized_qc_dir)
    } else {
      normalized_qc_dir
    }

    sample_plots <- get_qc_plots_for_sample(results_dir, input$qc_sample_select)
    if (length(sample_plots) > 0) {
      rv$qc_plots <- sample_plots
    }
  }, ignoreInit = TRUE)
  
  # QC Info Boxes
  output$qc_box_cells_before <- renderInfoBox({
    report <- qc_report_active()
    req(report)
    
    n_cells <- report$filtering_statistics$cells_initial
    
    infoBox(
      "Cells (Before QC)",
      format_number(n_cells),
      icon = icon("circle"),
      color = "blue"
    )
  })
  
  output$qc_box_cells_after <- renderInfoBox({
    report <- qc_report_active()
    req(report)
    
    n_cells <- report$filtering_statistics$cells_final
    
    infoBox(
      "Cells (After QC)",
      format_number(n_cells),
      icon = icon("check-circle"),
      color = "green"
    )
  })
  
  output$qc_box_genes_after <- renderInfoBox({
    report <- qc_report_active()
    req(report)
    
    n_genes <- report$filtering_statistics$genes_final
    
    infoBox(
      "Genes (After QC)",
      format_number(n_genes),
      icon = icon("dna"),
      color = "purple"
    )
  })
  
  output$qc_box_retention <- renderInfoBox({
    report <- qc_report_active()
    req(report)
    
    retention_pct <- report$filtering_statistics$cells_retained_pct
    
    infoBox(
      "Cell Retention",
      sprintf("%.1f%%", retention_pct),
      icon = icon("percentage"),
      color = "yellow"
    )
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
    
    # Create plot
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
      color = as.formula(paste0("~", input$umap_color_by)),
      text = ~paste("Cell:", cell_id),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = sprintf("UMAP colored by %s", input$umap_color_by),
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2"),
        hovermode = 'closest'
      )
    
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
              color = ~expression,
              colorscale = 'Viridis',
              showscale = TRUE,
              colorbar = list(title = "Expression")
            ),
            text = ~paste("Cell:", cell_id, "<br>Expression:", round(expression, 2)),
            hoverinfo = 'text'
          ) %>%
            layout(
              title = sprintf("Expression of %s", gene_name),
              xaxis = list(title = "UMAP 1"),
              yaxis = list(title = "UMAP 2"),
              hovermode = 'closest'
            )
          
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
              title = sprintf("Gene Set Score (%d genes)", length(gene_list)),
              xaxis = list(title = "UMAP 1"),
              yaxis = list(title = "UMAP 2"),
              hovermode = 'closest'
            )
          
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
  
  # Load DGE results
  observeEvent(input$btn_load_dge, {
    
    req(input$input_dge_dir)
    
    withProgress(message = 'Loading DGE results...', value = 0, {
      
      tryCatch({
        load_dge_results(input$input_dge_dir, show_progress = TRUE)
        
        showNotification(
          sprintf("Loaded %d contrasts successfully", length(rv_dge$contrasts)),
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error loading DGE results:", e$message),
          type = "error",
          duration = 10
        )
        rv_dge$dge_loaded <- FALSE
      })
    })
  })

  # Auto-load DGE results for integrated mode when data is loaded
  observeEvent(rv$data_loaded, {
    if (!isTRUE(rv$data_loaded) || !identical(input$data_mode, "integrated")) {
      return()
    }

    dge_dir <- input$input_dge_dir
    if (!rv_dge$dge_loaded && !is.null(dge_dir) && nzchar(dge_dir) && dir.exists(dge_dir)) {
      tryCatch({
        load_dge_results(dge_dir, show_progress = FALSE)
        showNotification(
          sprintf("Auto-loaded DGE results (%d contrasts)", length(rv_dge$contrasts)),
          type = "message",
          duration = 3
        )
      }, error = function(e) {
        showNotification(
          paste("DGE auto-load skipped:", e$message),
          type = "warning",
          duration = 5
        )
      })
    }
  }, ignoreInit = TRUE)
  
  # Volcano plot
  output$dge_volcano_plot <- renderPlot({
    req(rv_dge$dge_loaded)
    req(input$dge_contrast_select)
    
    # Get selected contrast data
    dge_df <- rv_dge$dge_results[[input$dge_contrast_select]]
    
    if (is.null(dge_df) || nrow(dge_df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data available for this contrast", cex = 1.5)
      return()
    }

    required_cols <- c("log2_fc", "pvalue_adj")
    missing_cols <- setdiff(required_cols, names(dge_df))
    if (length(missing_cols) > 0) {
      plot.new()
      text(0.5, 0.5, sprintf("Missing required DGE columns: %s", paste(missing_cols, collapse = ", ")), cex = 1.2)
      return()
    }

    dge_df <- dge_df %>%
      filter(!is.na(log2_fc), !is.na(pvalue_adj), pvalue_adj > 0)

    if (nrow(dge_df) == 0) {
      plot.new()
      text(0.5, 0.5, "No valid points available for volcano plot", cex = 1.2)
      return()
    }
    
    # Add significance column
    dge_df$significant <- with(dge_df, 
      abs(log2_fc) >= input$dge_logfc_threshold & 
      pvalue_adj < input$dge_pval_threshold
    )
    
    # Add direction column for coloring
    dge_df$direction <- ifelse(
      !dge_df$significant, "Not Significant",
      ifelse(dge_df$log2_fc > 0, "Upregulated", "Downregulated")
    )
    
    # Create volcano plot
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
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    # Add gene labels if requested
    if (input$dge_show_gene_names && input$dge_top_n_genes > 0) {
      
      # Get top N significant genes by p-value
      top_genes <- dge_df %>%
        filter(significant) %>%
        arrange(pvalue_adj) %>%
        head(input$dge_top_n_genes)
      
      if (nrow(top_genes) > 0) {
        p <- p + 
          ggrepel::geom_text_repel(
            data = top_genes,
            aes(label = gene),
            size = input$dge_gene_label_size,
            max.overlaps = 20,
            box.padding = 0.5,
            point.padding = 0.3
          )
      }
    }
    
    print(p)
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
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
      # Add gene labels if requested
      if (input$dge_show_gene_names && input$dge_top_n_genes > 0) {
        top_genes <- dge_df %>%
          filter(significant) %>%
          arrange(pvalue_adj) %>%
          head(input$dge_top_n_genes)
        
        if (nrow(top_genes) > 0) {
          p <- p + 
            ggrepel::geom_text_repel(
              data = top_genes,
              aes(label = gene),
              size = input$dge_gene_label_size,
              max.overlaps = 20
            )
        }
      }
      
      # Save to file
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
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
  # TAB 6: ANNOTATION STATION
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
      color = ~annotation,
      text = ~paste("Cell:", cell_id, "<br>Label:", annotation),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = plot_title,
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2"),
        hovermode = 'closest'
      )
    
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
    req(input$input_h5ad_path)
    
    isolate({
      
      withProgress(message = 'Saving annotation...', value = 0, {
        
        tryCatch({
          
          incProgress(0.3, detail = "Preparing to save...")
          
          # Determine save mode
          create_copy <- (input$annot_save_mode == "create_copy")
          
          incProgress(0.3, detail = "Writing to H5AD...")
          
          # Save annotation
          output_path <- save_annotation_to_h5ad(
            h5ad_path = input$input_h5ad_path,
            annotation_name = rv_annotation$annotation_name,
            labels = rv_annotation$labels,
            output_path = NULL,
            create_copy = create_copy
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
              input$input_h5ad_path,
              backed = input$input_backed_mode
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
