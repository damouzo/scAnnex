# scAnnex Dashboard - Server Logic
# Handles reactive data processing and plot generation

server <- function(input, output, session) {
  
  # ===========================================================================
  # REACTIVE VALUES
  # ===========================================================================
  
  rv <- reactiveValues(
    data_obj = NULL,
    qc_report = NULL,
    qc_plots = list(),
    data_loaded = FALSE,
    umap_color_choices = character(0)
  )
  
  # ===========================================================================
  # DATA LOADING
  # ===========================================================================
  
  # Load data when button clicked
  observeEvent(input$btn_load_data, {
    
    req(input$input_h5ad_path)
    
    withProgress(message = 'Loading data...', value = 0, {
      
      tryCatch({
        
        # Check file exists
        if (!file.exists(input$input_h5ad_path)) {
          stop(sprintf("File not found: %s", input$input_h5ad_path))
        }
        
        incProgress(0.2, detail = "Reading H5AD file")
        
        # Load H5AD data
        rv$data_obj <- load_h5ad_data(
          input$input_h5ad_path,
          backed = input$input_backed_mode
        )
        
        incProgress(0.4, detail = "Loading QC report")
        
        # Load QC report if directory exists
        if (dir.exists(input$input_qc_dir)) {
          rv$qc_report <- load_qc_report(input$input_qc_dir)
          rv$qc_plots <- get_qc_plots(input$input_qc_dir)
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
        "✓ Data loaded successfully\n\nDataset: %s cells × %s genes\nBacked mode: %s\nQC report: %s",
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
  
  # QC Info Boxes
  output$qc_box_cells_before <- renderInfoBox({
    req(rv$qc_report)
    
    n_cells <- rv$qc_report$filtering_statistics$cells_before
    
    infoBox(
      "Cells (Before QC)",
      format_number(n_cells),
      icon = icon("circle"),
      color = "blue"
    )
  })
  
  output$qc_box_cells_after <- renderInfoBox({
    req(rv$qc_report)
    
    n_cells <- rv$qc_report$filtering_statistics$cells_after
    
    infoBox(
      "Cells (After QC)",
      format_number(n_cells),
      icon = icon("check-circle"),
      color = "green"
    )
  })
  
  output$qc_box_genes_after <- renderInfoBox({
    req(rv$qc_report)
    
    n_genes <- rv$qc_report$filtering_statistics$genes_after
    
    infoBox(
      "Genes (After QC)",
      format_number(n_genes),
      icon = icon("dna"),
      color = "purple"
    )
  })
  
  output$qc_box_retention <- renderInfoBox({
    req(rv$qc_report)
    
    retention_pct <- rv$qc_report$filtering_statistics$cells_retained_pct
    
    infoBox(
      "Cell Retention",
      sprintf("%.1f%%", retention_pct),
      icon = icon("percentage"),
      color = "yellow"
    )
  })
  
  # QC Metrics Table
  output$qc_metrics_table <- renderDT({
    req(rv$qc_report)
    
    # Extract metrics
    metrics_before <- rv$qc_report$qc_metrics_before
    metrics_after <- rv$qc_report$qc_metrics_after
    
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
  
  # QC Thresholds Text
  output$qc_thresholds_text <- renderText({
    req(rv$qc_report)
    
    thresholds <- rv$qc_report$thresholds_applied
    
    if (is.character(thresholds) && thresholds == "manual") {
      return("Manual thresholds were used (not MAD-based)")
    }
    
    # Format thresholds
    threshold_lines <- c("MAD-based automatic thresholds:\n")
    
    for (metric in names(thresholds)) {
      vals <- thresholds[[metric]]
      threshold_lines <- c(
        threshold_lines,
        sprintf("  %s: [%.1f, %.1f]", metric, vals[1], vals[2])
      )
    }
    
    paste(threshold_lines, collapse = "\n")
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
    
    datatable(
      rv$data_obj$metadata,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        search = list(search = ''),
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
  
  # Gene expression status
  output$gene_search_status <- renderText({
    "Enter a gene name and click 'Plot Expression' to visualize."
  })
  
  # Gene expression UMAP (placeholder for now)
  output$gene_expression_umap <- renderPlotly({
    req(input$btn_plot_gene)
    
    isolate({
      req(input$gene_search_input)
      req(rv$data_loaded)
      req(rv$data_obj$umap_coords)
      
      gene_name <- trimws(input$gene_search_input)
      
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
    })
  })
  
}
