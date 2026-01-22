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
    
    n_cells <- rv$qc_report$filtering_statistics$cells_initial
    
    infoBox(
      "Cells (Before QC)",
      format_number(n_cells),
      icon = icon("circle"),
      color = "blue"
    )
  })
  
  output$qc_box_cells_after <- renderInfoBox({
    req(rv$qc_report)
    
    n_cells <- rv$qc_report$filtering_statistics$cells_final
    
    infoBox(
      "Cells (After QC)",
      format_number(n_cells),
      icon = icon("check-circle"),
      color = "green"
    )
  })
  
  output$qc_box_genes_after <- renderInfoBox({
    req(rv$qc_report)
    
    n_genes <- rv$qc_report$filtering_statistics$genes_final
    
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
  
  # QC Thresholds Table
  output$qc_thresholds_table <- renderDT({
    req(rv$qc_report)
    
    thresholds <- rv$qc_report$thresholds_applied
    
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
  # TAB 5: ANNOTATION STATION
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
