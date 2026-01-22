# scAnnex Dashboard - User Interface
# Defines the layout and UI components

ui <- dashboardPage(
  skin = "blue",
  
  # ==============================================================================
  # HEADER
  # ==============================================================================
  dashboardHeader(
    title = "scAnnex Dashboard",
    titleWidth = 300
  ),
  
  # ==============================================================================
  # SIDEBAR
  # ==============================================================================
  dashboardSidebar(
    width = 300,
    
    sidebarMenu(
      id = "sidebar_menu",
      
      menuItem(
        "Data Input",
        tabName = "tab_input",
        icon = icon("upload")
      ),
      
      menuItem(
        "QC Overview",
        tabName = "tab_qc",
        icon = icon("check-circle")
      ),
      
      menuItem(
        "Clustering & UMAP",
        tabName = "tab_clustering",
        icon = icon("project-diagram")
      ),
      
      menuItem(
        "Gene Expression",
        tabName = "tab_genes",
        icon = icon("dna")
      ),
      
      menuItem(
        "Annotation Station",
        tabName = "tab_annotation",
        icon = icon("tags")
      ),
      
      menuItem(
        "About",
        tabName = "tab_about",
        icon = icon("info-circle")
      )
    ),
    
    hr(),
    
    # Dataset info box
    div(
      style = "padding: 15px;",
      h5("Dataset Info", style = "color: white;"),
      verbatimTextOutput("sidebar_dataset_info", placeholder = TRUE)
    )
  ),
  
  # ==============================================================================
  # BODY
  # ==============================================================================
  dashboardBody(
    
    # Custom CSS
    tags$head(
      tags$style(HTML("
        .info-box { min-height: 90px; }
        .info-box-icon { height: 90px; line-height: 90px; }
        .info-box-content { padding-top: 10px; padding-bottom: 10px; }
        .box-title { font-size: 18px; font-weight: bold; }
        .shiny-notification { position: fixed; top: 50%; right: 50%; }
        
        /* QC Plot images - fit within box with max height */
        #qc_plot_before img, #qc_plot_after img {
          max-width: 100%;
          max-height: 500px;
          width: auto;
          height: auto;
          display: block;
          margin: 0 auto;
        }
      "))
    ),
    
    tabItems(
      
      # ========================================================================
      # TAB 1: Data Input
      # ========================================================================
      tabItem(
        tabName = "tab_input",
        
        fluidRow(
          box(
            title = "Load scRNA-seq Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            p("Select the H5AD file and QC results directory to visualize:"),
            
            fluidRow(
              column(
                width = 6,
                textInput(
                  "input_h5ad_path",
                  "H5AD File Path:",
                  value = DEFAULT_H5AD_FILE,
                  width = "100%"
                )
              ),
              column(
                width = 6,
                textInput(
                  "input_qc_dir",
                  "QC Results Directory:",
                  value = DEFAULT_QC_DIR,
                  width = "100%"
                )
              )
            ),
            
            checkboxInput(
              "input_backed_mode",
              "Use backed mode (recommended for >50k cells)",
              value = FALSE
            ),
            
            actionButton(
              "btn_load_data",
              "Load Data",
              icon = icon("play"),
              class = "btn-primary btn-lg"
            ),
            
            hr(),
            
            verbatimTextOutput("data_load_status")
          )
        )
      ),
      
      # ========================================================================
      # TAB 2: QC Overview
      # ========================================================================
      tabItem(
        tabName = "tab_qc",
        
        h2("Quality Control Overview"),
        
        # Summary boxes
        fluidRow(
          infoBoxOutput("qc_box_cells_before", width = 3),
          infoBoxOutput("qc_box_cells_after", width = 3),
          infoBoxOutput("qc_box_genes_after", width = 3),
          infoBoxOutput("qc_box_retention", width = 3)
        ),
        
        fluidRow(
          box(
            title = "QC Metrics Summary",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            
            DTOutput("qc_metrics_table")
          )
        ),
        
        fluidRow(
          box(
            title = "QC Thresholds Applied",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            DTOutput("qc_thresholds_table")
          )
        ),
        
        h3("QC Plots"),
        
        fluidRow(
          box(
            title = "Before Filtering",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            
            selectInput(
              "qc_plot_before_select",
              "Select Plot:",
              choices = c("Violin", "Scatter", "Distributions")
            ),
            
            imageOutput("qc_plot_before")
          ),
          
          box(
            title = "After Filtering",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            
            selectInput(
              "qc_plot_after_select",
              "Select Plot:",
              choices = c("Violin", "Scatter", "Distributions")
            ),
            
            imageOutput("qc_plot_after")
          )
        )
      ),
      
      # ========================================================================
      # TAB 3: Clustering & UMAP
      # ========================================================================
      tabItem(
        tabName = "tab_clustering",
        
        h2("Clustering & UMAP Visualization"),
        
        fluidRow(
          box(
            title = "UMAP Controls",
            status = "primary",
            solidHeader = TRUE,
            width = 3,
            
            selectInput(
              "umap_color_by",
              "Color by:",
              choices = c("batch", "sample_id", "condition"),
              selected = "batch"
            ),
            
            sliderInput(
              "umap_point_size",
              "Point size:",
              min = 1,
              max = 10,
              value = 5,
              step = 0.5
            ),
            
            sliderInput(
              "umap_opacity",
              "Opacity:",
              min = 0.1,
              max = 1,
              value = 1,
              step = 0.1
            )
          ),
          
          box(
            title = "Interactive UMAP",
            status = "info",
            solidHeader = TRUE,
            width = 9,
            
            plotlyOutput("umap_plot", height = "600px")
          )
        ),
        
        fluidRow(
          box(
            title = "Cell Metadata Table",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            DTOutput("metadata_table")
          )
        )
      ),
      
      # ========================================================================
      # TAB 4: Gene Expression
      # ========================================================================
      tabItem(
        tabName = "tab_genes",
        
        h2("Gene Expression Visualization"),
        
        fluidRow(
          box(
            title = "Gene Expression",
            status = "primary",
            solidHeader = TRUE,
            width = 3,
            
            textAreaInput(
              "gene_input",
              "Gene Name(s):",
              placeholder = "Single gene: CD3D\n\nMultiple genes (one per line):\nCD3D\nCD3E\nCD8A\nCD8B",
              rows = 8,
              width = "100%"
            ),
            
            actionButton(
              "btn_plot_genes",
              "Plot Expression",
              icon = icon("chart-line"),
              class = "btn-primary btn-block"
            ),
            
            br(),
            
            helpText(
              "Enter a single gene name for expression, or multiple genes (one per line) for gene set scoring."
            )
          ),
          
          box(
            title = "Gene Expression UMAP",
            status = "info",
            solidHeader = TRUE,
            width = 9,
            
            plotlyOutput("gene_expression_umap", height = "600px")
          )
        )
      ),
      
      # ========================================================================
      # TAB 5: Annotation Station
      # ========================================================================
      tabItem(
        tabName = "tab_annotation",
        
        h2("Annotation Station"),
        
        fluidRow(
          box(
            title = "Annotation Controls",
            status = "primary",
            solidHeader = TRUE,
            width = 3,
            
            textInput(
              "annot_name",
              "Annotation Name:",
              value = "custom_annotation",
              placeholder = "e.g., custom_first_annot"
            ),
            
            hr(),
            
            textAreaInput(
              "annot_rules",
              "Annotation Rules:",
              placeholder = "leiden_0.5,4,HSC\nleiden_0.5,2,T cells\nleiden_1.0,0,Monocytes",
              rows = 10,
              width = "100%"
            ),
            
            helpText(
              "Format: clustering_name,cluster_id,label",
              br(),
              "One rule per line. Later rules override earlier ones."
            ),
            
            hr(),
            
            h5("Visualization Settings"),
            
            selectInput(
              "annot_umap_select",
              "Select UMAP:",
              choices = c("X_umap"),
              selected = "X_umap"
            ),
            
            sliderInput(
              "annot_point_size",
              "Point size:",
              min = 1,
              max = 10,
              value = 5,
              step = 0.5
            ),
            
            sliderInput(
              "annot_opacity",
              "Opacity:",
              min = 0.1,
              max = 1,
              value = 1,
              step = 0.1
            ),
            
            hr(),
            
            actionButton(
              "btn_plot_annotation",
              "Plot",
              icon = icon("chart-area"),
              class = "btn-primary btn-block"
            ),
            
            br(),
            
            actionButton(
              "btn_save_annotation",
              "Save in Object",
              icon = icon("save"),
              class = "btn-success btn-block"
            ),
            
            br(),
            
            radioButtons(
              "annot_save_mode",
              "Save mode:",
              choices = c(
                "Overwrite original" = "overwrite",
                "Create new version" = "create_copy"
              ),
              selected = "create_copy"
            )
          ),
          
          box(
            title = "Custom Annotation UMAP",
            status = "info",
            solidHeader = TRUE,
            width = 9,
            
            plotlyOutput("annotation_umap", height = "600px"),
            
            br(),
            
            verbatimTextOutput("annotation_stats")
          )
        )
      ),
      
      # ========================================================================
      # TAB 6: About
      # ========================================================================
      tabItem(
        tabName = "tab_about",
        
        h2("About scAnnex Dashboard"),
        
        fluidRow(
          box(
            title = "Project Information",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            HTML("
              <h3>scAnnex Dashboard</h3>
              
              <p>
                Interactive visualization and exploration dashboard for single-cell RNA-seq data
                processed through the <strong>scAnnex</strong> Nextflow pipeline.
              </p>
              
              <h4>Dashboard Features:</h4>
              <ul>
                <li><strong>Quality Control Overview:</strong> Interactive QC metrics, thresholds tables, and before/after filtering plots</li>
                <li><strong>UMAP Visualization:</strong> Interactive exploration with customizable coloring by metadata (batch, sample, condition, etc.)</li>
                <li><strong>Automatic Cell-Type Annotation:</strong> Visualize CellTypist auto-annotations in the Clustering & UMAP tab</li>
                <li><strong>Gene Expression:</strong> Single gene expression visualization on UMAP</li>
                <li><strong>Gene Set Scoring:</strong> Calculate and visualize gene signature scores (0-1 normalized scale)</li>
                <li><strong>Annotation Station:</strong> Create custom cell-type annotations by assigning labels to clusters, with real-time visualization and H5AD export</li>
                <li><strong>Metadata Export:</strong> Filter and download cell metadata as CSV or Excel</li>
                <li><strong>Auto-Detection:</strong> Automatically finds H5AD and QC files in results directory</li>
              </ul>
              
              <h4>Pipeline Features:</h4>
              <ul>
                <li><strong>Unified Input:</strong> H5AD, RDS (Seurat), and MTX formats</li>
                <li><strong>Quality Control:</strong> MAD-based automatic thresholding with detailed attrition tracking</li>
                <li><strong>Doublet Detection:</strong> Scrublet integration for doublet removal</li>
                <li><strong>Normalization:</strong> Log-normalization and highly variable gene detection</li>
                <li><strong>Batch Integration:</strong> Harmony-based batch correction</li>
                <li><strong>Dimensionality Reduction:</strong> PCA and UMAP generation</li>
                <li><strong>Clustering:</strong> Leiden clustering with multiple resolutions</li>
                <li><strong>Cell Type Annotation:</strong> Automated annotation using CellTypist</li>
              </ul>
              
              <h4>Technology:</h4>
              <ul>
                <li><strong>Pipeline:</strong> Nextflow DSL2 + Python (Scanpy/AnnData)</li>
                <li><strong>Dashboard:</strong> R Shiny + reticulate + plotly</li>
                <li><strong>Execution:</strong> Conda environments or containerization (Docker/Apptainer)</li>
              </ul>
              
              <h4>Documentation & Source Code:</h4>
              <p>
                <a href='https://github.com/damouzo/scAnnex' target='_blank' style='font-size: 16px;'>
                  <i class='fa fa-github'></i> GitHub Repository: damouzo/scAnnex
                </a>
              </p>
              <p>
                Visit the repository for complete documentation, installation instructions, usage examples, and the latest updates.
              </p>
              
              <hr>
              
              <p style='color: #666;'>
                <em>Dashboard Version: 1.0.0 | Pipeline: scAnnex Nextflow</em>
              </p>
            ")
          )
        )
      )
    )
  )
)
