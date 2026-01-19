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
                  value = "/srv/shiny-server/data/normalized_integrated.h5ad",
                  width = "100%"
                )
              ),
              column(
                width = 6,
                textInput(
                  "input_qc_dir",
                  "QC Results Directory:",
                  value = "/srv/shiny-server/data/qc_results",
                  width = "100%"
                )
              )
            ),
            
            checkboxInput(
              "input_backed_mode",
              "Use backed mode (recommended for >50k cells)",
              value = TRUE
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
            
            verbatimTextOutput("qc_thresholds_text")
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
            
            imageOutput("qc_plot_before", height = "500px")
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
            
            imageOutput("qc_plot_after", height = "500px")
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
              value = 3,
              step = 0.5
            ),
            
            sliderInput(
              "umap_opacity",
              "Opacity:",
              min = 0.1,
              max = 1,
              value = 0.7,
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
            title = "Gene Search",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            p("Enter a gene name to visualize its expression on UMAP:"),
            
            fluidRow(
              column(
                width = 8,
                textInput(
                  "gene_search_input",
                  "Gene Name:",
                  placeholder = "e.g., CD3D, CD79A, MS4A1",
                  width = "100%"
                )
              ),
              column(
                width = 4,
                br(),
                actionButton(
                  "btn_plot_gene",
                  "Plot Expression",
                  icon = icon("chart-line"),
                  class = "btn-primary"
                )
              )
            ),
            
            verbatimTextOutput("gene_search_status")
          )
        ),
        
        fluidRow(
          box(
            title = "Gene Expression UMAP",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            plotlyOutput("gene_expression_umap", height = "600px")
          )
        )
      ),
      
      # ========================================================================
      # TAB 5: About
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
              <h3>scAnnex: Interactive Downstream Analysis for scRNA-seq</h3>
              
              <p>
                <strong>scAnnex</strong> is a production-ready Nextflow pipeline with an 
                interactive Shiny dashboard for single-cell RNA-seq analysis.
              </p>
              
              <h4>Key Features:</h4>
              <ul>
                <li><strong>Unified Input:</strong> Supports H5AD, RDS (Seurat), and MTX formats</li>
                <li><strong>Quality Control:</strong> MAD-based automatic thresholding</li>
                <li><strong>Batch Integration:</strong> Harmony-based batch correction</li>
                <li><strong>Interactive Visualization:</strong> Real-time UMAP and gene expression plots</li>
                <li><strong>Scalable:</strong> Backed H5AD mode for datasets >100k cells</li>
              </ul>
              
              <h4>Technology Stack:</h4>
              <ul>
                <li><strong>Pipeline:</strong> Nextflow + Python (Scanpy)</li>
                <li><strong>Dashboard:</strong> R Shiny + reticulate</li>
                <li><strong>Containerization:</strong> Docker / Apptainer</li>
              </ul>
              
              <h4>Documentation:</h4>
              <p>
                For full documentation, see <code>InitProject.md</code> in the project repository.
              </p>
              
              <hr>
              
              <p><em>Version: 0.1.0 (Phase 8 - Dashboard Implementation)</em></p>
            ")
          )
        )
      )
    )
  )
)
