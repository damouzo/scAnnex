# scAnnex Dashboard - User Interface (bslib page_navbar)

suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(DT)
    library(plotly)
})

# Helper: compact KPI card used in Overview tab
sa_kpi <- function(icon_name, label, value_id, bg, ns = identity) {
    tags$div(
        class = "sa-kpi",
        style = paste0("background-color:", bg, ";"),
        tags$div(class = "sa-kpi-icon", icon(icon_name)),
        tags$div(
            class = "sa-kpi-body",
            tags$div(class = "sa-kpi-label", label),
            tags$div(class = "sa-kpi-value",
                     textOutput(ns(value_id), inline = TRUE))
        )
    )
}

ui <- page_navbar(
    title = tags$span("scAnnex", style = "font-weight: 700; letter-spacing: 0.03em;"),
    theme = bs_theme(
        bootswatch  = "flatly",
        primary     = "#1565C0",
        secondary   = "#00897B",
        info        = "#0288D1",
        base_font   = font_google("Inter"),
        code_font   = font_google("Source Code Pro")
    ),
    header = tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "scannex.css")
    ),
    window_title = "scAnnex Dashboard",
    id = "main_navbar",

    # =========================================================================
    # Tab 1: Overview
    # =========================================================================
    nav_panel(
        title = "Overview",
        icon  = icon("home"),

        fluidRow(
            column(12,
                tags$div(
                    style = "margin-bottom: 8px;",
                    h4("scRNA-seq Analysis Overview",
                       style = "margin: 0; font-weight: 600; color: #1565C0;"),
                    tags$small(class = "text-muted",
                               textOutput("overview_results_dir_label", inline = TRUE))
                )
            )
        ),

        fluidRow(
            column(2, sa_kpi("layer-group",  "Samples",           "ov_n_samples",     "#1565C0")),
            column(2, sa_kpi("circle",       "Cells (before QC)", "ov_cells_before",  "#1976D2")),
            column(2, sa_kpi("check-circle", "Cells (after QC)",  "ov_cells_after",   "#00897B")),
            column(2, sa_kpi("percent",      "Avg Retention",     "ov_avg_retention", "#0288D1")),
            column(2, sa_kpi("dna",          "Avg Genes / Cell",  "ov_avg_genes",     "#5C6BC0")),
            column(2, sa_kpi("biohazard",    "Avg Median MT%",    "ov_avg_mt",        "#7B1FA2"))
        ),

        fluidRow(
            column(12,
                card(
                    card_header("Per-Sample QC Summary"),
                    DTOutput("overview_sample_table")
                )
            )
        ),

        fluidRow(
            column(12,
                card(
                    card_header(
                        tags$span(
                            "QC Metric Distributions by Sample",
                            tags$small(class = "text-muted ms-2",
                                       "(after QC filtering, one line per sample)")
                        )
                    ),
                    uiOutput("overview_density_status"),
                    fluidRow(
                        column(4, plotOutput("overview_density_mt",     height = "260px")),
                        column(4, plotOutput("overview_density_genes",  height = "260px")),
                        column(4, plotOutput("overview_density_counts", height = "260px"))
                    )
                )
            )
        )
    ),

    # =========================================================================
    # Tab 2: Quality Control Overview by Sample
    # =========================================================================
    nav_panel(
        title = tags$span(
            tags$span(class = "nav-label-long", "Quality Control Overview by Sample"),
            tags$span(class = "nav-label-short", "QC by Sample")
        ),
        icon = icon("check-circle"),

        fluidRow(column(12, h4("Quality Control Overview by Sample",
                               style = "font-weight: 600; color: #1565C0;"))),

        conditionalPanel(
            condition = "output.is_integrated",
            fluidRow(
                column(12,
                    card(
                        card_header("Sample Selector"),
                        selectInput("qc_sample_select", "View QC for Sample:",
                                    choices = c("summary"), selected = "summary",
                                    width = "40%")
                    )
                )
            )
        ),

        fluidRow(
            column(3, uiOutput("qc_box_cells_before")),
            column(3, uiOutput("qc_box_cells_after")),
            column(3, uiOutput("qc_box_genes_after")),
            column(3, uiOutput("qc_box_retention"))
        ),

        fluidRow(
            column(12,
                card(card_header("QC Metrics Summary"), DTOutput("qc_metrics_table"))
            )
        ),

        fluidRow(
            column(12,
                card(card_header("QC Thresholds Applied"), DTOutput("qc_thresholds_table"))
            )
        ),

        h5("QC Plots", style = "font-weight: 600; margin-top: 8px;"),

        fluidRow(
            column(6,
                card(
                    card_header("Before Filtering"),
                    selectInput("qc_plot_before_select", "Select Plot:",
                                choices = c("Violin", "Scatter", "Distributions")),
                    imageOutput("qc_plot_before", height = "520px")
                )
            ),
            column(6,
                card(
                    card_header("After Filtering"),
                    selectInput("qc_plot_after_select", "Select Plot:",
                                choices = c("Violin", "Scatter", "Distributions")),
                    imageOutput("qc_plot_after", height = "520px")
                )
            )
        )
    ),

    # =========================================================================
    # Tab 3: Clustering & UMAP
    # =========================================================================
    nav_panel(
        title = "Clustering & UMAP",
        icon  = icon("project-diagram"),

        fluidRow(column(12, h4("Clustering & UMAP Visualization",
                               style = "font-weight: 600; color: #1565C0;"))),

        fluidRow(
            column(3,
                card(
                    card_header("UMAP Controls"),
                    selectInput("umap_color_by", "Color by:",
                                choices = c("batch", "sample_id", "condition"),
                                selected = "batch"),
                    tags$div(
                        style = "display: flex; gap: 10px;",
                        tags$div(style = "flex:1;",
                                 numericInput("umap_point_size", "Point size:",
                                              value = 3, min = 0.5, max = 20, step = 0.5)),
                        tags$div(style = "flex:1;",
                                 numericInput("umap_opacity", "Opacity:",
                                              value = 0.8, min = 0.1, max = 1, step = 0.05))
                    ),
                    helpText("Legend marker size is fixed (8) regardless of point size.")
                )
            ),
            column(9,
                card(
                    card_header("Interactive UMAP"),
                    plotlyOutput("umap_plot", height = "600px")
                )
            )
        ),

        fluidRow(
            column(12,
                card(card_header("Cell Metadata Table"), DTOutput("metadata_table"))
            )
        )
    ),

    # =========================================================================
    # Tab 4: Gene Expression
    # =========================================================================
    nav_panel(
        title = "Gene Expression",
        icon  = icon("dna"),

        fluidRow(column(12, h4("Gene Expression Visualization",
                               style = "font-weight: 600; color: #1565C0;"))),

        fluidRow(
            column(3,
                card(
                    card_header("Settings"),
                    textAreaInput(
                        "gene_input", "Gene Name(s):",
                        placeholder = "Single gene: DDX41\n\nMultiple genes (one per line):\nDDX41\nDHX34\nRPS27",
                        rows = 8, width = "100%"
                    ),
                    actionButton("btn_plot_genes", "Plot Expression",
                                 icon = icon("chart-line"), class = "btn-primary w-100"),
                    helpText("Single gene: expression gradient. Multiple genes: gene set score (0-1)."),
                    tags$hr(),
                    tags$div(
                        class = "top-n-section",
                        tags$strong("Highlight Top Cells", style = "font-size: 0.88rem;"),
                        tags$p(class = "text-muted",
                               style = "font-size: 0.78rem; margin: 3px 0 6px 0;",
                               "Mark top N expressing cells in red; all others in gray. Set 0 for normal gradient."),
                        numericInput("gene_top_n_cells", "Top N cells to highlight:",
                                     value = 0, min = 0, max = 5000, step = 10)
                    )
                )
            ),
            column(9,
                card(
                    card_header("Expression UMAP"),
                    plotlyOutput("gene_expression_umap", height = "600px")
                )
            )
        )
    ),

    # =========================================================================
    # Tab 5: Differential Expression
    # =========================================================================
    nav_panel(
        title = tags$span(
            tags$span(class = "nav-label-long", "Differential Expression"),
            tags$span(class = "nav-label-short", "DGE")
        ),
        icon = icon("chart-line"),

        fluidRow(column(12, h4("Differential Gene Expression Analysis",
                               style = "font-weight: 600; color: #1565C0;"))),

        fluidRow(
            column(3,
                card(
                    card_header("Settings"),
                    selectInput("dge_contrast_select", "Select Contrast:",
                                choices = c("No contrasts loaded"), selected = NULL),
                    tags$hr(),
                    tags$h6("Filtering Options", style = "font-weight: 600; margin-bottom: 8px;"),
                    numericInput("dge_pval_threshold", "Adj. P-value threshold:",
                                 value = 0.05, min = 1e-10, max = 1, step = 0.001),
                    numericInput("dge_logfc_threshold", "Log2 FC threshold:",
                                 value = 0.25, min = 0, max = 10, step = 0.05),
                    numericInput("dge_top_n_genes", "Top N genes to label:",
                                 value = 10, min = 0, max = 200, step = 1),
                    sliderInput("dge_gene_label_size", "Label font size:",
                                min = 2, max = 8, value = 3, step = 0.5),
                    helpText("Labels balance log2FC extremes and significance equally."),
                    tags$hr(),
                    actionButton("btn_apply_dge", "Apply",
                                 icon = icon("play"), class = "btn-primary w-100"),
                    br(), br(),
                    downloadButton("btn_download_volcano", "Download Volcano (PNG)",
                                   class = "btn-outline-secondary w-100 btn-sm")
                )
            ),
            column(9,
                navset_card_tab(
                    nav_panel("Volcano Plot",
                              plotOutput("dge_volcano_plot", height = "600px")),
                    nav_panel("Significant Genes",
                              DTOutput("dge_significant_genes_table"),
                              br(),
                              fluidRow(
                                  column(5, downloadButton("btn_download_sig_genes",
                                                           "Significant Genes (CSV)",
                                                           class = "btn-sm btn-outline-secondary")),
                                  column(5, downloadButton("btn_download_all_results",
                                                           "All Results (CSV)",
                                                           class = "btn-sm btn-outline-secondary"))
                              ))
                )
            )
        )
    ),

    # =========================================================================
    # Tab 6: GSEA
    # =========================================================================
    nav_panel(
        title = "GSEA",
        icon  = icon("route"),

        fluidRow(column(12, h4("Gene Set Enrichment Analysis",
                               style = "font-weight: 600; color: #1565C0;"))),

        fluidRow(
            column(3,
                card(
                    card_header("Settings"),
                    selectInput("gsea_contrast_select", "Select Contrast:",
                                choices = c("No contrasts loaded"), selected = NULL),
                    selectInput("gsea_db_select", "Gene Set Database:",
                                choices = c(
                                    "GO Biological Process" = "go_bp",
                                    "GO Molecular Function" = "go_mf",
                                    "GO Cellular Component" = "go_cc",
                                    "KEGG"                  = "kegg",
                                    "Reactome"              = "reactome"
                                ),
                                selected = "go_bp"),
                    numericInput("gsea_n_pathways", "Top N pathways:",
                                 value = 10, min = 1, max = 50, step = 1),
                    numericInput("gsea_n_gseaplot", "GSEA plot: top N:",
                                 value = 5, min = 1, max = 20, step = 1),
                    numericInput("gsea_padj_cutoff", "Max adj. p-value:",
                                 value = 0.25, min = 0.001, max = 1, step = 0.01),
                    actionButton("btn_apply_gsea", "Apply",
                                 icon = icon("play"), class = "btn-primary w-100 mt-2")
                )
            ),
            column(9,
                uiOutput("gsea_main_panels")
            )
        )
    ),

    # =========================================================================
    # Tab 7: Annotation Station
    # =========================================================================
    nav_panel(
        title = tags$span(
            tags$span(class = "nav-label-long", "Annotation Station"),
            tags$span(class = "nav-label-short", "Annot.")
        ),
        icon = icon("tags"),

        fluidRow(column(12, h4("Annotation Station",
                               style = "font-weight: 600; color: #1565C0;"))),

        fluidRow(
            column(3,
                card(
                    card_header("Annotation Controls"),
                    textInput("annot_name", "Annotation Name:",
                              value = "custom_annotation",
                              placeholder = "e.g., custom_first_annot"),
                    tags$hr(),
                    textAreaInput("annot_rules", "Annotation Rules:",
                                  placeholder = "leiden_0.5,4,HSC\nleiden_0.5,2,T cells\nleiden_1.0,0,Monocytes",
                                  rows = 10, width = "100%"),
                    tags$p(class = "annotation-help",
                           "Format: clustering_name,cluster_id,label", br(),
                           "One rule per line. Later rules override earlier ones."),
                    tags$hr(),
                    h6("Visualization Settings", style = "font-weight: 600;"),
                    selectInput("annot_umap_select", "Select UMAP:",
                                choices = c("X_umap"), selected = "X_umap"),
                    tags$div(
                        style = "display: flex; gap: 10px;",
                        tags$div(style = "flex:1;",
                                 numericInput("annot_point_size", "Point size:",
                                              value = 3, min = 0.5, max = 20, step = 0.5)),
                        tags$div(style = "flex:1;",
                                 numericInput("annot_opacity", "Opacity:",
                                              value = 0.8, min = 0.1, max = 1, step = 0.05))
                    ),
                    tags$hr(),
                    actionButton("btn_plot_annotation", "Plot",
                                 icon = icon("chart-area"), class = "btn-primary w-100"),
                    br(),
                    actionButton("btn_save_annotation", "Save in Object",
                                 icon = icon("save"), class = "btn-success w-100"),
                    br(),
                    radioButtons("annot_save_mode", "Save mode:",
                                 choices = c("Overwrite original" = "overwrite",
                                             "Create new version" = "create_copy"),
                                 selected = "create_copy")
                )
            ),
            column(9,
                card(
                    card_header("Custom Annotation UMAP"),
                    plotlyOutput("annotation_umap", height = "560px"),
                    br(),
                    verbatimTextOutput("annotation_stats")
                )
            )
        )
    ),

    # =========================================================================
    # Tab 8: About
    # =========================================================================
    nav_panel(
        title = "About",
        icon  = icon("info-circle"),

        fluidRow(
            column(12,
                card(
                    card_header("About scAnnex"),
                    HTML(paste0(
                        "<h4>scAnnex Dashboard</h4>",
                        "<p>Interactive visualization and exploration dashboard for single-cell RNA-seq data ",
                        "processed through the <strong>scAnnex</strong> Nextflow DSL2 pipeline.</p>",
                        "<h5>Dashboard Features</h5>",
                        "<ul>",
                        "<li><strong>Overview:</strong> Multi-sample QC summary, KPI cards, per-sample statistics, density plots</li>",
                        "<li><strong>Quality Control by Sample:</strong> Per-sample QC metrics, thresholds, and before/after filtering plots</li>",
                        "<li><strong>Clustering &amp; UMAP:</strong> Interactive UMAP; numeric point size/opacity; fixed-size legend markers</li>",
                        "<li><strong>Gene Expression:</strong> Single gene or gene set scoring; top-N cell highlighting</li>",
                        "<li><strong>Differential Expression:</strong> Volcano plots with balanced LFC+significance labeling</li>",
                        "<li><strong>GSEA:</strong> Dotplot, Ridgeplot, GSEA running score, Pathview (KEGG), TreeDot tabs</li>",
                        "<li><strong>Annotation Station:</strong> Custom cluster-to-label annotation with H5AD export</li>",
                        "</ul>",
                        "<hr>",
                        "<p class='text-muted'><em>scAnnex v0.1.0 | Nextflow DSL2 + Python (Scanpy) + R Shiny + bslib</em></p>"
                    ))
                )
            )
        )
    ),

    nav_spacer(),
    nav_item(
        tags$small(class = "text-muted",
                   paste0("scAnnex v0.1.0 | ", format(Sys.Date(), "%Y")))
    )
)
