
Project: scAnnex
Subtitle: Interactive Downstream Analysis Automation for Single-Cell RNA-seq

Note on inspiration and structuring:
The idea, especially regarding the input module and the use of Python, is inspired by https://deepwiki.com/knollr/scRNAseq_pipeline_rapids_dashboard.
For the overall structuring, the goal is for scAnnex to be the natural next step after https://github.com/nf-core/scrnaseq. It is not the same tool, but it should be usable immediately after, following the same logic of structuring and best practices for Nextflow workflows.

Status: Architecture Design Phase (v1.0)

Primary Goal: To provide a seamless, containerized bridge between raw mapping outputs (e.g., nf-core/scrnaseq) and expert biological curation via an interactive dashboard.

1. System Architecture

The tool is split into two distinct layers to balance high-performance computing with high-flexibility visualization.

   - Orchestration Layer: Nextflow (DSL2). Manages data flow, resource allocation, and reproducibility.
   - Processing Engine (The "Body"): Python (Scanpy). Handles heavy computation, data unification, and batch correction.
   - Visualization Layer (The "Soul"): R (Shiny). Provides an interactive interface for biologists to perform manual curation and biomarker discovery.
   - Environment Management: Apptainer (Singularity). Every process must run within a predefined container to ensure portability across HPC environments.

**Critical Architectural Decision (Seqera Recommendation):**
Separate the Nextflow pipeline from the Shiny dashboard. Nextflow should produce outputs ready for consumption, and Shiny should be an independent application that reads those outputs. Do not attempt to launch Shiny from within Nextflow to avoid maintenance issues.

Recommended project structure:
```
scAnnex/
├── workflows/        # Nextflow pipeline
├── modules/          # Modular processes
├── subworkflows/     # Reusable subworkflows
├── bin/              # Executable Python/R scripts
├── dashboard/        # Independent Shiny app
├── containers/       # Container definitions
└── conf/             # Configurations
```

2. Pipeline Specifications (Nextflow)

2.1 Universal Input Module

The pipeline must accept a samplesheet.csv containing paths to three potential formats:
    - .h5ad: AnnData objects.
    - .rds: Seurat objects (converted via anndata2ri or reticulate).
    - 10x MTX: Directory containing matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz.

**Recommended samplesheet structure:**
```csv
sample_id,file_type,file_path,batch,condition
PBMC_1,h5ad,/path/to/sample1.h5ad,batch1,control
PBMC_2,rds,/path/to/sample2.rds,batch1,treated
PBMC_3,mtx,/path/to/sample3/,batch2,control
```

**Critical Input Validation:**
- Use nf-validation plugin (nf-core standard) for samplesheet validation
- Avoid absolute paths in samplesheet; use paths relative to work directory or S3 URIs
- Implement validateParameters() at workflow start

**RDS → H5AD Conversion Considerations:**
- anndata2ri has known limitations: may lose Seurat metadata (custom reductions), R factor to pandas category conversion can fail, and Seurat v5 objects (with layers) need special handling.
- **Recommended approach:** Use SeuratDisk + Python
  1. In R: SaveH5Seurat() → .h5seurat file
  2. In Python: Convert .h5seurat → .h5ad with full control

2.2 Processing Workflow (Modular Steps)

    Step 1: Unification: Standardize all inputs into a single AnnData structure.
    Step 2: Quality Control (QC): Calculate mitochondrial, ribosomal, and hemoglobin gene percentages.
    Step 3: Doublet Detection (Extra): Implement Scrublet to score and remove potential doublets.
        - Note: Scrublet is slow on datasets >50k cells. Consider alternatives: DoubletFinder or scDblFinder.
    Step 4: Normalization & Integration: Standard Log-normalization.
        - Harmony Integration: Use scanpy.external.pp.harmony (faster than harmonypy wrapper) if multiple samples/batches are detected.
        - Always generate pre/post integration plots and calculate quality metrics (kBET, LISI).
        - Save pre-integration versions for potential rollback.
        - Validate that batches have sufficient variation (Harmony can fail silently otherwise).
    Step 5: Dimensionality Reduction: Compute PCA and UMAP embeddings.
    Step 6: Auto-Annotation: Preliminary cluster naming based on reference marker lists (CSV).
        - Always mark as "preliminary" in results to avoid misleading interpretations.

**Process Structure Best Practices:**
Each Nextflow process should follow this template:
```groovy
process UNIFY_INPUT {
    tag "$meta.id"  // For log identification
    label 'process_low'  // Defined in conf/base.config
    container 'scanpy_container:1.0'  // Or use Wave
    
    input:
    tuple val(meta), path(input_file)
    
    output:
    tuple val(meta), path("${meta.id}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions
    
    script:
    def file_type = meta.file_type
    """
    unify_input.py \\
        --input ${input_file} \\
        --type ${file_type} \\
        --output ${meta.id}.h5ad
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """
}
```

**Version Emission (nf-core standard):**
Every process must emit versions.yml for reproducibility.

3. Dashboard Specifications (Shiny)

The dashboard must be optimized for speed. It should load a "Lightweight H5AD" (or a backed object) to prevent high RAM consumption on the server.

**Performance Architecture:**
Do not load the complete H5AD at Shiny startup. Instead, load only metadata and UMAP coordinates:
```r
# At app initialization:
umap_coords <- read.csv("umap_coordinates.csv")
metadata <- read.csv("cell_metadata.csv")
markers_list <- readRDS("marker_genes.rds")

# Load on demand when user requests FeaturePlot for gene X:
get_gene_expression <- function(gene_name) {
    # Read only that column from backed H5AD
    # Or use a pre-built SQLite database
}
```

**Large Dataset Handling (>100k cells):**
```python
# In Python processing:
import scanpy as sc

# For writing:
adata.write_h5ad('output.h5ad', compression='gzip')

# For reading in dashboard (BACKED MODE):
adata = sc.read_h5ad('output.h5ad', backed='r')
# Only reads metadata, matrices stay on disk
```

```r
# For Shiny with reticulate + anndata backed:
library(reticulate)
ad <- import("anndata")
adata <- ad$read_h5ad("output.h5ad", backed = "r")
```

**Alternative:** Export specific matrices from Nextflow:
- UMAP coordinates → lightweight CSV
- Gene expression subsets → small dense matrices
- Metadata → parquet or CSV

Tab 1: QC Overview
    - Visualizations: Violin plots and Scatter plots (nFeature vs nCount).
    - Interactive: Sliders to view how different thresholds would affect the cell count.

Tab 2: Clustering & Auto-Annot
    - Visualization: Interactive UMAP using plotly (recommended):
```r
library(plotly)
plot_ly(
    data = umap_coords,
    x = ~UMAP_1,
    y = ~UMAP_2,
    color = ~cluster,
    type = 'scatter',
    mode = 'markers',
    text = ~paste("Cell:", cell_id)
)
```
    - Feature: Toggle between "Leiden Clusters" and "Automatic Cell Type Labels."

Tab 3: Manual Curation & Biomarkers
    - Search Engine: User input for gene names. Returns a FeaturePlot (UMAP colored by expression).
    - AUCell Integration: Capability to upload a gene set (CSV/GMT) to calculate and visualize signature scores across clusters.
        - **CRITICAL:** Implement `withProgress()` in Shiny for custom signature calculations so users know the app is working and not frozen.
        - Provide default gene sets (.gmt or .csv) with basic signatures (Immune, Cell Cycle, Stress) so Tab 3 is not empty at startup.
    - Renaming Module: Allow users to manually rename clusters and export the final annotated object.

**Result Export Capabilities:**
Define what users can download from the Dashboard:
- **Essential:** Final .h5ad object with manual curation saved (most valuable output)
- UMAP visualizations as high-resolution PNG/PDF
- CSV with cluster renaming table (old → new labels)
- Cell metadata table with final annotations
- Marker gene tables per cluster
- QC summary report (HTML or PDF)

**Deployment Considerations:**
- Shiny consumes RAM even with backed mode if there are many users
- Consider ShinyProxy or Posit Connect for production deployments

4. Development Roadmap (PoC)

Phase    Milestone            Deliverables
Phase 1  The Input Bridge     Python script for RDS/H5AD/MTX conversion.
Phase 2  The Core Pipeline    Nextflow workflow with Scanpy + Harmony.
Phase 3  Interactive PoC      Basic Shiny app reading processed H5AD.
Phase 4  Deployment           Apptainer .sif files and HPC testing.

5. Technical Requirements & Best Practices

**Container Management (CRITICAL):**
Use profiles in nextflow.config for maximum flexibility across environments:
```groovy
profiles {
    docker {
        docker.enabled = true
        singularity.enabled = false
    }
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        docker.enabled = false
    }
    apptainer {
        apptainer.enabled = true
        apptainer.autoMounts = true
    }
    test {
        params.input = 'test_data/samplesheet_test.csv'
        params.outdir = 'test_results'
        params.max_memory = '6.GB'
        params.max_cpus = 2
    }
}
```

**Container Best Practices:**
- Use Wave (Seqera tool) for building containers on-demand with conda packages
- Specify EXACT package versions: scanpy==1.9.8, harmonypy==0.0.10
- Test containers locally with Docker before converting to Singularity/Apptainer

Version Control: Git-flow. Use feature branches.
Nextflow Config: Use a nextflow.config with profiles for local and apptainer.
Resource Management: Assign specific cpus and memory labels to processes in conf/base.config.
Data Handling: Use H5AD "backed" mode for matrices > 100k cells to ensure the Shiny app remains responsive.
Documentation:
    - Python: Google-style Docstrings.
    - R: Roxygen2 comments.
    - Project: Concise README.md and a progress.todo file in the root directory.

**Testing Strategy:**
Create tests for each module:
```groovy
// In tests/test_unify.nf
nextflow.enable.dsl=2

include { UNIFY_INPUT } from '../modules/unify_input.nf'

workflow {
    Channel
        .of([
            [id: 'test', file_type: 'h5ad'], 
            file('test_data/sample.h5ad')
        ])
        .set { ch_input }
    
    UNIFY_INPUT(ch_input)
}
```

Use pytest for Python scripts and testthat for R scripts.

6. Demo Dataset for Testing

Use the 10x Genomics 1k PBMC v3 dataset (free and standard):
```bash
# Download the dataset:
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5

# For integration testing, use:
# - PBMC 1k (batch1)
# - PBMC 1k artificial subset (batch2 with simulated batch effect)
```

Test samples:
    - Sample A: Processed as .h5ad.
    - Sample B: Processed as .rds. Purpose: Validate the unification module and test Harmony's batch correction capabilities.

Next Step for the Team: Initialize the GitHub repository and implement the UNIFY_INPUT module in Python. Refer to the knollr repository for initial input logic but adapt it to the scAnnex modular structure.

---

## 7. Pre-Development Checklist (Seqera Critical Recommendations)

Before writing any code, ensure the following are defined:

- [ ] Define schema.json for parameter validation (nf-validation)
- [ ] Choose container strategy (Wave vs custom Dockerfiles)
- [ ] Decide final output format (single H5AD? per-sample?)
- [ ] Plan AnnData metadata structure (.obs, .var, .uns conventions)
- [ ] Define cluster naming convention (leiden_res_0.5, etc.)
- [ ] Establish logging strategy (MultiQC? Custom report?)
- [ ] Define testing strategy (pytest for Python, testthat for R)
- [ ] Prepare default gene set collections (.gmt/.csv) for AUCell (Immune signatures, Cell Cycle, Stress response)
- [ ] Specify export formats and functionality for Dashboard (H5AD, plots, tables)
- [ ] Design progress indicators for long-running Shiny operations (withProgress())

## 8. Common Pitfalls to Avoid (Critical Warnings)

⚠️ **Do NOT:**
- Use absolute paths in samplesheet → use relative to work dir or S3 URIs
- Assume Harmony will work without validation → it can fail silently if no batch variation exists
- Rely solely on Scrublet for large datasets (>50k cells) → consider alternatives
- Present auto-annotation as definitive → always mark as "preliminary"
- Launch Shiny from within Nextflow → keep them separate
- Load complete H5AD in Shiny without backed mode → RAM issues with multiple users
- Ignore version tracking → every process must emit versions.yml

## 9. Metadata and Output Conventions

**AnnData Structure Standards:**
- `.obs`: Cell-level metadata (clusters, QC metrics, batch info)
- `.var`: Gene-level metadata (highly variable genes flags)
- `.uns`: Unstructured metadata (parameters, plots, color maps)
- `.obsm`: Multi-dimensional annotations (PCA, UMAP coordinates)
- `.layers`: Alternative data representations (raw counts, normalized)

**Cluster Naming:**
Use consistent naming convention: `leiden_res_0.5`, `leiden_res_1.0`, etc.

**Output Files to Generate:**
- Primary: `{sample_id}.h5ad` (full AnnData object with backed support)
- For Dashboard:
  - `umap_coordinates.csv` (lightweight, fast loading)
  - `cell_metadata.csv` (sample, cluster, QC metrics)
  - `marker_genes.rds` or `.csv` (top markers per cluster)
  - `qc_plots/` (pre-generated quality control figures)
  - `gene_sets/` (default signature collections: immune.gmt, cell_cycle.gmt, stress.gmt)
- User-exportable:
  - `{sample_id}_curated.h5ad` (with manual annotations)
  - `cluster_annotations.csv` (mapping table)
  - `umap_final.pdf/png` (publication-ready figures)

## 10. Integration Quality Metrics

After running Harmony or any batch correction:

**Required Validations:**
1. Generate pre/post integration UMAP plots
2. Calculate batch mixing metrics:
   - kBET (k-nearest neighbor batch effect test)
   - LISI (Local Inverse Simpson's Index)
3. Visual inspection of known marker genes
4. Check if biological signal is preserved (not over-corrected)

**Save checkpoint:** Always keep pre-integration AnnData for rollback if needed.