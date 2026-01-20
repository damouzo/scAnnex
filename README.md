# scAnnex

A production-ready Nextflow pipeline for automated single-cell RNA-seq analysis with integrated quality control, batch correction, and interactive visualization.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Overview

scAnnex implements a comprehensive workflow for downstream analysis of single-cell RNA-seq data, providing automated quality control with adaptive thresholds, batch integration using Harmony, and an interactive Shiny dashboard for data exploration. The pipeline is designed for reproducibility and scalability, supporting execution on local machines, HPC clusters, and cloud environments.

### Key Features

- **Flexible Input Formats**: Native support for H5AD (AnnData), RDS (Seurat), and 10X MTX formats
- **Adaptive Quality Control**: MAD-based automatic threshold calculation with manual override capability
- **Batch Integration**: Harmony-based correction with pre/post integration quality assessment
- **Polyglot Annotation**: Seamless integration of Python (CellTypist) and R (Azimuth) annotation tools
- **Interactive Dashboard**: Docker-based Shiny application with backed H5AD reading for large datasets
- **Production-Grade**: Comprehensive logging, reproducible environments, and standardized outputs

## Installation

### Prerequisites

1. **Nextflow** (version 23.04.0 or later)
   ```bash
   curl -s https://get.nextflow.io | bash
   mv nextflow /usr/local/bin/
   ```

2. **Container Engine** (choose one):
   - Docker (recommended for local execution)
   - Singularity/Apptainer (recommended for HPC)
   - Podman (alternative to Docker)

### Quick Test

Verify installation with the included test profile:

```bash
nextflow run main.nf -profile test,docker
```

This runs a minimal dataset (1,049 cells) through the complete pipeline in approximately 5-10 minutes.

## Quick Start

### Basic Analysis

Process a single H5AD file with default QC parameters:

```bash
nextflow run main.nf \
  --input sample_data.h5ad \
  --input_type h5ad \
  --outdir results \
  -profile docker
```

### Analysis with Batch Correction

Integrate multiple batches using Harmony:

```bash
nextflow run main.nf \
  --input sample_data.h5ad \
  --input_type h5ad \
  --run_integration true \
  --batch_key batch \
  --integration_method harmony \
  --outdir results \
  -profile docker
```

### With Automated Annotation

Enable CellTypist annotation:

```bash
nextflow run main.nf \
  --input sample_data.h5ad \
  --annotation.enabled true \
  --annotation.tools celltypist \
  --annotation.celltypist_model Immune_All_Low.pkl \
  --outdir results \
  -profile docker
```

## Core Parameters

### Input/Output

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--input` | path | *required* | Input file (H5AD, RDS, or MTX directory) |
| `--input_type` | string | `h5ad` | Format: `h5ad`, `rds`, or `mtx` |
| `--outdir` | path | `./results` | Output directory for all results |

### Quality Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--use_mad_thresholds` | boolean | `true` | Use MAD-based adaptive thresholds |
| `--mad_multiplier` | float | `5.0` | MAD threshold stringency (higher = more permissive) |
| `--min_genes` | integer | `200` | Minimum genes per cell (if MAD disabled) |
| `--min_counts` | integer | `500` | Minimum UMI counts per cell (if MAD disabled) |
| `--max_mito_percent` | float | `20.0` | Maximum mitochondrial content (%) |

### Batch Integration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_integration` | boolean | `false` | Enable Harmony batch correction |
| `--batch_key` | string | `batch` | Column name for batch variable in `.obs` |
| `--n_top_genes` | integer | `2000` | Highly variable genes for integration |
| `--n_pcs` | integer | `50` | Principal components to compute |
| `--harmony_theta` | float | `2.0` | Harmony diversity clustering penalty |
| `--harmony_max_iter` | integer | `10` | Maximum Harmony iterations |

### Annotation

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--annotation.enabled` | boolean | `false` | Enable automated cell type annotation |
| `--annotation.tools` | string | `celltypist` | Annotation tools: `celltypist`, `azimuth`, or `celltypist,azimuth` |
| `--annotation.celltypist_model` | string | `Immune_All_Low.pkl` | CellTypist model name |

## Pipeline Output

```
results/
├── qc/
│   ├── qc_filtered.h5ad              # Filtered cells passing QC
│   ├── qc_report.json                # QC metrics and thresholds
│   └── plots/                        # QC visualization PDFs
├── normalized_integrated/
│   ├── normalized_integrated.h5ad    # Batch-corrected data with UMAP
│   └── integration_report.json       # Integration quality metrics
├── annotation/                       # (if enabled)
│   ├── celltypist/
│   │   └── sample_celltypist.csv     # CellTypist predictions
│   ├── azimuth/                      # (if enabled)
│   │   └── sample_azimuth.csv        # Azimuth predictions
│   └── sample_annotated.h5ad         # Final annotated dataset
└── pipeline_info/
    ├── execution_report.html         # Nextflow execution report
    └── pipeline_dag.html             # Workflow DAG visualization
```

## Interactive Dashboard

The scAnnex dashboard provides an interactive web-based interface for exploring single-cell RNA-seq analysis results. Built with R Shiny, it offers real-time data visualization with support for both small datasets and large-scale analyses (>100k cells).

### Features

- **Data Input Tab**: Load H5AD files with automatic format detection and backed mode for large datasets
- **QC Overview Tab**: Interactive quality control metrics with before/after filtering comparison
- **Clustering & UMAP Tab**: 
  - Interactive UMAP visualization with WebGL acceleration
  - Color by cell types, QC metrics, or batch information
  - Adjustable point size and opacity
  - Cell metadata table with search and filtering
- **Gene Expression Tab**: 
  - Real-time gene expression visualization on UMAP
  - Continuous color scales (Viridis)
  - Expression value hover tooltips
- **About Tab**: Pipeline version and configuration information

### Quick Start

#### Option 1: Conda Environment (Recommended - No sudo required)

```bash
cd dashboard
./setup_dashboard.sh      # Auto-detects and sets up environment
./launch_dashboard.sh     # Launches on http://localhost:8888
```

#### Option 2: Docker

```bash
cd dashboard
docker build -t scannex-dashboard .
docker run -p 3838:3838 -v $(pwd)/../results:/data scannex-dashboard
# Access at: http://localhost:3838
```

#### Option 3: HPC with Apptainer/Singularity

```bash
cd dashboard
apptainer build scannex-dashboard.sif scannex-dashboard.def
apptainer run --bind ./results:/data scannex-dashboard.sif
```

### Dashboard Usage

1. **Load your data**: Enter the path to your `*_annotated.h5ad` file from pipeline results
2. **Explore QC metrics**: View filtering statistics and quality control plots
3. **Visualize clusters**: Interactive UMAP colored by cell types or metadata
4. **Inspect gene expression**: Search for marker genes and view their expression patterns

### Example Dataset

The dashboard has been tested with the included PBMC dataset:
- **935 cells** × **14,521 genes**
- **6 cell types** annotated:
  - Tcm/Naive helper T cells (36.6%)
  - Classical monocytes (23.7%)
  - Naive B cells (18.8%)
  - MAIT cells (13.0%)
  - CD16+ NK cells (5.6%)
  - Non-classical monocytes (2.2%)

### Memory Optimization

The dashboard automatically detects file size and uses backed mode (memory-efficient reading) for large datasets:
- **< 500 MB**: Full in-memory loading (fast, complete feature access)
- **≥ 500 MB**: Backed mode (low memory, suitable for HPC)

### Troubleshooting

For detailed setup instructions and troubleshooting:
- **Simple guide**: `dashboard/README_SIMPLE.md`
- **Quick reference**: `dashboard/QUICKSTART.md`
- **Manual setup**: `dashboard/MANUAL_LAUNCH.md`
- **WSL2 issues**: `dashboard/TROUBLESHOOTING_WSL2.md`

### Dashboard Requirements

**Conda environment** (recommended):
- R 4.3+
- Python 3.10+ with anndata, scanpy
- R packages: shiny, shinydashboard, plotly, DT, reticulate

All dependencies are managed automatically by `setup_dashboard.sh`.

## Execution Profiles

scAnnex supports multiple execution profiles:

| Profile | Description | Container Engine |
|---------|-------------|------------------|
| `docker` | Local execution with Docker | Docker |
| `singularity` | HPC execution with Singularity | Singularity/Apptainer |
| `test` | Minimal test dataset | None (combine with `docker` or `singularity`) |

### HPC Execution (SLURM)

```bash
nextflow run main.nf \
  --input sample_data.h5ad \
  --outdir results \
  -profile singularity \
  -process.executor slurm \
  -process.queue normal
```

## Architecture

scAnnex implements a modular Nextflow workflow with the following components:

1. **UNIFY_INPUT**: Converts RDS/MTX to standardized H5AD format
2. **QUALITY_CONTROL**: MAD-based adaptive QC with visualization
3. **DOUBLET_DETECTION**: Scrublet-based doublet scoring (optional)
4. **NORMALIZE_INTEGRATE**: Log-normalization and Harmony integration
5. **H5AD_TO_RDS**: H5AD → Seurat RDS bridge (for R-based annotation)
6. **AUTO_ANNOT_CELLTYPIST**: Python-based CellTypist annotation
7. **AUTO_ANNOT_AZIMUTH**: R-based Azimuth annotation (optional)
8. **AUTO_ANNOT_MERGE**: Consolidates multi-tool predictions

## Data Conventions

### AnnData Structure

All outputs follow standardized AnnData conventions:

```python
adata.obs            # Cell metadata (sample_id, batch, QC metrics, clusters)
adata.var            # Gene metadata (highly_variable, mean, dispersion)
adata.X              # Normalized expression matrix (log1p)
adata.layers['counts']  # Raw counts (preserved)
adata.obsm['X_pca']     # PCA coordinates (50 components)
adata.obsm['X_pca_harmony']  # Harmony-corrected PCA (if integration enabled)
adata.obsm['X_umap']    # UMAP coordinates (2D)
adata.uns['qc']         # QC parameters and thresholds
adata.uns['integration']  # Integration parameters
```

### Annotation Columns

Automated annotations are stored as:

```python
adata.obs['celltype_celltypist']  # CellTypist predictions
adata.obs['celltype_azimuth']     # Azimuth predictions (if enabled)
adata.obs['celltype_consensus']   # Merged consensus labels (if multi-tool)
```

## License

MIT License. See `LICENSE` file for details.

## Software Dependencies

- **Python**: scanpy (1.9+), anndata, numpy, pandas, harmony-pytorch
- **R**: Seurat, SeuratDisk, reticulate (for RDS conversion and Azimuth)
- **Nextflow**: 23.04.0+
- **Dashboard**: R Shiny, reticulate, plotly, DT

All dependencies are managed via Docker/Singularity containers.
