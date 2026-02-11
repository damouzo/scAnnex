# scAnnex Dashboard

Interactive R Shiny dashboard for visualizing and exploring scRNA-seq analysis results from the scAnnex pipeline.

---

## Quick Start

### Setup (one-time)

```bash
cd dashboard
./setup_dashboard.sh
```

This will automatically detect and set up the appropriate environment (Conda, Docker, or Singularity).

### Launch

```bash
cd dashboard
./launch_dashboard.sh [path/to/results]
```

The dashboard will be available at: http://localhost:3838

Press Ctrl+C to stop the dashboard.

---

## Requirements

The dashboard works with any of these options:

- **Conda/Mamba** (recommended, no sudo needed)
- **Docker** (for containers)
- **Apptainer/Singularity** (for HPC environments)

---

## Features

### QC Overview
- Summary metrics (cells before/after filtering)
- Interactive QC plots (violin, scatter, distributions)
- Cell attrition tracking

### Clustering
- Interactive UMAP visualization
- Color by metadata (batch, sample, condition)
- Multiple clustering resolutions

### Gene Expression
- Search genes and visualize expression on UMAP
- Interactive hover tooltips
- Expression distribution plots

### Annotation Station
- Rule-based cell type annotation
- Combine clustering and auto-annotations
- Export annotations back to H5AD

---

## File Structure

```
dashboard/
├── app.R                        # Main entry point
├── global.R                     # Functions and data loading
├── ui.R                         # User interface
├── server.R                     # Server logic
├── environment_dashboard.yml    # Conda environment
├── scannex-dashboard.def        # Singularity definition
├── setup_dashboard.sh           # Setup script
└── launch_dashboard.sh          # Launch script
```

---

## Data Requirements

The dashboard expects H5AD files from scAnnex with:

- `adata.obs['predicted_labels']` - Cell type annotations
- `adata.obsm['X_umap']` - UMAP coordinates
- `adata.obsm['X_pca']` - PCA coordinates
- QC metrics and clustering results

