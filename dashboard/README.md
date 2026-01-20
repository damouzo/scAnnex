# scAnnex Dashboard

Interactive R Shiny dashboard for visualizing scRNA-seq analysis results from the scAnnex pipeline.

## 📚 Complete Documentation

**For detailed setup instructions, troubleshooting, and usage guides, see:**

👉 **[Complete Dashboard Documentation](../docs/dashboard/README.md)**

### Quick Links
- [Quick Start Guide](../docs/dashboard/QUICKSTART.md) - 5-minute setup
- [Simple User Guide](../docs/dashboard/README_SIMPLE.md) - Step-by-step for beginners
- [WSL2 Troubleshooting](../docs/dashboard/TROUBLESHOOTING_WSL2.md)
- [Firewall Configuration](../docs/dashboard/FIREWALL_FIX.md)

---

## Quick Start

### Option 1: Conda (Recommended - No sudo required)
```bash
cd dashboard
./setup_dashboard.sh      # Auto-detects and sets up environment
./launch_dashboard.sh     # Launches on http://localhost:8888
```

### Option 2: Docker
```bash
cd dashboard
docker build -t scannex-dashboard .
docker run -p 3838:3838 -v $(pwd)/../results:/data scannex-dashboard
# Access at: http://localhost:3838
```

### Option 3: Apptainer/Singularity (HPC)
```bash
cd dashboard
apptainer build scannex-dashboard.sif scannex-dashboard.def
apptainer run --bind ./results:/data scannex-dashboard.sif
```

---

## Features

### Tab 1: Data Input
- Load H5AD files with backed mode support for large datasets (>100k cells)
- Configure QC results directory
- Real-time data loading status

### Tab 2: QC Overview
- Summary metrics (cells before/after, retention rates)
- Interactive QC metrics tables
- Before/after filtering plots (violin, scatter, distributions)
- MAD-based threshold visualization

### Tab 3: Clustering & UMAP
- Interactive UMAP plots with plotly (WebGL-accelerated)
- Customizable point size and opacity
- Color by metadata (batch, sample_id, condition, etc.)
- Searchable cell metadata table

### Tab 4: Gene Expression
- Gene search with expression visualization on UMAP
- Viridis color scale for expression levels
- Interactive hover tooltips

### Tab 5: About
- Project information and documentation

## Technology Stack

- **R**: 4.3.3
- **Shiny**: Interactive web application framework
- **reticulate**: Python integration for reading H5AD files
- **plotly**: Interactive plots with WebGL
- **DT**: Interactive data tables
- **Python**: scanpy/anndata for H5AD reading

---

## File Structure

```
dashboard/
├── app.R                        # Main Shiny app entry point
├── global.R                     # Global functions and data loading
├── ui.R                         # User interface layout
├── server.R                     # Server logic and reactive functions
├── environment_dashboard.yml    # Conda environment specification
├── Dockerfile                   # Docker container definition
├── scannex-dashboard.def        # Apptainer container definition
├── setup_dashboard.sh           # Auto-setup script
├── launch_dashboard.sh          # Universal launcher
└── launch_dashboard.slurm       # SLURM job template
```

---

## Data Requirements

The dashboard expects H5AD files from the scAnnex pipeline with:

```python
adata.obs['predicted_labels']    # Cell type annotations
adata.obs['celltypist_score']    # Confidence scores
adata.obsm['X_umap']              # UMAP coordinates
adata.obsm['X_pca']               # PCA coordinates
# QC metrics, clustering results, etc.
```

Example data path:
```
results/
├── auto/
│   └── SAMPLE_annotated.h5ad     # Main annotated dataset
└── qc/
    ├── qc_report.json            # QC metrics
    └── plots/                     # QC visualization plots
```

---

## Performance

### Memory Usage by Dataset Size

| Dataset Size | Mode | Memory | Load Time |
|--------------|------|--------|-----------|
| < 10k cells | In-memory | ~100-500 MB | 1-5 sec |
| 10-50k cells | In-memory | ~500 MB - 2 GB | 5-15 sec |
| 50-100k cells | Backed | ~100-200 MB | 10-30 sec |
| > 100k cells | Backed | ~200-500 MB | 30-60 sec |

The dashboard automatically optimizes memory based on file size (<500 MB = in-memory mode).

---

## Deployment Options Comparison

| Method | Setup Time | Requires | Best For |
|--------|------------|----------|----------|
| **Conda** | 5 min | Conda/Mamba | Local development, no sudo |
| **Docker** | 10 min | Docker | Production, web servers |
| **Apptainer** | 15 min | Apptainer | HPC clusters |

---

## Troubleshooting

### Quick Fixes

**Dashboard won't start:**
```bash
cd dashboard
./setup_dashboard.sh   # Recreate environment
./launch_dashboard.sh  # Try again
```

**Data won't load:**
```bash
# Verify file structure
python3 -c "import anndata; adata = anndata.read_h5ad('/path/to/file.h5ad'); print(adata.obsm.keys())"
```

**UMAP plot blank:**
- Disable backed mode for small files (<500 MB)
- Check browser console (F12) for errors

**For detailed troubleshooting:** See [Complete Documentation](../docs/dashboard/README.md)

---

## Development

### Testing
```bash
# Test all functionality
Rscript test_dashboard_full.R

# Launch in development mode
R -e "shiny::runApp('.', host='127.0.0.1', port=8888)"
```

### Adding Features
1. **New plot:** Add function to `global.R`
2. **New UI element:** Edit `ui.R`
3. **New logic:** Edit `server.R`
4. **Test:** Restart app

---

## Production Deployment

For production with multiple users:
- **ShinyProxy**: Per-user containerized sessions
- **Posit Connect**: Commercial with authentication  
- **Nginx**: Reverse proxy + load balancing

Example with resource limits:
```bash
docker run -d -p 3838:3838 \
  --name scannex-dashboard \
  --memory="4g" --cpus="2" \
  -v /data:/srv/shiny-server/data \
  scannex-dashboard:latest
```

---

## Support

- **Issues:** https://github.com/[username]/scAnnex/issues
- **Documentation:** [docs/dashboard/](../docs/dashboard/)
- **Email:** [Your contact]

---

**Version:** 1.0  
**Last Updated:** January 20, 2026  
**Status:** Production Ready
