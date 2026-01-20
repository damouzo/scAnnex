# scAnnex Dashboard Documentation

Complete documentation for the scAnnex interactive visualization dashboard.

## Quick Links

### Getting Started
- **[Quick Start Guide](QUICKSTART.md)** - Fast setup and launch (5 minutes)
- **[Simple User Guide](README_SIMPLE.md)** - Step-by-step beginner-friendly instructions
- **[Manual Launch Instructions](MANUAL_LAUNCH.md)** - Detailed setup for all deployment methods

### Troubleshooting
- **[WSL2 Troubleshooting](TROUBLESHOOTING_WSL2.md)** - Windows Subsystem for Linux issues
- **[Firewall Configuration](FIREWALL_FIX.md)** - Network access and firewall solutions

### Advanced Topics
- **[GitHub Actions Costs](GITHUB_ACTIONS_COSTS.md)** - CI/CD configuration and cost analysis

---

## Dashboard Overview

The scAnnex dashboard is an interactive R Shiny application for exploring single-cell RNA-seq analysis results. It provides:

- **Interactive UMAP visualization** with WebGL acceleration
- **Quality control metrics** and filtering statistics
- **Gene expression visualization** with real-time updates
- **Cell type annotations** with confidence scores
- **Metadata exploration** with search and filtering

### Key Features

✅ **Zero-configuration setup** - Auto-detects environment (Conda/Docker/Apptainer)  
✅ **Memory efficient** - Supports large datasets (>100k cells) with backed mode  
✅ **HPC compatible** - Works on SLURM clusters with Apptainer  
✅ **No sudo required** - Conda-based deployment works everywhere  
✅ **Reproducible** - Containerized environments with pinned dependencies

---

## Installation Methods

### Method 1: Conda (Recommended)

**Best for:** Local machines, development, no sudo access

```bash
cd dashboard
./setup_dashboard.sh      # Auto-detects and sets up environment
./launch_dashboard.sh     # Launches on http://localhost:8888
```

**Pros:**
- No container engine required
- Works everywhere (Windows/macOS/Linux)
- Fast startup
- Easy to debug

**Cons:**
- Requires Conda/Mamba installation
- Uses local storage for environment

---

### Method 2: Docker

**Best for:** Production deployments, web servers

```bash
cd dashboard
docker build -t scannex-dashboard .
docker run -p 3838:3838 -v $(pwd)/../results:/data scannex-dashboard
```

Access at: http://localhost:3838

**Pros:**
- Isolated environment
- Easy to deploy multiple instances
- Port mapping flexibility

**Cons:**
- Requires Docker installation (may need sudo)
- Larger image size
- Slower startup

---

### Method 3: Apptainer/Singularity (HPC)

**Best for:** HPC clusters, shared environments

```bash
cd dashboard
apptainer build scannex-dashboard.sif scannex-dashboard.def
apptainer run --bind ./results:/data scannex-dashboard.sif
```

Or submit SLURM job:
```bash
sbatch launch_dashboard.slurm
```

**Pros:**
- No root/sudo required
- Works on HPC clusters
- Single-file container (.sif)
- SLURM integration

**Cons:**
- Requires Apptainer/Singularity installation
- Longer build time
- Less common than Docker

---

## Usage Guide

### 1. Launch Dashboard

Choose your method:
```bash
# Conda
./launch_dashboard.sh

# Docker
./run_dashboard.sh run

# Apptainer
apptainer run scannex-dashboard.sif
```

### 2. Load Data

1. Navigate to **Data Input** tab
2. Enter path to your `*_annotated.h5ad` file:
   ```
   /path/to/results/auto/SAMPLE_annotated.h5ad
   ```
3. (Optional) Enter QC results directory:
   ```
   /path/to/results/qc/
   ```
4. Adjust **backed mode** setting:
   - ✅ **OFF** for datasets < 50k cells (faster)
   - ✅ **ON** for datasets > 50k cells (memory efficient)
5. Click **Load Data**

### 3. Explore QC Metrics

Navigate to **QC Overview** tab to view:
- Cell counts before/after filtering
- QC metric distributions
- Filtering thresholds applied
- Before/after comparison plots

### 4. Visualize Clusters

Navigate to **Clustering & UMAP** tab:

1. **Select coloring variable:**
   - `predicted_labels` - Cell types from CellTypist
   - `leiden_X` - Clustering resolution (X = 0.1-0.9)
   - `batch` - Batch information
   - QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`)

2. **Adjust visualization:**
   - Point size (1-10)
   - Opacity (0.1-1.0)

3. **Interact with plot:**
   - Hover for cell information
   - Zoom with scroll wheel
   - Pan by clicking and dragging
   - Box select for detailed view

4. **Explore metadata table:**
   - Search for specific cells
   - Filter by any column
   - Export filtered data

### 5. Visualize Gene Expression

Navigate to **Gene Expression** tab:

1. Enter a gene name (e.g., `CD3D`, `CD14`, `MS4A1`)
2. Click **Plot Expression**
3. View expression levels on UMAP with Viridis color scale
4. Hover to see exact expression values

**Marker genes to try:**
- **T cells**: CD3D, CD3E, CD4, CD8A
- **B cells**: CD79A, MS4A1, CD19
- **Monocytes**: CD14, FCGR3A (CD16)
- **NK cells**: NKG7, GNLY, NCAM1

---

## Dashboard Architecture

```
dashboard/
├── app.R                    # Main Shiny app entry point
├── global.R                 # Data loading functions, Python setup
├── server.R                 # Server logic (reactive data)
├── ui.R                     # User interface layout
├── environment_dashboard.yml # Conda environment specification
├── Dockerfile               # Docker container definition
├── scannex-dashboard.def    # Apptainer container definition
├── setup_dashboard.sh       # Auto-setup script
├── launch_dashboard.sh      # Universal launcher
└── launch_dashboard.slurm   # SLURM job template
```

### Data Flow

```
H5AD File → load_h5ad_data() → {
  adata (AnnData object)
  metadata (R data.frame)
  umap_coords (R data.frame)
  var_info (gene info)
} → Reactive Values → Plotly Visualizations
```

### Technologies Used

**Backend:**
- R Shiny (web framework)
- reticulate (R-Python bridge)
- anndata (Python, data loading)
- scanpy (Python, single-cell analysis)

**Frontend:**
- shinydashboard (UI framework)
- plotly (interactive plots with WebGL)
- DT (interactive tables)

---

## Performance Optimization

### Memory Usage

The dashboard automatically optimizes memory based on dataset size:

| Dataset Size | Mode | Memory Usage | Load Time |
|--------------|------|--------------|-----------|
| < 10k cells | In-memory | ~100-500 MB | 1-5 sec |
| 10-50k cells | In-memory | ~500 MB - 2 GB | 5-15 sec |
| 50-100k cells | Backed | ~100-200 MB | 10-30 sec |
| > 100k cells | Backed | ~200-500 MB | 30-60 sec |

### Speed Tips

1. **Use backed mode** for large datasets (>50k cells)
2. **Reduce point size** in UMAP (improves rendering)
3. **Lower opacity** to see overlapping points
4. **Use WebGL** rendering (enabled by default with `scattergl`)

---

## Troubleshooting

### Common Issues

#### Dashboard won't start

**Check:**
1. Conda environment activated: `conda activate scannex-dashboard`
2. Python path set: `echo $RETICULATE_PYTHON`
3. Required files present: `ls app.R global.R server.R ui.R`

**Solution:**
```bash
cd dashboard
./setup_dashboard.sh   # Recreate environment
./launch_dashboard.sh  # Launch again
```

#### Data won't load

**Check:**
1. File path is correct and absolute (not relative)
2. File exists: `ls -lh /path/to/file.h5ad`
3. File has UMAP coordinates in `.obsm['X_umap']`

**Solution:**
```bash
# Verify h5ad file structure
python3 -c "import anndata; adata = anndata.read_h5ad('/path/to/file.h5ad'); print(adata); print(adata.obsm.keys())"
```

#### UMAP plot doesn't display

**Check:**
1. UMAP coordinates exist in data
2. No NaN/Inf values in coordinates
3. Backed mode setting (try disabling for small files)

**Solution:**
- Disable backed mode for files < 500 MB
- Check browser console for JavaScript errors (F12)

#### Gene expression plot is blank

**Check:**
1. Gene name is correct (case-sensitive)
2. Gene exists in dataset
3. Gene has non-zero expression

**Solution:**
```bash
# List available genes
python3 -c "import anndata; adata = anndata.read_h5ad('/path/to/file.h5ad'); print(adata.var_names[:20])"
```

For more detailed troubleshooting:
- **[WSL2 Issues](TROUBLESHOOTING_WSL2.md)**
- **[Firewall Issues](FIREWALL_FIX.md)**
- **[Manual Setup](MANUAL_LAUNCH.md)**

---

## Development

### Running in Development Mode

```bash
cd dashboard
conda activate scannex-dashboard
R -e "shiny::runApp('.', host='127.0.0.1', port=8888)"
```

### Testing Changes

```bash
# Test data loading
Rscript test_dashboard_full.R

# Test specific functions
R -e "source('global.R'); load_h5ad_data('/path/to/file.h5ad')"
```

### Adding New Features

1. **New plot type:** Add function to `global.R`
2. **New UI element:** Add to `ui.R`
3. **New reactive behavior:** Add to `server.R`
4. **Test:** Restart app and verify functionality

---

## Citations

If you use the scAnnex dashboard in your research, please cite:

```bibtex
@software{scannex_dashboard,
  author = {scAnnex Development Team},
  title = {scAnnex Dashboard: Interactive Single-Cell RNA-seq Visualization},
  year = {2026},
  url = {https://github.com/[username]/scAnnex}
}
```

---

## Support

- **GitHub Issues:** https://github.com/[username]/scAnnex/issues
- **Documentation:** https://github.com/[username]/scAnnex/tree/main/docs/dashboard
- **Email:** [Your contact]

---

**Last Updated:** January 20, 2026  
**Dashboard Version:** 1.0  
**Status:** Production Ready
