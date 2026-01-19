# scAnnex Dashboard

Interactive R Shiny dashboard for visualizing scRNA-seq analysis results from the scAnnex pipeline.

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

## Docker Deployment (Recommended)

### Build the Docker image:
```bash
cd dashboard
docker build -t scannex-dashboard:latest .
```

### Run the dashboard:
```bash
# Basic run (uses default test data path)
docker run -d -p 3838:3838 \
  --name scannex-dashboard \
  -v $(pwd)/../test_data/analytical_core_results:/srv/shiny-server/data \
  scannex-dashboard:latest

# Access at: http://localhost:3838
```

### With custom data directory:
```bash
docker run -d -p 3838:3838 \
  --name scannex-dashboard \
  -v /path/to/your/results:/srv/shiny-server/data \
  scannex-dashboard:latest
```

### Stop the dashboard:
```bash
docker stop scannex-dashboard
docker rm scannex-dashboard
```

## Local Development (if R + packages installed)

### Install required R packages:
```R
install.packages(c(
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "reticulate",
  "plotly",
  "DT",
  "ggplot2",
  "viridis",
  "data.table",
  "jsonlite"
))
```

### Install Python packages:
```bash
pip3 install scanpy anndata h5py numpy pandas
```

### Run locally:
```bash
R -e "shiny::runApp('dashboard', port=3838, host='0.0.0.0')"
```

## File Structure

```
dashboard/
├── Dockerfile              # Docker container definition
├── app.R                   # Main Shiny app entry point
├── global.R                # Global functions and data loading
├── ui.R                    # User interface layout
├── server.R                # Server logic and reactive functions
├── modules/                # (Future) Modular Shiny components
└── www/                    # Static assets (CSS, images)
```

## Data Requirements

The dashboard expects the following data structure:

```
data/
├── qc_filtered.h5ad        # QC-filtered H5AD file
└── qc_results/             # QC results directory
    ├── qc_report.json      # QC metrics and thresholds
    ├── qc_before_*.png     # Pre-filtering plots
    └── qc_after_*.png      # Post-filtering plots
```

## Performance Considerations

### Backed Mode (>50k cells)
The dashboard uses AnnData's backed mode for large datasets:
- Metadata loaded into memory
- Expression matrices stay on disk
- Gene expression loaded on-demand

### WebGL Acceleration
UMAP plots use `scattergl` for hardware-accelerated rendering of large point clouds.

### Recommended Limits
- **Optimal**: <100k cells, interactive and fast
- **Acceptable**: 100k-500k cells, use backed mode
- **Large**: >500k cells, consider subsetting for visualization

## Customization

### Change default data path:
Edit the `DEFAULT_DATA_PATH` variable in `global.R` or set environment variable:
```bash
docker run -d -p 3838:3838 \
  -e SCANNEX_DATA_PATH=/custom/path \
  -v /host/path:/custom/path \
  scannex-dashboard:latest
```

### Add custom CSS:
Place CSS files in `www/` directory and reference in `ui.R`.

### Add new tabs:
1. Add `menuItem` in `ui.R` sidebar
2. Add `tabItem` in `ui.R` body
3. Add reactive logic in `server.R`

## Troubleshooting

### Dashboard won't start
- Check Docker logs: `docker logs scannex-dashboard`
- Verify data directory is mounted correctly
- Ensure H5AD file exists and is readable

### Can't read H5AD file
- Verify Python + scanpy are installed: `python3 -c "import scanpy; print(scanpy.__version__)"`
- Check file permissions
- Try loading in Python directly to test

### Plots not showing
- Check QC results directory contains PNG files
- Verify file paths in Data Input tab
- Check browser console for JavaScript errors

### Slow performance
- Enable backed mode for large datasets
- Reduce UMAP point size and opacity
- Consider subsetting data for visualization

## Production Deployment

For production use with multiple users, consider:

- **ShinyProxy**: Containerized per-user sessions
- **Posit Connect**: Commercial Shiny server with authentication
- **Nginx**: Reverse proxy for load balancing
- **Resource limits**: Set Docker memory/CPU limits

Example with resource limits:
```bash
docker run -d -p 3838:3838 \
  --name scannex-dashboard \
  --memory="4g" \
  --cpus="2" \
  -v /data:/srv/shiny-server/data \
  scannex-dashboard:latest
```

## License

See main scAnnex project LICENSE file.

## Contact

For issues and feature requests, see the main scAnnex repository.
