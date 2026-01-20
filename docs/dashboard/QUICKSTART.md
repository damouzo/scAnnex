# scAnnex Dashboard Quick Start Guide

This guide provides simple instructions for launching the interactive dashboard to explore your scRNA-seq results.

## Choose Your Environment

### Option 1: HPC with Apptainer/Singularity (RECOMMENDED for HPC)

**Best for:** Computing clusters, HPC systems, shared environments

1. **Build the container** (one-time setup):
   ```bash
   cd dashboard
   ./build_dashboard_container.sh apptainer
   ```

2. **Launch on SLURM cluster**:
   ```bash
   sbatch launch_dashboard.slurm
   ```

3. **Set up SSH tunnel from your local machine**:
   ```bash
   # Check the job log for the node and port
   cat dashboard_<jobid>.log
   
   # Create SSH tunnel (replace values from log)
   ssh -N -L 8888:<compute-node>:8888 username@hpc-login.edu
   ```

4. **Open in browser**:
   ```
   http://localhost:8888
   ```

### Option 2: Local with Docker

**Best for:** Personal computers, development environments

1. **Build the container** (one-time setup):
   ```bash
   cd dashboard
   ./build_dashboard_container.sh docker
   ```

2. **Run the dashboard**:
   ```bash
   docker run -d -p 3838:3838 \
     --name scannex-dashboard \
     -v $(pwd)/../results_slc_first_run:/data \
     scannex-dashboard:latest
   ```

3. **Open in browser**:
   ```
   http://localhost:3838
   ```

4. **Stop the dashboard**:
   ```bash
   docker stop scannex-dashboard
   docker rm scannex-dashboard
   ```

### Option 3: Manual Python Exploration (No Container Needed)

**Best for:** Quick data inspection without web interface

```python
import anndata as ad
import scanpy as sc

# Load your annotated data
adata = ad.read_h5ad('results_slc_first_run/auto/PBMC_TEST_annotated.h5ad')

# View cell type counts
print(adata.obs['predicted_labels'].value_counts())

# Plot UMAP colored by cell type
sc.pl.umap(adata, color='predicted_labels', legend_loc='on data')

# Plot gene expression
sc.pl.umap(adata, color=['CD3D', 'CD14', 'CD79A'], cmap='viridis')

# Export for other tools
adata.write_h5ad('annotated_for_export.h5ad')
```

## Dashboard Features

Once the dashboard is running, you can:

### Data Input Tab
- Load your H5AD file
- Specify QC directory for plots
- Enable backed mode for large datasets (>100k cells)

### QC Overview Tab
- Cell counts before/after filtering
- QC metric distributions
- Threshold visualizations

### Clustering & UMAP Tab
- Interactive UMAP plot
- Color by any metadata (cell types, batch, condition, etc.)
- Searchable cell metadata table
- Adjustable point size and opacity

### Gene Expression Tab
- Search for any gene
- Visualize expression on UMAP
- Interactive viridis color scale

## File Locations

Your successful run output:
```
results_slc_first_run/
├── auto/
│   ├── PBMC_TEST_annotated.h5ad     ← Use this for dashboard
│   └── PBMC_TEST_celltypist.csv     ← Cell type predictions
├── qc/                               ← QC plots and reports
├── doublet_detection/
├── standard/
└── unified_input/
```

## Troubleshooting

### Container won't build
- **Apptainer**: May need `--fakeroot` flag or sudo access
- **Docker**: Check Docker daemon is running: `docker ps`

### Can't access dashboard
- **HPC**: Verify SSH tunnel is running and port matches
- **Docker**: Check container is running: `docker ps`
- **Firewall**: Ensure port is open

### Dashboard loads but no data
- Check file path in "Data Input" tab
- Ensure H5AD file has required fields (UMAP coords, metadata)
- Check logs for Python/R errors

### Slow performance
- Enable "backed mode" for large datasets
- Reduce UMAP point size
- Check available memory

## Getting Help

- Dashboard logs: `dashboard_<jobid>.log` (SLURM) or `docker logs scannex-dashboard` (Docker)
- Python errors: Check if scanpy/anndata can read the file directly
- Container issues: Verify image built successfully

## Next Steps

After exploring your data:
- Export high-resolution plots for publications
- Subset interesting cell populations
- Integrate with additional tools (Seurat, CellxGene, etc.)
- Share the container definition for reproducibility

---

For more details, see `dashboard/README.md`
