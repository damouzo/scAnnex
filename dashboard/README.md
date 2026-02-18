# scAnnex Dashboard

Interactive R Shiny dashboard for visualizing and exploring scRNA-seq analysis results from the scAnnex pipeline.

---

## Quick Start

### Local Execution

```bash
cd dashboard
./launch_dashboard.sh [path/to/results]
```

The dashboard will be available at: http://localhost:3838

Press Ctrl+C to stop the dashboard.

### HPC Execution (SLURM)

**Recommended method** for HPC clusters:

```bash
cd dashboard
./launch_dashboard_hpc.sh [path/to/results]
```

This will:
1. Request a SLURM interactive job on a compute node
2. Display SSH tunnel instructions for your local machine
3. Launch the dashboard server

See [HPC Dashboard Guide](../docs/internal/HPC_DASHBOARD_GUIDE.md) for detailed instructions.

---

## Deployment Modes

### Mode 1: Local (Laptop/Workstation)

Best for local analysis on your own computer.

**Launch:**
```bash
./launch_dashboard.sh /path/to/results
```

**Access:** http://localhost:3838

**Requirements:** 
- Conda/Mamba installed
- Sufficient RAM for your dataset

### Mode 2: HPC Production (SLURM Interactive)

Best for HPC clusters with SLURM scheduler. Runs on dedicated compute node with configurable resources.

**Launch:**
```bash
# Basic (4 CPUs, 8GB RAM, 4 hours)
./launch_dashboard_hpc.sh /path/to/results

# Custom resources
./launch_dashboard_hpc.sh --cpus 8 --mem 16G --time 8:00:00 /path/to/results
```

**Access:** http://localhost:3838 (via SSH tunnel)

**Features:**
- Dedicated CPU and memory resources
- Auto-port detection
- Configurable time limits
- SSH tunnel instructions provided

### Mode 3: HPC Quick Testing

For quick testing on HPC without requesting a compute node. Not recommended for production use.

**Launch:**
```bash
./launch_dashboard.sh /path/to/results
# Confirm when prompted about HPC environment
```

---

## HPC Resource Recommendations

Choose resources based on your dataset size:

| Dataset Size | CPUs | Memory | Time Limit | Command |
|--------------|------|--------|------------|---------|
| < 10k cells | 2-4 | 4-8G | 2-4h | `./launch_dashboard_hpc.sh` |
| 10k-50k cells | 4-8 | 8-16G | 4-8h | `./launch_dashboard_hpc.sh --cpus 4 --mem 8G` |
| 50k-100k cells | 8-16 | 16-32G | 8-12h | `./launch_dashboard_hpc.sh --cpus 8 --mem 16G` |
| > 100k cells | 16+ | 32G+ | 12-24h | `./launch_dashboard_hpc.sh --cpus 16 --mem 32G --time 12:00:00` |

**Notes:**
- Dashboard is moderately CPU-intensive (1-4 cores for most operations)
- Memory requirements scale with dataset size
- Time limits prevent runaway jobs

---

## Requirements

The dashboard works with any of these options:

- **Conda/Mamba** (recommended, no sudo needed)
- **Docker** (for containers)
- **Apptainer/Singularity** (for HPC environments)

Environment setup is automatic on first launch.

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
├── app.R                          # Main entry point
├── global.R                       # Functions and data loading
├── ui.R                           # User interface
├── server.R                       # Server logic
├── environment_dashboard.yml      # Conda environment
├── launch_dashboard.sh            # Local/HPC quick launcher
├── launch_dashboard_hpc.sh        # HPC production launcher (SLURM)
├── test_ssh_tunnel.sh             # SSH tunnel diagnostic tool
└── README.md                      # This file
```

---

## Troubleshooting

### Dashboard Won't Load

1. **Check SSH tunnel** (HPC only):
   ```bash
   ./test_ssh_tunnel.sh
   ```

2. **Verify port availability**:
   - Default port: 3838
   - Script will auto-increment if busy (3839, 3840, etc.)

3. **Check dashboard logs**:
   ```bash
   cat dashboard.log
   ```

### Common Issues

**Issue:** "Port 3838 already in use"  
**Solution:** Script will automatically use next available port. Update browser URL accordingly.

**Issue:** "Cannot connect to http://localhost:3838"  
**Solution (HPC):** Ensure SSH tunnel is running. See SSH tunnel command in dashboard launcher output.

**Issue:** "No annotated .h5ad file found"  
**Solution:** Verify pipeline completed and results directory is correct.

---

## Data Requirements

The dashboard expects H5AD files from scAnnex with:

- `adata.obs['predicted_labels']` - Cell type annotations
- `adata.obsm['X_umap']` - UMAP coordinates
- `adata.obsm['X_pca']` - PCA coordinates
- QC metrics and clustering results

---

## Additional Documentation

- **HPC Deployment:** [HPC Dashboard Guide](../docs/internal/HPC_DASHBOARD_GUIDE.md)
- **Pipeline Documentation:** [Main README](../README.md)
- **Apocrita HPC:** https://docs.hpc.qmul.ac.uk

---

## Support

For issues or questions:

- **HPC clusters (Apocrita):** its-research-support@qmul.ac.uk
- **scAnnex pipeline:** Check GitHub issues
- **Dashboard problems:** Run `./test_ssh_tunnel.sh` for diagnostics

