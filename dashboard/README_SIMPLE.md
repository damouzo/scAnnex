# 🚀 scAnnex Dashboard - Zero-Setup Launch

**Interactive visualization for your scRNA-seq analysis results**

No sudo required • HPC-compatible • One command to launch

---

## Quickest Start (2 steps)

```bash
# 1. One-time setup (installs everything you need)
cd dashboard
./setup_dashboard.sh

# 2. Launch dashboard
./launch_dashboard.sh
```

That's it! Open your browser to `http://localhost:3838`

---

## What Gets Installed?

The setup script automatically detects and uses the **easiest option** available:

| Method | Time | Size | Sudo? | Best For |
|--------|------|------|-------|----------|
| **Conda** | ~5 min | ~1 GB | ❌ No | Most users (recommended) |
| **Docker** | ~2 min | ~500 MB | ⚠️ Yes* | Local dev |
| **Apptainer** | ~10 min | ~500 MB | ❌ No | HPC systems |

\* Docker Desktop doesn't need sudo on Windows/Mac

---

## For Your Users (GitHub Cloners)

When someone clones your repository:

```bash
git clone https://github.com/yourusername/scAnnex.git
cd scAnnex/dashboard

# First time only - setup (5-10 minutes)
./setup_dashboard.sh

# Every time - launch (instant)
./launch_dashboard.sh /path/to/results_directory
```

**That's the entire workflow!** No manual installs, no configuration files to edit.

---

## What If I Don't Have Conda/Docker/Apptainer?

The setup script will guide you! It recommends the easiest option:

1. **Conda/Mamba** (no sudo, works everywhere):
   ```bash
   curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   bash Miniforge3-Linux-x86_64.sh
   ```

2. Then run `./setup_dashboard.sh` again

---

## HPC SLURM Clusters

For running on compute nodes:

```bash
# 1. Setup once on login node
./setup_dashboard.sh

# 2. Edit launch_dashboard.slurm with your paths
vim launch_dashboard.slurm

# 3. Submit job
sbatch launch_dashboard.slurm

# 4. SSH tunnel from your local machine
ssh -N -L 8888:compute-node:8888 user@hpc.edu

# 5. Open browser: http://localhost:8888
```

---

## Pre-built Containers (For Sharing)

If you publish a release, GitHub Actions automatically builds containers:

- **Docker**: `ghcr.io/yourusername/scannex-dashboard:latest`
- **Apptainer**: Download from GitHub Releases

Users can pull these directly (no build needed):

```bash
# Docker
docker pull ghcr.io/yourusername/scannex-dashboard:latest

# Apptainer (on HPC)
wget https://github.com/yourusername/scAnnex/releases/latest/download/scannex-dashboard.sif
```

---

## Dashboard Features

Once running:

### 📊 Interactive UMAP
- Color by cell type, batch, condition, or any metadata
- Adjustable point size and transparency
- WebGL-accelerated for 100k+ cells

### 🧬 Gene Expression
- Search any gene
- Visualize expression on UMAP
- Viridis color scales

### 📈 QC Metrics
- Before/after filtering comparison
- Violin plots, scatter plots
- Threshold visualization

### 🔍 Cell Metadata
- Searchable, sortable table
- Filter by any column
- Export selections

---

## Troubleshooting

### "conda: command not found"
→ Install Miniforge (see "What If I Don't Have..." section above)

### "Port already in use"
→ Script auto-detects and uses next available port

### "Can't find .h5ad file"
→ Specify results directory: `./launch_dashboard.sh /path/to/results`

### Dashboard loads but shows no data
→ Check that your results directory contains `*annotated*.h5ad` file

---

## Architecture

```
launch_dashboard.sh
    ↓
Detects: Conda? Docker? Apptainer?
    ↓
Launches R Shiny server
    ↓
Uses Python (via reticulate) to read H5AD
    ↓
Interactive plots in browser
```

---

## File Structure

```
dashboard/
├── setup_dashboard.sh           ← Run once: sets up environment
├── launch_dashboard.sh          ← Run anytime: starts dashboard
├── launch_dashboard.slurm       ← For HPC: SLURM job script
├── environment_dashboard.yml    ← Conda dependencies
├── scannex-dashboard.def        ← Apptainer definition
├── Dockerfile                   ← Docker definition
├── app.R, ui.R, server.R       ← Shiny app code
└── global.R                     ← Functions & data loading
```

---

## For Developers

### Testing locally
```bash
conda activate scannex-dashboard
R -e "shiny::runApp('.', port=3838)"
```

### Building containers manually
```bash
# Docker
docker build -t scannex-dashboard .

# Apptainer
apptainer build scannex-dashboard.sif scannex-dashboard.def
```

### Modifying the dashboard
Edit `ui.R` (layout) and `server.R` (logic), then relaunch.

---

## Best Practices for Distribution

1. **Always provide pre-built containers** via GitHub Releases
2. **Document with screenshots** (add to wiki)
3. **Pin versions** in environment.yml for reproducibility
4. **Test on fresh system** before releasing
5. **Provide example data** for users to test with

---

## Comparison with Other Tools

| Tool | Setup | HPC | Interactive | Gene Search | Custom |
|------|-------|-----|-------------|-------------|---------|
| **scAnnex Dashboard** | 1 command | ✅ | ✅ | ✅ | ✅ |
| CellxGene | Manual install | ❌ | ✅ | ✅ | ❌ |
| UCSC Cell Browser | Complex | ⚠️ | ✅ | ✅ | ❌ |
| Loupe Browser | Download app | ❌ | ✅ | ✅ | ❌ |
| Manual (Scanpy) | Easy | ✅ | ❌ | ⚠️ | ✅ |

---

## Citation

If you use scAnnex in your research, please cite:

```
[Your Citation Here]
```

---

## Support

- 📖 **Full docs**: See `README.md` in parent directory
- 🐛 **Issues**: GitHub Issues
- 💬 **Questions**: GitHub Discussions

---

**Made with ❤️ for single-cell researchers who want tools that just work**
