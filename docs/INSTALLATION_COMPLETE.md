# scAnnex SLC - Environment Setup Complete

## ✅ Installation Summary

All components have been successfully installed and configured for production SLC (Single-Cell LifeCycle) analysis:

### 1. Miniforge Installation
- **Location:** `/home/damo/miniforge3`
- **Version:** Conda 25.11.0, Mamba 2.4.0
- **Status:** ✅ Installed and initialized

### 2. Conda Environment
- **Name:** `scannex-minimal`
- **Python:** 3.10.19
- **Key Packages:**
  - Scanpy: 1.9.8
  - AnnData: 0.10.5.post1
  - CellTypist: 1.6.2
  - Scrublet: 0.2.3
  - HarmonyPy: 0.0.9
- **Status:** ✅ Created and verified

### 3. Project Migration
- **Original:** `/mnt/c/Users/danie/OneDrive/Documentos/CodeProjects/scAnnex` (OneDrive/WSL)
- **New Location:** `/home/damo/scAnnex` (Native Linux filesystem)
- **Status:** ✅ Migrated (42.7 MB transferred)
- **Performance:** ~20x faster I/O operations

### 4. Hardware Configuration
- **RAM:** 16GB total
- **Max Memory (Pipeline):** 8GB
- **Max CPUs:** 4
- **Profile:** `laptop` (optimized for WSL2)
- **Status:** ✅ Configured in `conf/low_memory.config`

### 5. Test Data
- **Format:** 10x MTX
- **Location:** `/home/damo/scAnnex/test_data/mtx/filtered_feature_bc_matrix/`
- **Samplesheet:** `/home/damo/scAnnex/test_data/samplesheet_slc_test.csv`
- **Status:** ✅ Available and ready

---

## 🚀 Launch Command

To run the full SLC pipeline (UNIFY_INPUT → QC → DOUBLET → PROCESSING → CELLTYPIST):

```bash
cd /home/damo/scAnnex
./run_slc_pipeline.sh
```

### Alternative: Direct Nextflow Command

```bash
cd /home/damo/scAnnex
export PATH="/home/damo/miniforge3/bin:$PATH"
source /home/damo/miniforge3/bin/activate scannex-minimal

nextflow run main.nf \
  -profile conda,laptop \
  --input test_data/samplesheet_slc_test.csv \
  --outdir results_slc_first_run \
  --run_doublet_detection true \
  --doublet_removal true \
  --run_auto_annotation true \
  --annotation_method celltypist \
  --celltypist_model 'Immune_All_Low.pkl' \
  --run_integration false \
  -with-report results_slc_first_run/pipeline_report.html \
  -with-timeline results_slc_first_run/timeline.html \
  -with-dag results_slc_first_run/dag.html \
  -resume
```

---

## 📊 Pipeline Workflow

The SLC pipeline will execute these steps:

```
1. UNIFY_INPUT
   ├─ Convert MTX → H5AD
   ├─ Add sample metadata (batch1, control)
   └─ Output: *_unified.h5ad

2. QUALITY_CONTROL
   ├─ Quantile-based filtering (10th-90th percentile)
   ├─ Mitochondrial gene filtering (<20%)
   ├─ Generate QC plots
   └─ Output: *_filtered.h5ad + QC metrics

3. DOUBLET_DETECTION
   ├─ Scrublet algorithm
   ├─ Expected doublet rate: 5%
   ├─ Remove detected doublets
   └─ Output: *_doublet_scored.h5ad + plots

4. STANDARD_PROCESSING
   ├─ Normalization (target_sum=10000)
   ├─ HVG selection (2000 genes)
   ├─ PCA (50 components)
   ├─ UMAP embedding
   ├─ Multi-resolution Leiden clustering (0.1, 0.3, 0.5, 0.7, 0.9)
   └─ Output: *_processed.h5ad + UMAP coordinates + metadata

5. AUTO_ANNOT_CELLTYPIST
   ├─ Model: Immune_All_Low.pkl
   ├─ Majority voting enabled
   ├─ Cell type prediction
   └─ Output: *_celltypist.h5ad + annotation CSV
```

---

## 🔍 Monitoring Progress

### Real-time Monitoring
```bash
# In another terminal, watch the work directory
watch -n 5 "ls -lh /home/damo/scAnnex/work/ | tail -20"
```

### Check Logs
```bash
# View Nextflow log
tail -f /home/damo/scAnnex/.nextflow.log

# View process-specific logs
find work/ -name ".command.log" | xargs tail -n 20
```

---

## 📁 Expected Output Structure

```
results_slc_first_run/
├── unified/
│   └── PBMC_TEST_unified.h5ad
├── qc/
│   ├── PBMC_TEST_filtered.h5ad
│   ├── qc_metrics.csv
│   ├── qc_before_violin.png
│   ├── qc_before_scatter.png
│   └── qc_after_violin.png
├── doublets/
│   ├── PBMC_TEST_doublet_scored.h5ad
│   ├── doublet_scores.csv
│   ├── doublet_histogram.png
│   └── doublet_umap.png
├── standard_processing_results/
│   ├── PBMC_TEST_processed.h5ad
│   ├── pca_variance.png
│   ├── clustering_multi_resolution.png
│   ├── umap_coordinates.csv
│   └── cell_metadata.csv
├── celltypist/
│   ├── PBMC_TEST_celltypist.h5ad
│   ├── celltypist_predictions.csv
│   └── celltypist_decision_matrix.csv
├── pipeline_info/
│   ├── execution_timeline.html
│   ├── execution_report.html
│   ├── execution_trace.txt
│   └── pipeline_dag.html
├── pipeline_report.html
├── timeline.html
└── dag.html
```

---

## ⚙️ Configuration Details

### Resource Limits (laptop profile)
- `max_memory`: 8 GB
- `max_cpus`: 4
- `max_time`: 8 hours

### Process-Specific Resources
- `process_single`: 1 CPU, 2 GB RAM
- `process_low`: 2 CPUs, 3 GB RAM (QC, Doublet Detection)
- `process_medium`: 3 CPUs, 4 GB RAM (Normalization)
- `process_high`: 4 CPUs, 6 GB RAM (Clustering, UMAP)

### SLC Parameters
- **QC:** Quantile-based (10th-90th percentile) + max_mito=20%
- **Doublet Rate:** 5%
- **Normalization:** Log-normalization, target_sum=10000
- **HVG:** 2000 genes
- **PCA:** 50 components
- **UMAP:** min_dist=0.5, n_neighbors=15
- **Clustering:** Leiden, resolutions=[0.1, 0.3, 0.5, 0.7, 0.9]
- **Annotation:** CellTypist (Immune_All_Low.pkl) with majority voting

---

## 🐛 Troubleshooting

### If Pipeline Fails

**Check error logs:**
```bash
cat results_slc_first_run/pipeline_info/execution_trace.txt
```

**Resume from checkpoint:**
```bash
./run_slc_pipeline.sh test_data/samplesheet_slc_test.csv results_slc_first_run resume
```

**Memory issues:**
If you see OOM (Out of Memory) errors:
1. Reduce `n_pcs` to 30: `--n_pcs 30`
2. Reduce `n_top_genes` to 1500: `--n_top_genes 1500`
3. Increase swap space (if possible)

**Conda environment issues:**
```bash
# Recreate environment
mamba env remove -n scannex-minimal
mamba env create -f /home/damo/scAnnex/env/scanpy-minimal.yml
```

---

## 📚 Next Steps After Pipeline Completion

1. **Inspect Results:**
   ```bash
   cd results_slc_first_run
   python -c "import anndata; adata = anndata.read_h5ad('celltypist/PBMC_TEST_celltypist.h5ad'); print(adata)"
   ```

2. **View Reports:**
   - Open `pipeline_report.html` in browser (via WSL: `explorer.exe pipeline_report.html`)
   - Check `timeline.html` for performance metrics
   - Review `dag.html` for workflow visualization

3. **Interactive Analysis:**
   - Use the dashboard (requires Docker)
   - Or use Jupyter notebook:
     ```bash
     conda activate scannex-minimal
     jupyter lab
     ```

4. **Production Deployment:**
   - Build custom containers (see `docs/CONTAINER_STRATEGY.md`)
   - Test on HPC with Singularity
   - Setup automated CI/CD

---

## 📞 Support

For issues or questions:
- Check `docs/` directory for detailed documentation
- Review `scAnnex_execution.todo` for implementation status
- See `docs/CONTAINER_STRATEGY.md` for production deployment

---

**Installation Date:** 2026-01-20
**Pipeline Version:** 0.1.0 (SLC)
**Environment:** WSL2 Ubuntu 18.04, 16GB RAM
