# scAnnex SLC Quick Start

## 🚀 Run the SLC Pipeline in 3 Steps

### Step 1: Download Test Data
```bash
python bin/download_test_data.py --output-dir test_data
```

### Step 2: Run the Pipeline
```bash
nextflow run main.nf \
  --input test_data/samplesheet.csv \
  --outdir results/ \
  -profile docker
```

### Step 3: Explore Results
```bash
# Cell Attrition Log
cat results/qc_results/cell_attrition_log.txt

# UMAP Coordinates (for dashboard)
head results/standard_processing_results/umap_coordinates.csv

# Multi-resolution clustering plots
ls results/standard_processing_results/*.png
```

---

## 📊 SLC Key Features

✅ **Quantile-based QC filtering** (10th-90th percentile)  
✅ **Cell Attrition Log** - Track where cells go  
✅ **Multi-resolution clustering** (0.1, 0.3, 0.5, 0.7, 0.9)  
✅ **CellTypist auto-annotation** - Immune_All_Low model  
✅ **Dashboard-ready outputs** - CSV files optimized for Shiny  

---

## 🎯 SLC Configuration

Edit `nextflow.config` to customize:

```groovy
params {
    // QC - Quantile filtering
    use_quantile_filtering     = true
    feature_quantile_low       = 0.10
    feature_quantile_high      = 0.90
    
    // Doublet Detection
    run_doublet_detection      = true
    doublet_removal            = true
    
    // Clustering resolutions
    clustering_resolutions     = '0.1,0.3,0.5,0.7,0.9'
    
    // Auto-annotation
    run_auto_annotation        = true
    celltypist_model           = 'Immune_All_Low.pkl'
}
```

---

## 📁 SLC Outputs

```
results/
├── qc_results/
│   ├── cell_attrition_log.csv       # Cell filtering breakdown
│   ├── cell_attrition_log.txt       # Human-readable log
│   ├── qc_before_violin.png         # Pre-filtering QC
│   └── qc_after_violin.png          # Post-filtering QC
├── standard_processing_results/
│   ├── umap_coordinates.csv         # Dashboard-ready UMAP
│   ├── cell_metadata.csv            # All annotations
│   ├── clustering_multi_resolution.png
│   └── umap_leiden_res*.png         # Per-resolution plots
└── final.h5ad                       # Complete processed data
```

---

## 🎨 Next: Dashboard

The SLC pipeline generates dashboard-ready outputs:

```r
# In Shiny dashboard
umap_data <- read.csv("results/standard_processing_results/umap_coordinates.csv")
metadata <- read.csv("results/standard_processing_results/cell_metadata.csv")

# Interactive UMAP
plot_ly(umap_data, x = ~UMAP_1, y = ~UMAP_2, color = ~leiden_0.5)
```

---

## 📖 Full Documentation

See `docs/summary_2026-01-20.md` for complete SLC implementation details.

---

**Strategy:** SLC (Simple, Lovable, Complete)  
**Version:** v1.0  
**Date:** January 20, 2026
