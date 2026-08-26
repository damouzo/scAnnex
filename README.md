<div align="center">

<table>
<tr>
<td width="300" align="center">
<img src="docs/images/Logo.png" alt="scAnnex Logo" width="280"/>
</td>
<td align="center">

## **scAnnex**
<sub><span style="color:red">**Currently under active development**</span></sub>

**Automated scRNA-seq analysis.**  
From raw counts to insights.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-44A833.svg?labelColor=000000)](https://docs.conda.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[Quick Start](#quick-start) •
[Documentation](#documentation) •
[Dashboard](#dashboard) •
[Examples](#examples)

</td>
</tr>
</table>

</div>

---

## Demo

https://github.com/user-attachments/assets/7465ba9e-d781-47de-9917-c4f73c89f8a8



---

## What it does

scAnnex automates the complete workflow for single-cell RNA-seq analysis:

- **Nextflow scRNA Analaysis Pipeline** — End-to-end pipeline in nextflow dsl2
- **Interactive Dashboard** — Real-time exploration with R Shiny
- **Quality Control** — Adaptive filtering with MAD-based automatic thresholds (5 MAD) plus a hard mitochondrial cap
- **Doublet Detection** — Scrublet integration for automated doublet removal
- **Normalization** — Log-normalization and highly variable gene selection
- **Batch Correction** — Harmony integration for multi-sample datasets
- **Clustering** — Multi-resolution Leiden clustering
- **Cell Annotation** — CellTypist + Azimuth (bone marrow) + SingleR (Novershtern) + scType
- **Annotation Station** — Define cell types your way with rule-based annotation
- **Differential Gene Expression** — Pseudo-bulk DESeq2 per cell type (paper-ready, avoids pseudo-replication) plus Wilcoxon cluster markers
- **Gene Set Enrichment Analysis** — GSEA on pseudo-bulk contrast tables


## Quick Start

### Test with Demo Data

First, generate the synthetic H5AD files:

```bash
python data_demo/generate_demo.py
```

Then run the pipeline:

```bash
bash data_demo/run_command.sh
```

### Run with Your Data

```bash
nextflow run main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --outdir results \
  --max_memory '8.GB'
```

**Samplesheet format** (`samplesheet.csv`):
```csv
sample_id,file_type,file_path,batch,condition
sample1,h5ad,data/sample1.h5ad,batch1,control
sample2,h5ad,data/sample2.h5ad,batch1,treated
```

### Launch the dashboard

```bash
# 1. One-time setup (installs everything you need)
cd dashboard
./setup_dashboard.sh

# 2. Launch dashboard
./launch_dashboard.sh
```

---

## Profiles

Choose the right execution profile for your environment:

| Profile | Use Case | Command |
|---------|----------|---------|
| **singularity** | HPC clusters (recommended) | `-profile singularity` |
| **conda** | Works everywhere (fallback) | `-profile conda` |
| **apocrita** | QMUL HPC with Singularity | `-profile apocrita,singularity` |
| **slurm** | Generic SLURM HPC | `-profile slurm,singularity` |


> Memory allocation adapts automatically. Processes use 60-75% of `--max_memory` to prevent failures.

---

## Parameters

### Essential

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Samplesheet CSV path | *required* |
| `--outdir` | Output directory | `./results` |
| `--max_memory` | Maximum memory available | `128.GB` |

### Quality Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--use_mad_thresholds` | Use MAD-based automatic thresholds | `true` |
| `--mad_multiplier` | MAD multiplier | `5.0` |
| `--use_quantile_filtering` | Use quantile-based thresholds (alternative) | `false` |
| `--feature_quantile_low` | Lower percentile for genes | `0.10` |
| `--feature_quantile_high` | Upper percentile for genes | `0.90` |
| `--max_mito_percent` | Max mitochondrial percentage (hard cap) | `10` |

### Processing

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--run_doublet_detection` | Detect and remove doublets | `true` |
| `--run_integration` | Enable batch correction | `false` |
| `--batch_key` | Batch column name | `null` |
| `--run_auto_annotation` | Annotate cell types | `true` |
| `--celltypist_models` | CellTypist models | `Immune_All_Low.pkl,Adult_cHSPCs_Illumina.pkl` |
| `--azimuth_refs` | Azimuth reference (bone marrow) | `bonemarrowref` |
| `--singler_refs` | SingleR reference | `NovershternHematopoieticData` |

### Pseudobulk DGE (DESeq2, per cell type)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--run_pseudobulk_dge` | Run pseudo-bulk DESeq2 DGE per cell type | `true` |
| `--pseudobulk_min_cells` | Min cells per sample x cell type aggregate | `10` |
| `--pseudobulk_groupby` | obs column for aggregates | `cell_type` |
| `--pseudobulk_control_group` | Reference condition (defaults to `dge_reference`) | `null` |
| `--pseudobulk_padj_cutoff` | Adjusted p-value cutoff | `0.05` |
| `--pseudobulk_shrink` | LFC shrinkage type | `apeglm` |

### Clustering

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--clustering_resolutions` | Leiden resolutions | `0.1,0.3,0.5,0.7,0.9` |
| `--n_neighbors` | Number of neighbors for UMAP | `15` |
| `--n_pcs` | Principal components | `50` |

---

# Dashboard

Explore your data in real time. No coding required.

**Quality Control**
- Visualize metrics across all samples
- Identify filtering thresholds
- Track cell attrition

**Clustering & UMAP**
- Interactive scatter plots with WebGL
- Color by any metadata or clustering
- Zoom, pan, export

**Gene Expression**
- Search any gene
- Overlay expression on UMAP
- See distribution instantly

**Annotation Station**
- Define cell types with simple rules
- Combine clustering and auto-annotations
- Preview changes in real time
- Save annotations directly to H5AD

---

## Examples

### Demo: 6-sample IFN-beta stimulation (PBMC)

Synthetic dataset: 3 Control + 3 Treatment samples, batch-corrected with Harmony, DGE and GSEA included.

```bash
# Generate synthetic H5AD files
python data_demo/generate_demo.py

# Run full pipeline (pseudo-bulk DESeq2 DGE + GSEA included)
nextflow run main.nf \
  -profile conda \
  --input data_demo/samplesheet.csv \
  --contrasts_file data_demo/contrasts.csv \
  --outdir results_demo \
  --run_integration true \
  --batch_key sample_id \
  --dge_groupby condition \
  --dge_reference Control \
  --run_gsea true \
  --organism human \
  --max_memory '8.GB'

# Launch dashboard
cd dashboard && bash launch_dashboard.sh
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.



## Citation


```
scAnnex: Automated Nextflow pipeline for single-cell RNA-seq analysis. https://github.com/damouzo/scAnnex
```


