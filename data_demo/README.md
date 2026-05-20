# scAnnex Demo Dataset

Synthetic PBMC dataset simulating the output of nf-core/scrnaseq (CellRanger → mtx_conversions).
Designed to exercise the full scAnnex pipeline without using patient or sensitive data.

## Dataset

| Property | Value |
|---|---|
| Samples | 6 (3 Control + 3 Treatment) |
| Cells per sample | ~1,000 |
| Genes | ~2,000 human HGNC symbols |
| Cell types | T CD4+, T CD8+, B cells, NK cells, Monocytes |
| Biological scenario | IFN-beta stimulation of PBMCs |
| DE signal | ~50 IFN-response genes upregulated in Treatment (log2FC 1.5–3.0) |
| Expected GSEA pathways | HALLMARK_INTERFERON_ALPHA_RESPONSE, HALLMARK_INTERFERON_GAMMA_RESPONSE |

## Setup

**Step 1 — Generate the H5AD files** (run once, requires the scAnnex conda environment):

```bash
# Activate the scAnnex environment first
conda activate scannex   # adjust name/path to match your installation

python data_demo/generate_demo.py
```

This creates 6 files in `data_demo/`:

```
Ctrl_PBMC_C1_filtered_matrix.h5ad
Ctrl_PBMC_C2_filtered_matrix.h5ad
Ctrl_PBMC_C3_filtered_matrix.h5ad
Treat_PBMC_T1_filtered_matrix.h5ad
Treat_PBMC_T2_filtered_matrix.h5ad
Treat_PBMC_T3_filtered_matrix.h5ad
```

**Step 2 — Run the pipeline:**

```bash
# Local (conda environment) — run from anywhere
bash data_demo/run_command.sh

# HPC / SLURM — run from anywhere
sbatch data_demo/submit.sh
```

Both scripts locate the project root automatically via `BASH_SOURCE[0]`, so they can be called from any directory.

## Expected results

`results_demo/` will contain:

- **QC reports** — per-sample cell counts, mitochondrial fraction, doublets
- **Clustering** — UMAP with cell type annotations (CellTypist: Immune_All_Low model)
- **Volcano plot** — ISG15, MX1, IFIT1, IFIT2, IFIT3, IFI44L clearly upregulated in Treatment
- **GSEA** — enrichment in interferon alpha and gamma response hallmarks

## Files

| File | Description |
|---|---|
| `generate_demo.py` | Script to synthesise the 6 H5AD files |
| `samplesheet.csv` | Pipeline input (6 samples, relative paths) |
| `contrasts.csv` | DGE contrast definition (Treatment vs Control) |
| `run_command.sh` | Local execution with conda profile |
| `submit.sh` | SLURM submission for HPC clusters |
