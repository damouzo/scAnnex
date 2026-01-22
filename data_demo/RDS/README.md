# RDS Demo Data

This directory should contain a Seurat RDS file for testing the scAnnex pipeline.

## Quick Setup

Run the generation script from the `data_demo/` directory:

```bash
cd data_demo
Rscript generate_rds.R
```

This will create `pbmc_1k.rds` from the 10xMTX data.

## Requirements

- R (≥ 4.0)
- Seurat package

Install Seurat in R:
```R
install.packages('Seurat')
```

## File Details

- **File**: pbmc_1k.rds
- **Format**: Seurat object (RDS)
- **Source**: PBMC 1k v3 from 10x Genomics
- **Cells**: ~1,000 PBMCs
- **Metadata**: sample_id, batch, condition, percent.mt

## Usage

Use the provided `samplesheet.csv` to run the scAnnex pipeline:

```bash
nextflow run main.nf --input data_demo/RDS/samplesheet.csv --outdir results_rds
```

## Note

If you don't have R/Seurat installed, you can:
1. Use the H5AD or 10xMTX demo data instead
2. Generate the RDS file on another system and copy it here
3. Skip RDS testing if not using this input format
