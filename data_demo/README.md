# scAnnex Demo Data

This directory contains demo datasets in multiple formats for testing the scAnnex pipeline.

## Dataset Details

- **Source**: PBMC 1k v3 from 10x Genomics
- **Cells**: ~1,000 PBMCs (peripheral blood mononuclear cells)
- **Species**: Human
- **Technology**: 10x Chromium v3

## Quick Start

### Test with H5AD format

```bash
nextflow run main.nf \
  --input data_demo/H5AD/samplesheet.csv \
  --outdir results_h5ad
```

### Test with 10xMTX format

```bash
nextflow run main.nf \
  --input data_demo/10xMTX/samplesheet.csv \
  --outdir results_mtx
```

### Test with RDS format

First, generate the RDS file (requires R + Seurat):
```bash
cd data_demo
Rscript generate_rds.R
cd ..
```

Then run the pipeline:
```bash
nextflow run main.nf \
  --input data_demo/RDS/samplesheet.csv \
  --outdir results_rds
```

## Samplesheet Format

Each format includes a `samplesheet.csv` with the following structure:

```csv
sample_id,file_type,file_path,batch,condition
PBMC_1k,h5ad,data_demo/H5AD/pbmc_1k.h5ad,batch1,control
```

Required columns:
- **sample_id**: Unique identifier for the sample
- **file_type**: Format type (h5ad, mtx, or rds)
- **file_path**: Path to data file/directory
- **batch**: Batch identifier for integration
- **condition**: Experimental condition

## Regenerating Demo Files

### H5AD from MTX
```bash
cd data_demo
python3 generate_h5ad.py
```

### RDS from MTX
```bash
cd data_demo
Rscript generate_rds.R
```

Requirements:
- **H5AD**: Python 3.8+, scanpy, pandas
- **RDS**: R 4.0+, Seurat

## Pipeline Options

Common parameters:
```bash
nextflow run main.nf \
  --input <samplesheet.csv> \
  --outdir <output_directory> \
  --skip_doublet_detection \
  --skip_integration \
  --min_genes 200 \
  --min_cells 3
```

For full options, see:
```bash
nextflow run main.nf --help
```

## Expected Output

The pipeline will generate:
- Quality control reports
- Normalized and integrated data
- Cell type annotations (if enabled)
- Interactive dashboard
- Analysis reports

Output structure:
```
results/
├── pipeline_info/
├── unify_input/
├── quality_control/
├── standard_processing/
├── normalize_integrate/
├── auto_annot/
└── dashboard/
```

## Notes

- The 10xMTX format is the original data format
- H5AD is generated from 10xMTX
- RDS must be generated separately (requires R/Seurat)
- All three formats contain the same cells and should produce equivalent results
