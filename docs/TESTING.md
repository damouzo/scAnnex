# Testing Guide

scAnnex includes comprehensive testing capabilities for both developers and users.

## Quick Start

### Test the Pipeline

Run a complete pipeline test with demo data:

```bash
nextflow run main.nf \
  -profile test,docker \
  --outdir test_results
```

This validates the entire pipeline from input to output using the included 1k PBMC dataset.

### Test with Your Data

Choose your input format:

**H5AD format:**
```bash
nextflow run main.nf \
  --input data_demo/H5AD/samplesheet.csv \
  --outdir results
```

**10x MTX format:**
```bash
nextflow run main.nf \
  --input data_demo/10xMTX/samplesheet.csv \
  --outdir results
```

**Seurat RDS format:**
```bash
# Generate RDS file first
cd data_demo/RDS && Rscript generate_rds.R && cd ../..

# Run pipeline
nextflow run main.nf \
  --input data_demo/RDS/samplesheet.csv \
  --outdir results
```

## Demo Data

The `data_demo/` directory contains example datasets in all supported formats:

```
data_demo/
├── H5AD/           # AnnData format
├── 10xMTX/         # 10x Genomics format
└── RDS/            # Seurat format
```

Each format includes:
- Demo data file (~1k cells, PBMC dataset)
- Sample samplesheet.csv
- Generation scripts (where applicable)

See `data_demo/README.md` for complete documentation.

## Development Testing

### Module Tests

Individual pipeline modules can be tested separately:

**Test QC and Integration:**
```bash
cd tests
./test_analytical_core.sh
```

**Quick Integration Test:**
```bash
cd tests
./test_integration_quick.sh
```

Results are saved to `tests/results/`

### Inspect Outputs

Examine H5AD file contents:

```bash
python tests/inspect_output.py path/to/output.h5ad
```

### Validate Outputs

Verify H5AD structure compliance:

```bash
python bin/validate_output.py path/to/output.h5ad
```

## CI/CD Testing

Automated tests run on every push via GitHub Actions:

- Configuration validation
- Pipeline execution tests
- Code style checks

See `.github/workflows/ci.yml` for details.

## Expected Outputs

A successful pipeline run produces:

```
results/
├── unify_input/          # Standardized input files
├── quality_control/      # QC metrics and filtered data
├── standard_processing/  # Normalized and processed data
├── normalize_integrate/  # Batch-corrected data (if enabled)
├── auto_annot/          # Cell type annotations (if enabled)
└── pipeline_info/       # Execution reports and logs
```

## Troubleshooting

**Pipeline fails:**
- Check `.nextflow.log`
- Review `results/pipeline_info/execution_trace.txt`
- Verify input data format

**Dashboard doesn't load results:**
- Ensure output H5AD contains `X_umap` in `.obsm`
- Verify required metadata columns exist
- Check dashboard logs for errors

**Out of memory:**
- Reduce `--n_top_genes` (default: 2000)
- Decrease `--n_pcs` (default: 50)
- Use `profile laptop` for smaller datasets

For more help, see `docs/Troubleshooting.md`
