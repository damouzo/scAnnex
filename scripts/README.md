# Utility Scripts

Helper scripts for running and managing scAnnex.

## Scripts

### `run_slc_pipeline.sh`

Launch the full scAnnex pipeline with pre-configured settings.

**Usage:**
```bash
scripts/run_slc_pipeline.sh [samplesheet] [outdir] [resume]
```

**Examples:**
```bash
# With defaults
scripts/run_slc_pipeline.sh

# Custom samplesheet
scripts/run_slc_pipeline.sh data/my_samples.csv results_custom

# Resume previous run
scripts/run_slc_pipeline.sh data/my_samples.csv results_custom resume
```

**Features:**
- Activates conda environment automatically
- Verifies required packages
- Displays configuration before launch
- Uses optimized profile for local execution

### `verify_environment.sh`

Verify that your environment is correctly configured.

**Usage:**
```bash
scripts/verify_environment.sh
```

**Checks:**
- Conda/Mamba installation
- scannex environment exists
- Required Python packages (scanpy, anndata, celltypist)
- Nextflow installation
- Demo data availability

**Example output:**
```
[1/5] Checking Conda/Mamba installation...
  ✓ Mamba found: 1.5.1
[2/5] Checking scannex environment...
  ✓ Environment 'scannex' exists
[3/5] Testing package imports...
  ✓ Scanpy 1.9.3
  ✓ AnnData 0.9.2
  ✓ CellTypist 1.6.0
[4/5] Checking Nextflow...
  ✓ Nextflow found: 23.04.0
[5/5] Checking demo data...
  ✓ Samplesheet found
  ✓ MTX data found (8.9M)

✓✓✓ ALL CHECKS PASSED - READY TO RUN PIPELINE ✓✓✓
```

## When to Use

**Use `run_slc_pipeline.sh` when:**
- You want a simplified launch experience
- Running locally with conda
- Need consistent settings across runs

**Use `verify_environment.sh` when:**
- Setting up a new system
- Troubleshooting environment issues
- After updating dependencies

**For direct pipeline control, use Nextflow directly:**
```bash
nextflow run main.nf --input samplesheet.csv --outdir results
```

See `docs/GETTING_STARTED.md` for complete documentation.
