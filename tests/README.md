# scAnnex Tests

Development test scripts for validating individual pipeline modules.

## Structure

```
tests/
├── test_analytical_core.sh       # QC + Integration module tests
├── test_integration_quick.sh     # Quick integration test
├── inspect_output.py             # H5AD inspection utility
├── samplesheet_slc_test.csv      # Legacy test samplesheet
└── README.md                     # This file
```

## Purpose

These tests are for **pipeline developers** to validate core functionality during development.

**Users testing the pipeline should use `data_demo/` instead.**

## Running Tests

### Prerequisites

- Docker installed and running
- scAnnex environment configured

### Module Tests

**Test QC and Integration modules:**

```bash
cd tests
./test_analytical_core.sh
```

Results saved to `tests/results/analytical_core/`

**Quick integration test:**

```bash
cd tests
./test_integration_quick.sh
```

### Inspect H5AD Files

```bash
python tests/inspect_output.py path/to/file.h5ad
```

## For Users

For pipeline testing and examples, see `data_demo/`:

```bash
# Run pipeline with demo data
nextflow run main.nf --input data_demo/H5AD/samplesheet.csv --outdir results
```

See `data_demo/README.md` for complete documentation.

## CI/CD Integration

These tests can be integrated with GitHub Actions. See `.github/workflows/ci.yml` for automated testing configuration.
