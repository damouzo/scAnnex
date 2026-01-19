# scAnnex Test Data

This directory contains test data and scripts for validating the scAnnex UNIFY_INPUT module.

## Directory Structure

```
test_data/
├── mtx/                              # 10x MTX format (PBMC 1k)
│   └── filtered_feature_bc_matrix/
├── h5ad/                             # H5AD format (to be generated)
├── rds/                              # RDS/Seurat format (to be generated)
├── outputs/                          # Test outputs
├── samplesheet_test.csv             # Test samplesheet
├── create_test_files.py             # Generate H5AD test file
├── create_test_rds.R                # Generate RDS test file
├── validate_output.py               # Validate output structure
└── run_tests.sh                     # Automated test runner
```

## Quick Start

### Option 1: Automated Testing (Requires Docker)

```bash
cd test_data
./run_tests.sh
```

This will:
1. Create test H5AD and RDS files using containers
2. Run UNIFY_INPUT on all three formats (MTX, H5AD, RDS)
3. Validate outputs against scAnnex conventions

### Option 2: Manual Testing

#### Step 1: Create Test Files

**H5AD file (requires Python with scanpy):**
```bash
cd test_data
python create_test_files.py
```

**RDS file (requires R with Seurat + SeuratDisk):**
```bash
cd test_data
Rscript create_test_rds.R
```

#### Step 2: Test Individual Conversions

**Test MTX → H5AD:**
```bash
python ../bin/unify_input.py \
    --input test_data/mtx/filtered_feature_bc_matrix \
    --input-type mtx \
    --output test_data/outputs/PBMC_MTX_unified.h5ad \
    --sample-id PBMC_MTX \
    --batch batch1 \
    --condition control
```

**Test H5AD → H5AD (with metadata integration):**
```bash
python ../bin/unify_input.py \
    --input test_data/h5ad/pbmc_100cells.h5ad \
    --input-type h5ad \
    --output test_data/outputs/PBMC_H5AD_unified.h5ad \
    --sample-id PBMC_H5AD \
    --batch batch1 \
    --condition treated
```

**Test RDS → H5AD (requires R with Seurat/SeuratDisk):**
```bash
python ../bin/unify_input.py \
    --input test_data/rds/pbmc_seurat.rds \
    --input-type rds \
    --output test_data/outputs/PBMC_RDS_unified.h5ad \
    --sample-id PBMC_RDS \
    --batch batch2 \
    --condition control
```

#### Step 3: Validate Output

```bash
python test_data/validate_output.py test_data/outputs/PBMC_MTX_unified.h5ad
```

Expected validation output:
- ✓ Required .obs columns: sample_id, batch, condition
- ✓ 'counts' layer present with raw counts
- ✓ scannex metadata in .uns
- ✓ No duplicate gene names

## Test Data Sources

### MTX Data
- **Source**: 10x Genomics PBMC 1k v3
- **URL**: https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/
- **Format**: filtered_feature_bc_matrix (MTX)
- **Cells**: ~1,000 PBMCs
- **Downloaded**: Automatically by test setup

### H5AD Test File
- **Generated from**: MTX PBMC 1k subset (100 cells)
- **Purpose**: Test H5AD input and metadata integration
- **Created by**: `create_test_files.py`

### RDS Test File
- **Generated from**: MTX PBMC 1k data
- **Format**: Seurat object (v4 or v5)
- **Purpose**: Test RDS→h5seurat→h5ad conversion pipeline
- **Created by**: `create_test_rds.R`

## Validation Checklist

For each unified output file, verify:

1. **File Creation**
   - [ ] Output .h5ad file created
   - [ ] File size reasonable (>1 MB for full PBMC)

2. **Metadata in .obs**
   - [ ] `sample_id` column present
   - [ ] `batch` column present
   - [ ] `condition` column present

3. **Data Layers**
   - [ ] `counts` layer contains raw counts
   - [ ] `.X` matrix present

4. **Conversion Metadata (.uns)**
   - [ ] `scannex.conversion` dict present
   - [ ] `input_type`, `sample_id`, `timestamp` recorded
   - [ ] `scannex.versions` dict with package versions

5. **Gene/Cell Indexing**
   - [ ] No duplicate gene names (or unique suffixes)
   - [ ] Cell barcodes preserved or generated

## Troubleshooting

### Python Package Errors
```bash
pip install scanpy numpy pandas anndata
```

### R Package Errors
```R
install.packages(c("Seurat", "remotes"))
remotes::install_github("mojaveazure/seurat-disk")
```

### RDS Conversion Fails
- Check R script output in bin/convert_rds_to_h5seurat.R
- Verify Seurat object structure
- Check SeuratDisk installation

### Memory Issues
- Reduce dataset size for testing
- Use backed mode for large H5AD files

## Next Steps After Testing

Once UNIFY_INPUT tests pass:

1. **Update scAnnex_execution.todo** with test results
2. **Commit test data and outputs** (optional, exclude large files)
3. **Proceed to Dashboard Implementation** (Phase 8)
4. **Test full pipeline** with Nextflow test profile

## Nextflow Pipeline Test

To test within the Nextflow pipeline:

```bash
nextflow run main.nf \
    --input test_data/samplesheet_test.csv \
    --outdir results_test \
    -profile docker,test
```

This will process all three input formats through the complete pipeline.
