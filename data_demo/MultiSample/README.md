# Multi-Sample Public Demo Data

This directory contains scripts and examples to generate a small public multi-sample dataset
for DGE testing, based on GEO `GSE96583` (Kang IFN-beta PBMC).

The generator can create one format at a time: `H5AD`, `MTX`, or `RDS`.

## Quick Start

Generate H5AD demo (default):

```bash
python data_demo/MultiSample/download_public_multisample_demo.py \
  --format h5ad \
  --max-cells-per-condition 600
```

Generate MTX demo:

```bash
python data_demo/MultiSample/download_public_multisample_demo.py \
  --format mtx \
  --max-cells-per-condition 600
```

Generate RDS demo:

```bash
python data_demo/MultiSample/download_public_multisample_demo.py \
  --format rds \
  --max-cells-per-condition 600
```

## Output Structure

Each format is generated in its own subfolder:

- `data_demo/MultiSample/H5AD/`
- `data_demo/MultiSample/MTX/`
- `data_demo/MultiSample/RDS/`

Each subfolder contains:

- `samplesheet.csv` (format-specific)
- `contrasts_example.csv` (format-specific)
- Data files for 4 samples (`control_batch1`, `control_batch2`, `treated_batch1`, `treated_batch2`)

## Run Pipeline

H5AD:

```bash
nextflow run main.nf \
  --input data_demo/MultiSample/H5AD/samplesheet.csv \
  --contrasts_file data_demo/MultiSample/H5AD/contrasts_example.csv \
  --outdir results_multisample_h5ad \
  --run_dge \
  -profile conda
```

MTX:

```bash
nextflow run main.nf \
  --input data_demo/MultiSample/MTX/samplesheet.csv \
  --contrasts_file data_demo/MultiSample/MTX/contrasts_example.csv \
  --outdir results_multisample_mtx \
  --run_dge \
  -profile conda
```

RDS:

```bash
nextflow run main.nf \
  --input data_demo/MultiSample/RDS/samplesheet.csv \
  --contrasts_file data_demo/MultiSample/RDS/contrasts_example.csv \
  --outdir results_multisample_rds \
  --run_dge \
  -profile conda
```

## Notes

- Root-level `samplesheet_example.csv` and `contrasts_examples.csv` are templates only.
- For actual runs, use the files generated inside `H5AD/`, `MTX/`, or `RDS/`.
- The generator balances control/treated cells before splitting into 4 pseudo-replicates.
