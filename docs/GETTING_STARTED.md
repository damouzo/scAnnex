# Getting Started

scAnnex is a comprehensive single-cell RNA-seq analysis pipeline. Get results in minutes.

## Quick Start

### 1. Clone and Setup

```bash
git clone https://github.com/yourusername/scAnnex.git
cd scAnnex
```

### 2. Install Dependencies

**Option A: Conda (Recommended)**

```bash
mamba env create -f env/scanpy.yml
conda activate scannex
```

**Option B: Docker**

No installation needed. Use `-profile docker` when running the pipeline.

### 3. Run with Demo Data

```bash
nextflow run main.nf \
  --input data_demo/H5AD/samplesheet.csv \
  --outdir results \
  -profile conda
```

That's it. Your analysis runs automatically.

**Important:** Always use `-profile conda` (or `docker`/`singularity`) to ensure dependencies are available.

## Input Formats

scAnnex supports three formats:

| Format | Description | Example |
|--------|-------------|---------|
| **H5AD** | AnnData files | `.h5ad` |
| **10x MTX** | 10x Genomics output | `filtered_feature_bc_matrix/` |
| **RDS** | Seurat objects | `.rds` |

See `data_demo/` for examples of each format.

## Pipeline Profiles

Choose the profile that matches your system:

```bash
# Local machine with conda
nextflow run main.nf -profile conda,laptop

# Local machine with Docker
nextflow run main.nf -profile docker

# HPC with Singularity
nextflow run main.nf -profile singularity,cluster
```

## Configuration

### Basic Options

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --min_genes 200 \
  --min_cells 3
```

### Quality Control

```bash
--use_mad_thresholds    # Automatic filtering with MAD
--mad_multiplier 5.0    # MAD threshold (default: 5.0)
```

### Integration

```bash
--run_integration       # Enable batch correction
--batch_key batch       # Batch column name
--integration_method harmony
```

### Cell Type Annotation

```bash
--run_auto_annotation              # Enable CellTypist
--celltypist_model Immune_All_Low.pkl
```

For complete options, run:
```bash
nextflow run main.nf --help
```

## Creating a Samplesheet

The samplesheet defines your input data:

```csv
sample_id,file_type,file_path,batch,condition
PBMC_1,h5ad,data/sample1.h5ad,batch1,control
PBMC_2,mtx,data/sample2_mtx/,batch1,treated
```

Required columns:
- `sample_id` — Unique identifier
- `file_type` — Format (h5ad, mtx, or rds)
- `file_path` — Path to data
- `batch` — Batch identifier
- `condition` — Experimental condition

## Outputs

Results are organized by analysis stage:

```
results/
├── unify_input/          # Standardized inputs
├── quality_control/      # QC reports and filtered data
├── standard_processing/  # Normalized data
├── normalize_integrate/  # Integrated data
├── auto_annot/          # Annotated cells
└── pipeline_info/       # Execution logs
```

## View Results

Launch the interactive dashboard:

```bash
cd dashboard
./launch_dashboard.sh
```

Open your browser to `http://localhost:3838`

## Next Steps

- **Advanced configuration**: See `docs/CONFIGURATION.md`
- **Test your data**: Try different input formats in `data_demo/`
- **Use the dashboard**: Launch with `dashboard/launch_dashboard.sh`

## Requirements

**Minimum:**
- 8 GB RAM
- 4 CPU cores
- 10 GB disk space

**Recommended:**
- 16 GB RAM
- 8 CPU cores
- 50 GB disk space

**Software:**
- Nextflow ≥ 23.04
- Docker or Conda/Mamba
- (HPC only) Singularity ≥ 3.7

## Support

Need help? Check `docs/Troubleshooting.md` or open an issue on GitHub.
