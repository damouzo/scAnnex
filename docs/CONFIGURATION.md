# Pipeline Configuration

Advanced configuration options for scAnnex.

## Quality Control

### Quantile-Based Filtering (Recommended)

Automatically filter cells based on data distribution:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --use_quantile_filtering true \
  --feature_quantile_low 0.10 \
  --feature_quantile_high 0.90
```

**Parameters:**
- `--use_quantile_filtering` - Enable quantile filtering (default: true)
- `--feature_quantile_low` - Lower percentile for features (default: 0.10)
- `--feature_quantile_high` - Upper percentile for features (default: 0.90)
- `--count_quantile_low` - Lower percentile for counts (default: 0.10)
- `--count_quantile_high` - Upper percentile for counts (default: 0.90)

### Hard Thresholds

Set explicit cutoffs:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --min_genes 200 \
  --min_cells 3 \
  --max_mito_percent 20
```

### MAD-Based Filtering

Automatic outlier detection using median absolute deviation:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --use_mad_thresholds true \
  --mad_multiplier 5.0
```

## Batch Correction

Integrate multiple samples or batches:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --run_integration \
  --batch_key batch \
  --integration_method harmony \
  --harmony_theta 2.0
```

**Integration methods:**
- `harmony` - Harmony algorithm (recommended)
- `scanorama` - Scanorama integration
- `bbknn` - Batch-balanced k-nearest neighbors

## Clustering

### Multi-Resolution Clustering

Test multiple resolutions automatically:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --clustering_resolutions '0.1,0.3,0.5,0.7,0.9' \
  --default_clustering_resolution 0.5
```

### Clustering Algorithm

Choose between Leiden or Louvain:

```bash
--clustering_method leiden    # Default, recommended
--clustering_method louvain   # Alternative
```

## Cell Type Annotation

### CellTypist

Automatic annotation with pre-trained models:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --run_auto_annotation \
  --celltypist_model Immune_All_Low.pkl \
  --celltypist_majority_voting true
```

**Available models:**
- `Immune_All_Low.pkl` - Pan-immune, low resolution
- `Immune_All_High.pkl` - Pan-immune, high resolution
- Custom models from CellTypist repository

### Marker-Based Annotation

Use custom marker lists:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --annotation_method marker_based \
  --marker_list markers.csv
```

## Resource Management

### Memory Limits

Set maximum memory usage:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --max_memory '8.GB'   # For laptops
  --max_memory '128.GB' # For servers
```

### CPU Limits

Control CPU allocation:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --max_cpus 4   # For laptops
  --max_cpus 16  # For workstations
```

### Combined Profiles

Use preset configurations:

```bash
# Laptop profile (low memory)
nextflow run main.nf -profile conda,laptop --max_memory '8.GB'

# HPC profile (high resources)
nextflow run main.nf -profile singularity,cluster
```

## Interactive Dashboard

### Configuration

Control dashboard behavior:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --enable_dashboard true \
  --dashboard_port 3838 \
  --dashboard_host localhost
```

### Launching

After pipeline completion:

```bash
cd dashboard
./launch_dashboard.sh
```

Access at `http://localhost:3838`

## Execution Profiles

### Local Execution

**With Conda:**
```bash
nextflow run main.nf -profile conda --input samplesheet.csv
```

**With Docker:**
```bash
nextflow run main.nf -profile docker --input samplesheet.csv
```

### HPC Execution

**With Singularity:**
```bash
nextflow run main.nf -profile singularity --input samplesheet.csv
```

**With Slurm:**
```bash
nextflow run main.nf \
  -profile singularity \
  --input samplesheet.csv \
  -c conf/slurm.config
```

## Output Organization

Configure output structure:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir my_results \
  --publish_dir_mode copy    # or 'symlink', 'link'
```

## Troubleshooting

### Pipeline Fails with Memory Error

Reduce memory requirements:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --max_memory '8.GB' \
  -profile conda,laptop
```

### Too Slow on Laptop

Use the laptop profile:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  -profile conda,laptop
```

### Containers Not Found

Use conda instead:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  -profile conda
```

## Complete Parameter List

View all available options:

```bash
nextflow run main.nf --help
```

Or check `nextflow_schema.json` for parameter descriptions.
