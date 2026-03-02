# Multi-Sample Demo Data for DGE Testing

This directory contains demo data for testing differential gene expression (DGE) analysis with multiple samples.

## Quick Start

### Generate Demo Data

Run the generation script (requires conda environment with scanpy):

```bash
bash data_demo/generate_multisample_demo.sh
```

This will create:
- 4 H5AD files (2 control + 2 treated samples)
- `samplesheet.csv` for pipeline input
- Samples split across 2 batches

### Run Pipeline with DGE

```bash
nextflow run main.nf \
  --input data_demo/MultiSample/samplesheet.csv \
  --outdir results_multisample \
  --run_dge \
  --contrasts_file data_demo/MultiSample/contrasts_example.csv \
  --run_integration \
  --batch_key batch \
  -profile conda \
  --max_memory 8.GB
```

## Sample Structure

| Sample ID | Batch | Condition | Cells |
|-----------|-------|-----------|-------|
| control_batch1 | batch1 | control | ~250 |
| control_batch2 | batch2 | control | ~250 |
| treated_batch1 | batch1 | treated | ~250 |
| treated_batch2 | batch2 | treated | ~250 |

Total: ~1000 cells (split from pbmc_1k.h5ad)

## Contrasts File

The `contrasts_example.csv` defines one simple contrast:

```csv
contrast_id,variable,group1,group2
treated_vs_control,condition,treated,control
```

This compares all "treated" samples vs all "control" samples.

## Advanced Contrasts

You can create more complex contrasts. See `docs/contrasts_schema.md` for examples:

- Filter by cluster/cell type
- Multiple contrasts in one file
- Batch-specific comparisons

## Expected Output

After running the pipeline, you'll find DGE results in:

```
results_multisample/
├── dge/
│   ├── dge_results/
│   │   ├── dge_results_combined.csv          # All contrasts combined
│   │   ├── treated_vs_control_all.csv        # Full results
│   │   ├── treated_vs_control_sig.csv        # Significant genes only
│   │   ├── treated_vs_control_top50.csv      # Top 50 genes
│   │   └── plots/
│   │       ├── volcano_treated_vs_control.png
│   │       └── volcano_treated_vs_control.pdf
│   └── *.h5ad                                 # Updated AnnData with DGE results
```

## Integration Metrics

When running with `--run_integration`, integration metrics will be calculated and saved in:

```
results_multisample/normalized/integration_results/
├── integration_metrics.json    # kBET, LISI, Silhouette scores
└── integration_plots/          # Visualization of integration quality
```

## Notes

- Demo data is a random split of pbmc_1k.h5ad (not real biological replicates)
- Use this for testing pipeline functionality, not for biological conclusions
- For real analysis, use proper biological replicates with batch information
