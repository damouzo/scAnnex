#!/bin/bash
# Test scAnnex pipeline with DGE and integration metrics
# Profile: Apocrita HPC with Singularity containers

set -e

echo "Testing scAnnex pipeline with DGE functionality..."
echo "Profile: apocrita + singularity"
echo ""

# Configuration
INPUT="data_demo/MultiSample/samplesheet.csv"
CONTRASTS="data_demo/MultiSample/contrasts_example.csv"
OUTDIR="/data/BCI-KRP/projects/scAnnex/results_dge_test"
WORKDIR="/gpfs/scratch/$USER/scannex_work_dge_test"

# Clean previous work directory (optional)
if [ -d "$WORKDIR" ]; then
    echo "Cleaning previous work directory: $WORKDIR"
    rm -rf "$WORKDIR"
fi

# Run pipeline
nextflow run main.nf \
    -profile apocrita,singularity \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --run_dge \
    --contrasts_file "$CONTRASTS" \
    --run_integration \
    --batch_key batch \
    --calculate_integration_metrics \
    --max_memory 16.GB \
    --max_cpus 4 \
    -w "$WORKDIR" \
    -resume

echo ""
echo "Pipeline completed!"
echo "Results in: $OUTDIR"
echo ""
echo "Check DGE results:"
echo "  - DGE tables: $OUTDIR/dge/dge_results/"
echo "  - DGE plots: $OUTDIR/dge/dge_results/plots/"
echo ""
echo "Check integration metrics:"
echo "  - Metrics: $OUTDIR/normalized/integration_results/integration_metrics.json"
echo "  - Plots: $OUTDIR/normalized/integration_results/"
