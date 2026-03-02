#!/bin/bash
# Generate multi-sample demo data for DGE testing
# This script uses the scannex-dashboard conda environment

set -e

echo "Generating multi-sample demo data..."

# Use the scannex-dashboard conda environment Python
PYTHON_BIN=~/.conda/envs/scannex-dashboard/bin/python

# Check if conda env exists
if [ ! -f "$PYTHON_BIN" ]; then
    echo "ERROR: scannex-dashboard conda environment not found"
    echo "Please create it first or update PYTHON_BIN in this script"
    exit 1
fi

# Create output directory
mkdir -p data_demo/MultiSample/H5AD

# Generate base H5AD if it doesn't exist
if [ ! -f data_demo/H5AD/pbmc_1k.h5ad ]; then
    echo "Generating base H5AD file from MTX data..."
    $PYTHON_BIN data_demo/H5AD/generate_h5ad.py
fi

# Run the Python script to generate samples
$PYTHON_BIN bin/create_multisample_demo.py \
    --input data_demo/H5AD/pbmc_1k.h5ad \
    --output-dir data_demo/MultiSample/H5AD \
    --n-samples 4 \
    --seed 42

echo ""
echo "Multi-sample demo data created successfully!"
echo "Location: data_demo/MultiSample/H5AD/"
echo "Samplesheet: data_demo/MultiSample/H5AD/samplesheet.csv"
echo "Contrasts: data_demo/MultiSample/contrasts_example.csv"
