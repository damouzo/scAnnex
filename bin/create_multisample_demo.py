#!/usr/bin/env python3
"""
Create multi-sample demo data from existing PBMC 1k dataset.
Splits into 4 samples: 2 control + 2 treated, across 2 batches.
"""

import argparse
import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Enable writing of nullable strings (required for anndata >= 0.11)
ad.settings.allow_write_nullable_strings = True

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Generate multi-sample demo data for DGE testing")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file (e.g., pbmc_1k.h5ad)")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory for demo samples")
    parser.add_argument("--n-samples", type=int, default=4, help="Number of samples to create (default: 4)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")

    args = parser.parse_args()

    # Set random seed
    np.random.seed(args.seed)

    # Load input data
    logger.info(f"Loading input data: {args.input}")
    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define sample metadata
    sample_metadata = [
        {"sample_id": "control_batch1", "batch": "batch1", "condition": "control"},
        {"sample_id": "control_batch2", "batch": "batch2", "condition": "control"},
        {"sample_id": "treated_batch1", "batch": "batch1", "condition": "treated"},
        {"sample_id": "treated_batch2", "batch": "batch2", "condition": "treated"},
    ]

    # Shuffle cells and split into equal parts
    n_cells_per_sample = adata.n_obs // args.n_samples
    cell_indices = np.random.permutation(adata.n_obs)

    # Create samplesheet rows
    samplesheet_rows = []

    for i, meta in enumerate(sample_metadata[: args.n_samples]):
        # Get cell indices for this sample
        start_idx = i * n_cells_per_sample
        if i == args.n_samples - 1:
            # Last sample gets remaining cells
            end_idx = adata.n_obs
        else:
            end_idx = (i + 1) * n_cells_per_sample

        sample_indices = cell_indices[start_idx:end_idx]

        # Subset AnnData
        adata_subset = adata[sample_indices, :].copy()

        # Add metadata
        adata_subset.obs["sample_id"] = meta["sample_id"]
        adata_subset.obs["batch"] = meta["batch"]
        adata_subset.obs["condition"] = meta["condition"]

        # Save to file
        output_file = output_dir / f"{meta['sample_id']}.h5ad"
        logger.info(f"Saving {adata_subset.n_obs} cells to {output_file}")
        adata_subset.write_h5ad(output_file)

        # Add to samplesheet
        samplesheet_rows.append(
            {
                "sample_id": meta["sample_id"],
                "file_type": "h5ad",
                "file_path": str(output_file),
                "batch": meta["batch"],
                "condition": meta["condition"],
            }
        )

    # Create samplesheet
    samplesheet_df = pd.DataFrame(samplesheet_rows)
    samplesheet_file = output_dir / "samplesheet.csv"
    logger.info(f"Writing samplesheet to {samplesheet_file}")
    samplesheet_df.to_csv(samplesheet_file, index=False)

    logger.info("Multi-sample demo data creation complete")
    logger.info(f"Created {args.n_samples} samples in {output_dir}")


if __name__ == "__main__":
    main()
