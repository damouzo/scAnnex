#!/usr/bin/env python3
"""
Merge Multiple AnnData Samples

Combines multiple H5AD files into a single AnnData object,
preserving metadata from all samples.

Author: scAnnex Development Team
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List

import scanpy as sc
import anndata as ad

# Enable writing of nullable strings
if hasattr(ad, 'settings') and hasattr(ad.settings, 'allow_write_nullable_strings'):
    ad.settings.allow_write_nullable_strings = True

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge multiple H5AD samples into a single AnnData object"
    )
    parser.add_argument(
        "--inputs",
        type=str,
        nargs='+',
        required=True,
        help="Input H5AD files to merge"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output merged H5AD file"
    )
    return parser.parse_args()


def merge_samples(input_files: List[str], output_file: str):
    """Merge multiple AnnData objects.
    
    Args:
        input_files: List of H5AD file paths
        output_file: Output file path
    """
    logger.info(f"Merging {len(input_files)} H5AD files...")
    
    # Load all samples
    adatas = []
    for i, file_path in enumerate(input_files):
        logger.info(f"  Loading sample {i+1}/{len(input_files)}: {file_path}")
        adata = sc.read_h5ad(file_path)
        logger.info(f"    {adata.n_obs} cells x {adata.n_vars} genes")
        adatas.append(adata)
    
    # Build per-sample keys for unique barcode prefixing.
    # 10x barcodes are identical across samples; without a key the merged
    # object has duplicate obs_names, which causes failures in downstream
    # R tools (read.csv / Seurat).  Using keys + index_unique produces
    # barcodes like  "HealthyDonor_MNC_H2_AAACCTGAGAAACCAT-1"  that are
    # unique AND traceable back to their source sample.
    sample_keys = []
    for i, (file_path, adata) in enumerate(zip(input_files, adatas)):
        if 'sample_id' in adata.obs.columns:
            key = str(adata.obs['sample_id'].iloc[0])
        else:
            key = Path(file_path).stem
        sample_keys.append(key)
        logger.info(f"    Barcode prefix for sample {i+1}: {key}")

    # Concatenate all samples
    logger.info("Concatenating samples...")
    merged = ad.concat(
        adatas,
        keys=sample_keys,
        index_unique="_",
        axis=0,
        join='outer',
        merge='unique',
        uns_merge='unique',
        fill_value=0
    )
    
    logger.info(f"Merged dataset: {merged.n_obs} cells x {merged.n_vars} genes")

    # Save merged dataset
    logger.info(f"Saving merged dataset to: {output_file}")
    merged.write_h5ad(output_file, compression='gzip')
    
    logger.info("✓ Merge complete")


def main():
    """Main entry point."""
    args = parse_args()
    
    try:
        merge_samples(args.inputs, args.output)
    except Exception as e:
        logger.error(f"Error during merge: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
