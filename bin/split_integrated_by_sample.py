#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc
import anndata as ad

if hasattr(ad, 'settings') and hasattr(ad.settings, 'allow_write_nullable_strings'):
    ad.settings.allow_write_nullable_strings = True


def parse_args():
    parser = argparse.ArgumentParser(description="Split integrated H5AD by sample key")
    parser.add_argument("--input", required=True, help="Integrated H5AD path")
    parser.add_argument("--sample-key", default="sample_id", help="obs column with sample IDs")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)

    if args.sample_key not in adata.obs.columns:
        return

    sample_ids = sorted(adata.obs[args.sample_key].astype(str).unique())
    for sample_id in sample_ids:
        sample_mask = adata.obs[args.sample_key].astype(str) == sample_id
        sample_adata = adata[sample_mask].copy()
        sample_adata.write_h5ad(output_dir / f"{sample_id}_integrated.h5ad", compression='gzip')


if __name__ == "__main__":
    main()
