#!/usr/bin/env python3
"""
Quality control filtering and metrics calculation for scRNA-seq data
"""

import argparse
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform quality control on scRNA-seq data"
    )
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output", type=str, required=True, help="Output H5AD file")
    parser.add_argument("--metrics", type=str, required=True, help="Output metrics CSV")
    parser.add_argument("--min-genes", type=int, default=200, help="Minimum genes per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum cells per gene")
    parser.add_argument("--max-genes", type=int, default=None, help="Maximum genes per cell")
    parser.add_argument("--max-counts", type=int, default=None, help="Maximum counts per cell")
    parser.add_argument("--max-mito", type=float, default=20, help="Maximum mitochondrial %")
    return parser.parse_args()


def calculate_qc_metrics(adata):
    """Calculate QC metrics including mitochondrial percentage"""
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    return adata


def plot_qc_metrics(adata, prefix="qc"):
    """Generate QC violin plots"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Plot total counts
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes, show=False)
    
    plt.tight_layout()
    plt.savefig(f"{prefix}_violin.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Scatter plot: counts vs genes
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', 
                  color='pct_counts_mt', ax=ax, show=False)
    plt.savefig(f"{prefix}_scatter.png", dpi=300, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Initial: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Calculate QC metrics
    print("Calculating QC metrics...")
    adata = calculate_qc_metrics(adata)
    
    # Plot before filtering
    print("Generating QC plots...")
    plot_qc_metrics(adata, prefix="qc_before")
    
    # Apply filters
    print("Applying filters...")
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    if args.max_genes:
        adata = adata[adata.obs.n_genes_by_counts < args.max_genes, :]
    if args.max_counts:
        adata = adata[adata.obs.total_counts < args.max_counts, :]
    adata = adata[adata.obs.pct_counts_mt < args.max_mito, :]
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    print(f"After filtering: {n_cells_after} cells x {n_genes_after} genes")
    print(f"Removed: {n_cells_before - n_cells_after} cells, {n_genes_before - n_genes_after} genes")
    
    # Save metrics
    metrics_df = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].describe()
    metrics_df.to_csv(args.metrics)
    
    # Save filtered data
    print(f"Saving filtered data to: {args.output}")
    adata.write_h5ad(args.output)
    
    print("✓ Quality control complete")


if __name__ == "__main__":
    main()
