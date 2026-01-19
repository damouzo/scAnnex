#!/usr/bin/env python3
"""
Doublet detection using Scrublet
"""

import argparse
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Detect doublets with Scrublet")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output", type=str, required=True, help="Output H5AD file")
    parser.add_argument("--scores", type=str, required=True, help="Output scores CSV")
    parser.add_argument("--expected-doublet-rate", type=float, default=0.05,
                       help="Expected doublet rate")
    return parser.parse_args()


def run_scrublet(adata, expected_doublet_rate=0.05):
    """Run Scrublet doublet detection"""
    print("Running Scrublet doublet detection...")
    
    scrub = scr.Scrublet(
        adata.X,
        expected_doublet_rate=expected_doublet_rate
    )
    
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=30
    )
    
    # Add results to adata
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets
    
    # Summary
    n_doublets = predicted_doublets.sum()
    pct_doublets = 100 * n_doublets / len(predicted_doublets)
    print(f"Predicted doublets: {n_doublets} ({pct_doublets:.2f}%)")
    
    return adata, scrub


def plot_doublet_results(adata, scrub):
    """Generate doublet detection plots"""
    # Histogram of doublet scores
    fig, ax = plt.subplots(figsize=(8, 5))
    scrub.plot_histogram()
    plt.savefig("doublet_histogram.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP colored by doublet score (if UMAP exists)
    if 'X_umap' in adata.obsm:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color='doublet_score', ax=ax, show=False)
        plt.savefig("doublet_umap.png", dpi=300, bbox_inches='tight')
        plt.close()


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Run Scrublet
    adata, scrub = run_scrublet(adata, args.expected_doublet_rate)
    
    # Generate plots
    print("Generating plots...")
    plot_doublet_results(adata, scrub)
    
    # Save scores
    adata.obs[['doublet_score', 'predicted_doublet']].to_csv(args.scores)
    
    # Save data with doublet annotations
    print(f"Saving results to: {args.output}")
    adata.write_h5ad(args.output)
    
    print("✓ Doublet detection complete")


if __name__ == "__main__":
    main()
