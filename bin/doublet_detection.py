#!/usr/bin/env python3
"""
Doublet detection using Scrublet (SLC Enhanced)
"""

import argparse
import json
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Detect doublets with Scrublet (SLC)")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output", type=str, required=True, help="Output H5AD file")
    parser.add_argument("--scores", type=str, required=True, help="Output scores CSV")
    parser.add_argument("--expected-doublet-rate", type=float, default=0.05,
                       help="Expected doublet rate")
    parser.add_argument("--remove-doublets", action="store_true",
                       help="Remove detected doublets (SLC toggle)")
    parser.add_argument("--save-attrition-log", action="store_true",
                       help="Save doublet attrition log")
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
    n_cells_initial = adata.n_obs
    print(f"Loaded: {n_cells_initial} cells x {adata.n_vars} genes")
    
    # Run Scrublet
    adata, scrub = run_scrublet(adata, args.expected_doublet_rate)
    
    # Generate plots
    print("Generating plots...")
    plot_doublet_results(adata, scrub)
    
    # Save scores
    adata.obs[['doublet_score', 'predicted_doublet']].to_csv(args.scores)
    
    # SLC: Remove doublets if requested
    n_doublets_predicted = adata.obs['predicted_doublet'].sum()
    if args.remove_doublets:
        print(f"\n🚫 Removing {n_doublets_predicted} predicted doublets...")
        adata = adata[~adata.obs['predicted_doublet']].copy()
        n_cells_final = adata.n_obs
        print(f"Cells after doublet removal: {n_cells_final}")
    else:
        print(f"\n📊 Doublets detected but NOT removed (remove-doublets=False)")
        n_cells_final = n_cells_initial
    
    # SLC: Save attrition log
    if args.save_attrition_log:
        attrition_data = {
            'step': 'Doublet Detection',
            'filter_type': 'scrublet',
            'threshold': f'expected_rate={args.expected_doublet_rate}',
            'cells_before': n_cells_initial,
            'cells_after': n_cells_final,
            'doublets_detected': int(n_doublets_predicted),
            'doublets_removed': int(n_doublets_predicted) if args.remove_doublets else 0,
            'pct_of_initial': round((n_doublets_predicted / n_cells_initial * 100), 2) if n_cells_initial > 0 else 0
        }
        
        with open('doublet_attrition.json', 'w') as f:
            json.dump(attrition_data, f, indent=2)
        print(f"✓ Saved doublet attrition log")
    
    # Save data with doublet annotations
    print(f"\nSaving results to: {args.output}")
    adata.write_h5ad(args.output)
    
    print("✓ Doublet detection complete")


if __name__ == "__main__":
    main()
