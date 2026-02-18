#!/usr/bin/env python3
"""
Trajectory inference and pseudotime analysis
"""

import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

# Enable writing of nullable strings (required for anndata >= 0.11)
# Only set if available (older versions don't have this attribute)
if hasattr(ad, 'settings') and hasattr(ad.settings, 'allow_write_nullable_strings'):
    ad.settings.allow_write_nullable_strings = True


def parse_args():
    parser = argparse.ArgumentParser(description="Trajectory inference analysis")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output", type=str, required=True, help="Output H5AD file")
    parser.add_argument("--root-cell", type=str, default=None,
                       help="Root cell index for pseudotime calculation")
    parser.add_argument("--method", type=str, default="paga",
                       choices=["paga", "dpt"],
                       help="Trajectory method")
    return parser.parse_args()


def run_paga(adata):
    """Run PAGA trajectory inference"""
    print("Running PAGA trajectory inference...")
    
    # PAGA requires clustering
    if 'leiden' not in adata.obs.columns:
        print("Running Leiden clustering first...")
        sc.tl.leiden(adata, resolution=1.0)
    
    # Run PAGA
    sc.tl.paga(adata, groups='leiden')
    
    # Plot PAGA graph
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.paga(adata, ax=ax, show=False)
    plt.savefig("paga_graph.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Initialize PAGA positions for UMAP
    sc.pl.paga(adata, show=False)
    sc.tl.umap(adata, init_pos='paga')
    
    return adata


def run_diffusion_pseudotime(adata, root_cell=None):
    """Run diffusion pseudotime analysis"""
    print("Running diffusion pseudotime...")
    
    if root_cell is None:
        # Use first cell as root
        root_cell = 0
        print(f"Using cell {root_cell} as root")
    
    # Set root cell
    adata.uns['iroot'] = int(root_cell)
    
    # Compute DPT
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    
    return adata


def plot_trajectory_results(adata):
    """Generate trajectory visualization plots"""
    # UMAP colored by pseudotime
    if 'dpt_pseudotime' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color='dpt_pseudotime', ax=ax, show=False)
        plt.savefig("umap_pseudotime.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # PAGA + UMAP overlay
    if 'paga' in adata.uns:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.paga_compare(adata, ax=ax, show=False)
        plt.savefig("paga_umap_overlay.png", dpi=300, bbox_inches='tight')
        plt.close()


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Run trajectory analysis
    if args.method == "paga":
        adata = run_paga(adata)
    
    # Run diffusion pseudotime
    adata = run_diffusion_pseudotime(adata, args.root_cell)
    
    # Generate plots
    print("Generating plots...")
    plot_trajectory_results(adata)
    
    # Save results
    print(f"Saving results to: {args.output}")
    adata.write_h5ad(args.output)
    
    print("✓ Trajectory analysis complete")


if __name__ == "__main__":
    main()
