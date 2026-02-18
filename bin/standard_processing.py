#!/usr/bin/env python3
"""
Standard Processing - SLC Pipeline
Scanpy standard workflow with multi-resolution clustering

Features:
- Normalize & Log1p
- HVG selection
- Scale
- PCA
- Neighbors
- UMAP
- Multi-resolution Leiden clustering (0.1, 0.3, 0.5, 0.7, 0.9)

Author: scAnnex SLC Development Team
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
import anndata as ad

# Enable writing of nullable strings (required for anndata >= 0.11)
ad.settings.allow_write_nullable_strings = True

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Standard Scanpy processing with multi-resolution clustering (SLC)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input H5AD file (post-QC, post-doublet)"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output H5AD file with clustering"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="standard_processing_results",
        help="Output directory for plots"
    )
    
    # Normalization
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Target sum for normalization (default: 10000)"
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Number of highly variable genes (default: 2000)"
    )
    
    # PCA
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of principal components (default: 50)"
    )
    
    # Neighbors
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for graph construction (default: 15)"
    )
    
    # UMAP
    parser.add_argument(
        "--umap-min-dist",
        type=float,
        default=0.5,
        help="UMAP min_dist parameter (default: 0.5)"
    )
    
    # Clustering
    parser.add_argument(
        "--clustering-method",
        type=str,
        default="leiden",
        choices=["leiden", "louvain"],
        help="Clustering method (default: leiden)"
    )
    parser.add_argument(
        "--clustering-resolutions",
        type=str,
        default="0.1,0.3,0.5,0.7,0.9",
        help="Comma-separated clustering resolutions (default: 0.1,0.3,0.5,0.7,0.9)"
    )
    parser.add_argument(
        "--default-resolution",
        type=float,
        default=0.5,
        help="Default resolution for visualization (default: 0.5)"
    )
    
    return parser.parse_args()


def normalize_and_hvg(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    n_top_genes: int = 2000
) -> sc.AnnData:
    """Normalize and select highly variable genes.
    
    Args:
        adata: Input AnnData object
        target_sum: Target sum for normalization
        n_top_genes: Number of HVGs
        
    Returns:
        Processed AnnData object
    """
    logger.info("=" * 70)
    logger.info("STEP 1: Normalization & Feature Selection")
    logger.info("=" * 70)
    
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize
    logger.info(f"Normalizing to target_sum={target_sum}...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    
    # Log transform
    logger.info("Log-transforming (log1p)...")
    sc.pp.log1p(adata)
    
    # Store normalized data
    adata.layers['normalized'] = adata.X.copy()
    
    # Highly variable genes
    logger.info(f"Selecting {n_top_genes} highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor='seurat',
        batch_key=None
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    logger.info(f"  ✓ Selected {n_hvg} highly variable genes")
    
    return adata


def scale_and_pca(
    adata: sc.AnnData,
    n_pcs: int = 50
) -> sc.AnnData:
    """Scale data and compute PCA.
    
    Args:
        adata: Input AnnData object
        n_pcs: Number of principal components
        
    Returns:
        Processed AnnData object
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 2: Scaling & PCA")
    logger.info("=" * 70)
    
    # Scale (on HVGs only)
    logger.info("Scaling data...")
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    logger.info(f"Computing PCA ({n_pcs} components)...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    
    # Report variance explained
    var_ratio = adata.uns['pca']['variance_ratio']
    cumsum_var = np.cumsum(var_ratio)
    
    logger.info(f"  PC1-10 variance explained: {cumsum_var[9]:.2%}")
    logger.info(f"  PC1-30 variance explained: {cumsum_var[29]:.2%}")
    logger.info(f"  PC1-50 variance explained: {cumsum_var[min(49, len(cumsum_var)-1)]:.2%}")
    
    return adata


def neighbors_and_umap(
    adata: sc.AnnData,
    n_neighbors: int = 15,
    n_pcs: int = 50,
    umap_min_dist: float = 0.5
) -> sc.AnnData:
    """Compute neighborhood graph and UMAP.
    
    Args:
        adata: Input AnnData object
        n_neighbors: Number of neighbors
        n_pcs: Number of PCs to use
        umap_min_dist: UMAP min_dist
        
    Returns:
        Processed AnnData object
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 3: Neighbors & UMAP")
    logger.info("=" * 70)
    
    # Compute neighbors
    logger.info(f"Computing neighborhood graph (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Compute UMAP
    logger.info(f"Computing UMAP (min_dist={umap_min_dist})...")
    sc.tl.umap(adata, min_dist=umap_min_dist)
    logger.info("  ✓ UMAP coordinates stored in .obsm['X_umap']")
    
    return adata


def multi_resolution_clustering(
    adata: sc.AnnData,
    method: str = 'leiden',
    resolutions: List[float] = [0.1, 0.3, 0.5, 0.7, 0.9],
    default_resolution: float = 0.5
) -> sc.AnnData:
    """Run clustering at multiple resolutions (SLC).
    
    Args:
        adata: Input AnnData object
        method: Clustering method ('leiden' or 'louvain')
        resolutions: List of resolutions to test
        default_resolution: Default resolution for visualization
        
    Returns:
        Processed AnnData object with clustering at all resolutions
    """
    logger.info("\n" + "=" * 70)
    logger.info(f"STEP 4: Multi-Resolution {method.capitalize()} Clustering (SLC)")
    logger.info("=" * 70)
    
    clustering_results = {}
    
    for res in resolutions:
        logger.info(f"Running {method} at resolution={res}...")
        
        if method == 'leiden':
            sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
            n_clusters = len(adata.obs[f'leiden_{res}'].unique())
        else:  # louvain
            sc.tl.louvain(adata, resolution=res, key_added=f'louvain_{res}')
            n_clusters = len(adata.obs[f'louvain_{res}'].unique())
        
        # Convert float key to string for H5AD compatibility
        clustering_results[str(res)] = n_clusters
        logger.info(f"  ✓ Resolution {res}: {n_clusters} clusters")
    
    # Set default clustering column
    default_key = f'{method}_{default_resolution}'
    if default_key in adata.obs.columns:
        adata.obs[method] = adata.obs[default_key].copy()
        logger.info(f"\n✓ Default clustering set to {default_key} ({clustering_results[str(default_resolution)]} clusters)")
    
    # Save clustering summary
    adata.uns['clustering_summary'] = {
        'method': method,
        'resolutions': resolutions,
        'default_resolution': default_resolution,
        'n_clusters_per_resolution': clustering_results
    }
    
    return adata


def generate_plots(
    adata: sc.AnnData,
    output_dir: Path,
    method: str = 'leiden',
    resolutions: List[float] = [0.1, 0.3, 0.5, 0.7, 0.9]
):
    """Generate standard processing plots.
    
    Args:
        adata: Processed AnnData object
        output_dir: Output directory
        method: Clustering method
        resolutions: List of resolutions
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 5: Generating Plots")
    logger.info("=" * 70)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. PCA variance plot
    logger.info("  - PCA variance explained...")
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    plt.savefig(output_dir / "pca_variance.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Multi-resolution clustering UMAPs
    logger.info("  - Multi-resolution clustering UMAPs...")
    n_res = len(resolutions)
    fig, axes = plt.subplots(1, n_res, figsize=(5 * n_res, 4))
    
    if n_res == 1:
        axes = [axes]
    
    for idx, res in enumerate(resolutions):
        key = f'{method}_{res}'
        if key in adata.obs.columns:
            sc.pl.umap(
                adata,
                color=key,
                ax=axes[idx],
                title=f'Resolution {res}',
                legend_loc='on data',
                legend_fontsize=6,
                show=False
            )
    
    plt.tight_layout()
    plt.savefig(output_dir / f"clustering_multi_resolution.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    ✓ Saved: clustering_multi_resolution.png")
    
    # 3. Individual resolution UMAPs
    for res in resolutions:
        key = f'{method}_{res}'
        if key in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.umap(
                adata,
                color=key,
                ax=ax,
                title=f'{method.capitalize()} Clustering (res={res})',
                show=False
            )
            plt.savefig(output_dir / f"umap_{method}_res{res}.png", dpi=300, bbox_inches='tight')
            plt.close()
    
    logger.info(f"  ✓ All plots saved to: {output_dir}/")


def export_metadata(
    adata: sc.AnnData,
    output_dir: Path,
    method: str = 'leiden',
    resolutions: List[float] = [0.1, 0.3, 0.5, 0.7, 0.9]
):
    """Export UMAP coordinates and metadata for dashboard.
    
    Args:
        adata: Processed AnnData object
        output_dir: Output directory
        method: Clustering method
        resolutions: List of resolutions
    """
    logger.info("\n" + "=" * 70)
    logger.info("STEP 6: Exporting Metadata for Dashboard")
    logger.info("=" * 70)
    
    # UMAP coordinates
    umap_df = pd.DataFrame(
        adata.obsm['X_umap'],
        columns=['UMAP_1', 'UMAP_2'],
        index=adata.obs.index
    )
    
    # Add clustering results
    for res in resolutions:
        key = f'{method}_{res}'
        if key in adata.obs.columns:
            umap_df[key] = adata.obs[key].values
    
    # Add QC metrics if available
    qc_cols = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    for col in qc_cols:
        if col in adata.obs.columns:
            umap_df[col] = adata.obs[col].values
    
    umap_path = output_dir / "umap_coordinates.csv"
    umap_df.to_csv(umap_path)
    logger.info(f"  ✓ Saved UMAP coordinates: umap_coordinates.csv")
    
    # Cell metadata
    metadata_cols = list(adata.obs.columns)
    metadata_df = adata.obs[metadata_cols].copy()
    metadata_path = output_dir / "cell_metadata.csv"
    metadata_df.to_csv(metadata_path)
    logger.info(f"  ✓ Saved cell metadata: cell_metadata.csv")


def main():
    """Main execution function."""
    args = parse_args()
    
    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("scAnnex Standard Processing Pipeline (SLC)")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Output Directory: {output_dir}")
    
    try:
        # Load data
        logger.info("\nLoading data...")
        adata = sc.read_h5ad(args.input)
        logger.info(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Parse resolutions
        resolutions = [float(r) for r in args.clustering_resolutions.split(',')]
        logger.info(f"\nClustering resolutions: {resolutions}")
        logger.info(f"Default resolution: {args.default_resolution}")
        
        # Run standard processing pipeline
        adata = normalize_and_hvg(adata, args.target_sum, args.n_top_genes)
        adata = scale_and_pca(adata, args.n_pcs)
        adata = neighbors_and_umap(adata, args.n_neighbors, args.n_pcs, args.umap_min_dist)
        adata = multi_resolution_clustering(
            adata,
            args.clustering_method,
            resolutions,
            args.default_resolution
        )
        
        # Generate plots
        generate_plots(adata, output_dir, args.clustering_method, resolutions)
        
        # Export metadata for dashboard
        export_metadata(adata, output_dir, args.clustering_method, resolutions)
        
        # Save processed data
        logger.info(f"\nSaving processed data to: {args.output}")
        adata.write_h5ad(args.output, compression='gzip')
        
        logger.info("\n" + "=" * 70)
        logger.info("✓ Standard Processing Complete")
        logger.info("=" * 70)
        logger.info(f"\nOutputs:")
        logger.info(f"  Processed data: {args.output}")
        logger.info(f"  Plots: {output_dir}/")
        logger.info(f"  Metadata: {output_dir}/umap_coordinates.csv")
        logger.info(f"  Metadata: {output_dir}/cell_metadata.csv")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    # Import pandas here for metadata export
    import pandas as pd
    main()
