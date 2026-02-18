#!/usr/bin/env python3
"""
Enhanced Normalization and Batch Integration

This module implements normalization and Harmony-based batch integration
following scAnnex specifications (InitProject.md Phase 4).

Features:
- Log normalization with configurable target sum
- Harmony integration using scanpy.external.pp.harmony_integrate
- Integration quality metrics (ASW, batch mixing entropy)
- Before/after UMAP visualization
- Pre-integration checkpoints for rollback

Author: scAnnex Development Team
"""

import argparse
import json
import logging
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

# Suppress warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def parse_args():
    """Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Enhanced normalization and batch integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Normalization only (no batches)
  %(prog)s --input qc_filtered.h5ad --output normalized.h5ad

  # With Harmony integration
  %(prog)s --input qc_filtered.h5ad --output integrated.h5ad \\
    --batch-key batch --run-integration

  # Save pre-integration checkpoint
  %(prog)s --input qc_filtered.h5ad --output integrated.h5ad \\
    --batch-key batch --run-integration \\
    --save-checkpoint pre_integration.h5ad
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input H5AD file (QC-filtered)"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output H5AD file (normalized/integrated)"
    )
    
    # Normalization parameters
    parser.add_argument(
        "--normalization-method",
        type=str,
        default="log",
        choices=["log", "scran"],
        help="Normalization method (default: log)"
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Target sum for normalization (default: 10000)"
    )
    
    # Feature selection
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Number of highly variable genes (default: 2000)"
    )
    
    # PCA parameters
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of principal components (default: 50)"
    )
    
    # Integration parameters
    parser.add_argument(
        "--run-integration",
        action="store_true",
        help="Run batch integration (requires --batch-key)"
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        default="batch",
        help="Column in .obs for batch labels (default: batch)"
    )
    parser.add_argument(
        "--integration-method",
        type=str,
        default="harmony",
        choices=["harmony"],
        help="Integration method (default: harmony)"
    )
    parser.add_argument(
        "--harmony-theta",
        type=float,
        default=2.0,
        help="Harmony diversity clustering penalty (default: 2.0)"
    )
    parser.add_argument(
        "--harmony-max-iter",
        type=int,
        default=10,
        help="Harmony maximum iterations (default: 10)"
    )
    
    # UMAP parameters
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for UMAP (default: 15)"
    )
    parser.add_argument(
        "--umap-min-dist",
        type=float,
        default=0.5,
        help="UMAP minimum distance (default: 0.5)"
    )
    
    # Output options
    parser.add_argument(
        "--output-dir",
        type=str,
        default="integration_results",
        help="Output directory for plots and metrics"
    )
    parser.add_argument(
        "--save-checkpoint",
        type=str,
        default=None,
        help="Save pre-integration checkpoint H5AD"
    )
    
    return parser.parse_args()


def normalize_data(
    adata: sc.AnnData,
    method: str = "log",
    target_sum: float = 1e4
) -> sc.AnnData:
    """Normalize count data.
    
    Args:
        adata: Input AnnData object
        method: Normalization method (log or scran)
        target_sum: Target sum for normalization
        
    Returns:
        Normalized AnnData object
    """
    logger.info(f"Normalizing data (method: {method}, target_sum: {target_sum})...")
    
    # Ensure raw counts are preserved in layers
    if 'counts' not in adata.layers:
        logger.info("  Storing raw counts in .layers['counts']")
        adata.layers['counts'] = adata.X.copy()
    
    if method == "log":
        # Standard log normalization
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        logger.info(f"  ✓ Log normalization complete")
        
    elif method == "scran":
        # Scran normalization
        try:
            import scanpy.external as sce
            sce.pp.scran_normalization(adata)
            logger.info(f"  ✓ Scran normalization complete")
        except ImportError:
            logger.warning("  ⚠ scran not available, falling back to log normalization")
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
    
    return adata


def select_highly_variable_genes(
    adata: sc.AnnData,
    n_top_genes: int = 2000
) -> sc.AnnData:
    """Select highly variable genes.
    
    Args:
        adata: Normalized AnnData object
        n_top_genes: Number of HVGs to select
        
    Returns:
        AnnData with HVG annotations in .var
    """
    logger.info(f"Selecting {n_top_genes} highly variable genes...")
    
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor='seurat',
        subset=False  # Keep all genes, mark HVGs
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    logger.info(f"  ✓ Identified {n_hvg} highly variable genes")
    
    return adata


def scale_and_pca(
    adata: sc.AnnData,
    n_pcs: int = 50,
    max_value: float = 10.0
) -> sc.AnnData:
    """Scale data and compute PCA.
    
    Args:
        adata: AnnData with HVGs selected
        n_pcs: Number of principal components
        max_value: Maximum value for scaling (clip outliers)
        
    Returns:
        AnnData with PCA in .obsm['X_pca']
    """
    logger.info(f"Scaling data and computing PCA ({n_pcs} components)...")
    
    # Scale using HVGs only
    sc.pp.scale(adata, max_value=max_value)
    
    # Compute PCA
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    
    # Log variance explained
    var_ratio = adata.uns['pca']['variance_ratio']
    cumsum = np.cumsum(var_ratio)
    logger.info(f"  ✓ PCA complete")
    logger.info(f"    PC1 variance: {var_ratio[0]:.3f}")
    logger.info(f"    First 10 PCs: {cumsum[9]:.3f}")
    logger.info(f"    First 30 PCs: {cumsum[29]:.3f}")
    
    return adata


def validate_batch_variation(
    adata: sc.AnnData,
    batch_key: str
) -> bool:
    """Validate that batch variation exists before integration.
    
    Per InitProject.md: Harmony can fail silently without variation.
    
    Args:
        adata: AnnData object
        batch_key: Batch column in .obs
        
    Returns:
        True if sufficient batch variation, False otherwise
    """
    logger.info(f"Validating batch variation in '{batch_key}'...")
    
    if batch_key not in adata.obs.columns:
        logger.error(f"  ✗ Batch key '{batch_key}' not found in .obs")
        return False
    
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)
    
    logger.info(f"  Found {n_batches} batches: {list(batches)}")
    
    if n_batches < 2:
        logger.warning(f"  ⚠ Only {n_batches} batch(es) - integration not needed")
        return False
    
    # Check cell counts per batch
    batch_counts = adata.obs[batch_key].value_counts()
    logger.info(f"  Cells per batch:")
    for batch, count in batch_counts.items():
        logger.info(f"    {batch}: {count} cells")
    
    # Warn if batch imbalance
    min_cells = batch_counts.min()
    max_cells = batch_counts.max()
    if max_cells / min_cells > 10:
        logger.warning(f"  ⚠ Large batch imbalance (ratio: {max_cells/min_cells:.1f})")
    
    return True


def run_harmony_integration(
    adata: sc.AnnData,
    batch_key: str,
    theta: float = 2.0,
    max_iter: int = 10
) -> sc.AnnData:
    """Run Harmony batch integration.
    
    Uses scanpy.external.pp.harmony_integrate as recommended in InitProject.md.
    
    Args:
        adata: AnnData with PCA computed
        batch_key: Batch column in .obs
        theta: Diversity clustering penalty (higher = more aggressive)
        max_iter: Maximum iterations
        
    Returns:
        AnnData with Harmony coordinates in .obsm['X_pca_harmony']
    """
    logger.info(f"Running Harmony integration...")
    logger.info(f"  Parameters: theta={theta}, max_iter={max_iter}")
    
    try:
        import scanpy.external as sce
        
        # Run Harmony
        sce.pp.harmony_integrate(
            adata,
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_pca_harmony',
            theta=theta,
            max_iter_harmony=max_iter
        )
        
        logger.info(f"  ✓ Harmony integration complete")
        logger.info(f"    Integrated PCA stored in .obsm['X_pca_harmony']")
        
        return adata
        
    except ImportError:
        logger.error("  ✗ harmonypy not available")
        logger.error("    Install with: pip install harmonypy")
        raise
    except Exception as e:
        logger.error(f"  ✗ Harmony integration failed: {e}")
        raise


def compute_umap(
    adata: sc.AnnData,
    n_neighbors: int = 15,
    min_dist: float = 0.5,
    use_rep: str = 'X_pca'
) -> sc.AnnData:
    """Compute UMAP embedding.
    
    Args:
        adata: AnnData object
        n_neighbors: Number of neighbors
        min_dist: Minimum distance
        use_rep: Representation to use from .obsm
        
    Returns:
        AnnData with UMAP in .obsm['X_umap']
    """
    logger.info(f"Computing UMAP (n_neighbors={n_neighbors}, min_dist={min_dist}, use_rep={use_rep})...")
    
    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)
    
    # Compute UMAP
    sc.tl.umap(adata, min_dist=min_dist)
    
    logger.info(f"  ✓ UMAP complete, stored in .obsm['X_umap']")
    
    return adata


def calculate_integration_metrics(
    adata: sc.AnnData,
    batch_key: str,
    use_rep: str = 'X_pca'
) -> Dict:
    """Calculate integration quality metrics.
    
    Implements ASW (Average Silhouette Width) and batch mixing entropy.
    
    Args:
        adata: Integrated AnnData
        batch_key: Batch column
        use_rep: Representation to use
        
    Returns:
        Dict of metric: value
    """
    logger.info(f"Calculating integration metrics...")
    
    metrics = {}
    
    try:
        from sklearn.metrics import silhouette_score
        
        # Calculate ASW for batch (lower is better for integration)
        if use_rep in adata.obsm:
            X = adata.obsm[use_rep]
            batch_labels = adata.obs[batch_key].values
            
            asw_batch = silhouette_score(X, batch_labels, metric='euclidean')
            metrics['asw_batch'] = float(asw_batch)
            logger.info(f"  ASW (batch): {asw_batch:.3f} (lower = better mixing)")
        
        # Calculate simple batch mixing entropy
        # For each cell's k nearest neighbors, measure batch diversity
        if 'neighbors' in adata.uns:
            from scipy.sparse import csr_matrix
            from scipy.stats import entropy
            
            # Get neighbor graph
            connectivities = adata.obsp['connectivities']
            batch_labels = adata.obs[batch_key].astype('category')
            batch_codes = batch_labels.cat.codes.values
            n_batches = len(batch_labels.cat.categories)
            
            # Calculate batch mixing for each cell
            mixing_scores = []
            for i in range(adata.n_obs):
                # Get neighbors
                neighbors_idx = connectivities[i].nonzero()[1]
                if len(neighbors_idx) == 0:
                    continue
                
                # Count batches in neighborhood
                neighbor_batches = batch_codes[neighbors_idx]
                batch_counts = np.bincount(neighbor_batches, minlength=n_batches)
                batch_probs = batch_counts / batch_counts.sum()
                
                # Calculate entropy (higher = better mixing)
                mix_entropy = entropy(batch_probs)
                mixing_scores.append(mix_entropy)
            
            mean_entropy = np.mean(mixing_scores)
            max_entropy = np.log(n_batches)  # Maximum possible entropy
            normalized_entropy = mean_entropy / max_entropy if max_entropy > 0 else 0
            
            metrics['batch_mixing_entropy'] = float(mean_entropy)
            metrics['batch_mixing_entropy_normalized'] = float(normalized_entropy)
            logger.info(f"  Batch mixing entropy: {mean_entropy:.3f} / {max_entropy:.3f}")
            logger.info(f"  Normalized: {normalized_entropy:.3f} (higher = better, max = 1.0)")
        
    except ImportError as e:
        logger.warning(f"  ⚠ Could not calculate all metrics: {e}")
    except Exception as e:
        logger.warning(f"  ⚠ Metric calculation failed: {e}")
    
    return metrics


def plot_before_after_umap(
    adata: sc.AnnData,
    adata_before: sc.AnnData,
    batch_key: str,
    output_dir: Path
) -> None:
    """Generate before/after UMAP comparison plots.
    
    Args:
        adata: Integrated AnnData (with X_pca_harmony)
        adata_before: Pre-integration AnnData (with X_pca)
        batch_key: Batch column
        output_dir: Output directory
    """
    logger.info("Generating before/after UMAP comparison...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Before integration - colored by batch
    sc.pl.umap(adata_before, color=batch_key, ax=axes[0, 0], show=False, title='Before Integration (Batch)')
    
    # After integration - colored by batch
    sc.pl.umap(adata, color=batch_key, ax=axes[0, 1], show=False, title='After Integration (Batch)')
    
    # Before integration - colored by sample_id if available
    if 'sample_id' in adata_before.obs.columns:
        sc.pl.umap(adata_before, color='sample_id', ax=axes[1, 0], show=False, title='Before Integration (Sample)')
    else:
        axes[1, 0].axis('off')
    
    # After integration - colored by sample_id
    if 'sample_id' in adata.obs.columns:
        sc.pl.umap(adata, color='sample_id', ax=axes[1, 1], show=False, title='After Integration (Sample)')
    else:
        axes[1, 1].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / "integration_before_after_umap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved: integration_before_after_umap.png")
    
    # Also create side-by-side batch comparison
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    sc.pl.umap(adata_before, color=batch_key, ax=axes[0], show=False, 
               title='PCA (Before Integration)', legend_loc='on data', legend_fontsize=10)
    sc.pl.umap(adata, color=batch_key, ax=axes[1], show=False,
               title='Harmony (After Integration)', legend_loc='on data', legend_fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / "integration_batch_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved: integration_batch_comparison.png")


def save_integration_report(
    metrics: Dict,
    output_dir: Path,
    adata_before: sc.AnnData,
    adata_after: sc.AnnData,
    batch_key: str
) -> None:
    """Save integration quality report.
    
    Args:
        metrics: Integration metrics
        output_dir: Output directory
        adata_before: Pre-integration AnnData
        adata_after: Post-integration AnnData
        batch_key: Batch column
    """
    logger.info("Saving integration report...")
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'batch_key': batch_key,
        'n_cells': int(adata_after.n_obs),
        'n_genes': int(adata_after.n_vars),
        'n_batches': int(adata_after.obs[batch_key].nunique()),
        'batch_sizes': adata_after.obs[batch_key].value_counts().to_dict(),
        'integration_metrics': metrics,
        'pca_variance': {
            'before_integration': {
                'pc1': float(adata_before.uns['pca']['variance_ratio'][0]),
                'cumsum_10': float(np.cumsum(adata_before.uns['pca']['variance_ratio'])[9]),
                'cumsum_30': float(np.cumsum(adata_before.uns['pca']['variance_ratio'])[29])
            }
        }
    }
    
    # Save JSON report
    report_path = output_dir / "integration_report.json"
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"  ✓ Saved: integration_report.json")


def main():
    """Main execution function."""
    args = parse_args()
    
    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("scAnnex Enhanced Normalization & Integration")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Results directory: {output_dir}")
    
    try:
        # Load data
        logger.info("\nLoading data...")
        adata = sc.read_h5ad(args.input)
        logger.info(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Step 1: Normalization
        logger.info("\n" + "─" * 70)
        logger.info("STEP 1: Normalization")
        logger.info("─" * 70)
        adata = normalize_data(adata, args.normalization_method, args.target_sum)
        
        # Step 2: Feature selection
        logger.info("\n" + "─" * 70)
        logger.info("STEP 2: Highly Variable Genes Selection")
        logger.info("─" * 70)
        adata = select_highly_variable_genes(adata, args.n_top_genes)
        
        # Step 3: Scaling and PCA
        logger.info("\n" + "─" * 70)
        logger.info("STEP 3: Scaling and PCA")
        logger.info("─" * 70)
        adata = scale_and_pca(adata, args.n_pcs)
        
        # Step 4: Batch integration (optional)
        if args.run_integration:
            logger.info("\n" + "─" * 70)
            logger.info("STEP 4: Batch Integration (Harmony)")
            logger.info("─" * 70)
            
            # Validate batch variation
            if not validate_batch_variation(adata, args.batch_key):
                logger.warning("⚠ Skipping integration due to insufficient batch variation")
                args.run_integration = False
            else:
                # Save checkpoint before integration
                if args.save_checkpoint:
                    logger.info(f"Saving pre-integration checkpoint: {args.save_checkpoint}")
                    adata.write_h5ad(args.save_checkpoint, compression='gzip')
                
                # Keep copy for before/after comparison
                adata_before = adata.copy()
                
                # Compute UMAP before integration
                logger.info("Computing UMAP before integration...")
                adata_before = compute_umap(adata_before, args.n_neighbors, args.umap_min_dist, use_rep='X_pca')
                
                # Run Harmony
                adata = run_harmony_integration(adata, args.batch_key, args.harmony_theta, args.harmony_max_iter)
                
                # Compute UMAP after integration (using Harmony PCs)
                logger.info("Computing UMAP after integration...")
                adata = compute_umap(adata, args.n_neighbors, args.umap_min_dist, use_rep='X_pca_harmony')
                
                # Calculate integration metrics
                logger.info("\n" + "─" * 70)
                logger.info("STEP 5: Integration Quality Metrics")
                logger.info("─" * 70)
                metrics_before = calculate_integration_metrics(adata_before, args.batch_key, use_rep='X_pca')
                metrics_after = calculate_integration_metrics(adata, args.batch_key, use_rep='X_pca_harmony')
                
                all_metrics = {
                    'before': metrics_before,
                    'after': metrics_after
                }
                
                logger.info("\nComparison:")
                if 'asw_batch' in metrics_before and 'asw_batch' in metrics_after:
                    improvement = metrics_before['asw_batch'] - metrics_after['asw_batch']
                    logger.info(f"  ASW improvement: {improvement:+.3f} (positive = better mixing)")
                
                # Generate before/after plots
                logger.info("\n" + "─" * 70)
                logger.info("STEP 6: Visualization")
                logger.info("─" * 70)
                plot_before_after_umap(adata, adata_before, args.batch_key, output_dir)
                
                # Save integration report
                save_integration_report(all_metrics, output_dir, adata_before, adata, args.batch_key)
        else:
            # No integration - just compute UMAP on PCA
            logger.info("\n" + "─" * 70)
            logger.info("STEP 4: UMAP (No Integration)")
            logger.info("─" * 70)
            adata = compute_umap(adata, args.n_neighbors, args.umap_min_dist, use_rep='X_pca')
        
        # Save final data
        logger.info(f"\nSaving processed data to: {args.output}")
        adata.write_h5ad(args.output, compression='gzip')
        
        logger.info("\n" + "=" * 70)
        logger.info("✓ Normalization and Integration Complete")
        logger.info("=" * 70)
        
        logger.info(f"\nFinal data:")
        logger.info(f"  Cells: {adata.n_obs}")
        logger.info(f"  Genes: {adata.n_vars}")
        logger.info(f"  .obsm keys: {list(adata.obsm.keys())}")
        logger.info(f"  .layers keys: {list(adata.layers.keys())}")
        
        logger.info(f"\nOutputs:")
        logger.info(f"  Processed data: {args.output}")
        if args.run_integration:
            logger.info(f"  Integration plots: {output_dir}/")
            logger.info(f"  Integration report: {output_dir}/integration_report.json")
        if args.save_checkpoint:
            logger.info(f"  Pre-integration checkpoint: {args.save_checkpoint}")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
