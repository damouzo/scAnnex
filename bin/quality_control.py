#!/usr/bin/env python3
"""
Enhanced Quality Control for scRNA-seq data

This module implements comprehensive QC metrics calculation and filtering
following scAnnex specifications (InitProject.md Phase 2).

Features:
- Multi-gene set QC: Mitochondrial, Ribosomal, Hemoglobin
- MAD-based automatic threshold detection
- Publication-quality visualizations
- Detailed QC reporting

Author: scAnnex Development Team
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

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
plt.rcParams['font.size'] = 10


def parse_args():
    """Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Enhanced quality control for scRNA-seq data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Automatic MAD-based thresholds
  %(prog)s --input data.h5ad --output qc_filtered.h5ad \\
    --qc-dir qc_results/ --use-mad-thresholds

  # Manual thresholds
  %(prog)s --input data.h5ad --output qc_filtered.h5ad \\
    --min-genes 200 --max-mito 20 --max-ribo 40

  # Skip filtering (metrics only)
  %(prog)s --input data.h5ad --output data.h5ad \\
    --qc-dir qc_results/ --skip-filtering
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input H5AD file"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output H5AD file (filtered)"
    )
    parser.add_argument(
        "--qc-dir",
        type=str,
        default="qc_results",
        help="Output directory for QC plots and metrics"
    )
    
    # Filtering thresholds
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum genes per cell (default: 200)"
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=3,
        help="Minimum cells per gene (default: 3)"
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=None,
        help="Maximum genes per cell (default: auto via MAD)"
    )
    parser.add_argument(
        "--max-counts",
        type=int,
        default=None,
        help="Maximum counts per cell (default: auto via MAD)"
    )
    parser.add_argument(
        "--max-mito",
        type=float,
        default=None,
        help="Maximum mitochondrial %% (default: auto via MAD)"
    )
    parser.add_argument(
        "--max-ribo",
        type=float,
        default=None,
        help="Maximum ribosomal %% (default: no filter)"
    )
    parser.add_argument(
        "--max-hb",
        type=float,
        default=None,
        help="Maximum hemoglobin %% (default: no filter)"
    )
    
    # MAD-based thresholds
    parser.add_argument(
        "--use-mad-thresholds",
        action="store_true",
        help="Use MAD-based automatic thresholds"
    )
    parser.add_argument(
        "--mad-threshold",
        type=float,
        default=5.0,
        help="MAD multiplier for outlier detection (default: 5.0)"
    )
    
    # Options
    parser.add_argument(
        "--skip-filtering",
        action="store_true",
        help="Skip filtering, only calculate metrics and plots"
    )
    
    return parser.parse_args()


def calculate_mad_thresholds(
    adata: sc.AnnData,
    nmads: float = 5.0
) -> Dict[str, Tuple[float, float]]:
    """Calculate MAD-based thresholds for QC metrics.
    
    Uses median absolute deviation (MAD) to identify outliers:
    threshold = median ± nmads * MAD
    
    Args:
        adata: AnnData object with calculated QC metrics
        nmads: Number of MADs for threshold (default: 5.0)
        
    Returns:
        Dict of metric: (lower_threshold, upper_threshold)
    """
    logger.info(f"Calculating MAD-based thresholds (nmads={nmads})...")
    
    thresholds = {}
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    
    # Add ribosomal and hemoglobin if available
    if 'pct_counts_ribo' in adata.obs.columns:
        metrics.append('pct_counts_ribo')
    if 'pct_counts_hb' in adata.obs.columns:
        metrics.append('pct_counts_hb')
    
    for metric in metrics:
        values = adata.obs[metric]
        median = np.median(values)
        mad = np.median(np.abs(values - median))
        
        # For percentage metrics, use upper bound only
        # For counts, use both bounds
        if metric.startswith('pct_'):
            lower = 0  # Don't filter low percentages
            upper = median + nmads * mad
        else:
            lower = max(0, median - nmads * mad)
            upper = median + nmads * mad
        
        thresholds[metric] = (lower, upper)
        
        logger.info(f"  {metric}:")
        logger.info(f"    Median: {median:.2f}, MAD: {mad:.2f}")
        logger.info(f"    Threshold: [{lower:.2f}, {upper:.2f}]")
    
    return thresholds


def calculate_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """Calculate comprehensive QC metrics.
    
    Calculates metrics for:
    - Mitochondrial genes (MT-)
    - Ribosomal genes (RPL, RPS)
    - Hemoglobin genes (HB)
    
    Args:
        adata: Input AnnData object
        
    Returns:
        AnnData object with QC metrics in .obs
    """
    logger.info("Calculating comprehensive QC metrics...")
    
    # Identify gene sets
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
    
    # Log gene set sizes
    n_mt = adata.var['mt'].sum()
    n_ribo = adata.var['ribo'].sum()
    n_hb = adata.var['hb'].sum()
    
    logger.info(f"  Identified gene sets:")
    logger.info(f"    Mitochondrial (MT-): {n_mt} genes")
    logger.info(f"    Ribosomal (RPL/RPS): {n_ribo} genes")
    logger.info(f"    Hemoglobin (HB): {n_hb} genes")
    
    if n_mt == 0:
        logger.warning("  ⚠ No mitochondrial genes found - check gene naming")
    if n_ribo == 0:
        logger.warning("  ⚠ No ribosomal genes found")
    if n_hb == 0:
        logger.info("  ℹ No hemoglobin genes found (normal for non-blood samples)")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # Rename columns for clarity
    adata.obs.rename(columns={
        'pct_counts_mt': 'pct_counts_mt',
        'pct_counts_ribo': 'pct_counts_ribo',
        'pct_counts_hb': 'pct_counts_hb'
    }, inplace=True, errors='ignore')
    
    # Log summary statistics
    logger.info("  QC Metrics Summary:")
    logger.info(f"    n_genes_by_counts: {adata.obs['n_genes_by_counts'].median():.0f} (median)")
    logger.info(f"    total_counts: {adata.obs['total_counts'].median():.0f} (median)")
    logger.info(f"    pct_counts_mt: {adata.obs['pct_counts_mt'].median():.2f}% (median)")
    if 'pct_counts_ribo' in adata.obs.columns:
        logger.info(f"    pct_counts_ribo: {adata.obs['pct_counts_ribo'].median():.2f}% (median)")
    if 'pct_counts_hb' in adata.obs.columns:
        logger.info(f"    pct_counts_hb: {adata.obs['pct_counts_hb'].median():.2f}% (median)")
    
    return adata


def plot_qc_metrics(
    adata: sc.AnnData,
    output_dir: Path,
    prefix: str = "qc",
    thresholds: Optional[Dict] = None
) -> None:
    """Generate publication-quality QC plots.
    
    Creates:
    - Violin plots for all QC metrics
    - Scatter plots (counts vs genes, colored by MT%)
    - Distribution plots with threshold lines
    
    Args:
        adata: AnnData object with QC metrics
        output_dir: Directory for output plots
        prefix: Filename prefix
        thresholds: Optional dict of metric thresholds to plot
    """
    logger.info(f"Generating QC plots ({prefix})...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine which metrics to plot
    base_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    extra_metrics = []
    
    if 'pct_counts_ribo' in adata.obs.columns:
        extra_metrics.append('pct_counts_ribo')
    if 'pct_counts_hb' in adata.obs.columns:
        extra_metrics.append('pct_counts_hb')
    
    all_metrics = base_metrics + extra_metrics
    n_metrics = len(all_metrics)
    
    # 1. Violin plots
    fig, axes = plt.subplots(1, n_metrics, figsize=(5 * n_metrics, 4))
    if n_metrics == 1:
        axes = [axes]
    
    for ax, metric in zip(axes, all_metrics):
        sns.violinplot(data=adata.obs, y=metric, ax=ax, color='lightblue')
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_xlabel('')
        
        # Add threshold lines if provided
        if thresholds and metric in thresholds:
            lower, upper = thresholds[metric]
            if upper < np.inf:
                ax.axhline(upper, color='red', linestyle='--', linewidth=1.5, 
                          label=f'Upper: {upper:.1f}')
            if lower > 0 and not metric.startswith('pct_'):
                ax.axhline(lower, color='red', linestyle='--', linewidth=1.5,
                          label=f'Lower: {lower:.1f}')
            ax.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_dir / f"{prefix}_violin.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  ✓ Saved: {prefix}_violin.png")
    
    # 2. Scatter plot: total_counts vs n_genes_by_counts
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(
        adata.obs['total_counts'],
        adata.obs['n_genes_by_counts'],
        c=adata.obs['pct_counts_mt'],
        cmap='viridis',
        alpha=0.6,
        s=5
    )
    ax.set_xlabel('Total Counts', fontsize=12)
    ax.set_ylabel('Number of Genes', fontsize=12)
    ax.set_title('QC Metrics Scatter Plot', fontsize=14)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Mitochondrial %', fontsize=10)
    
    # Add threshold lines
    if thresholds:
        if 'total_counts' in thresholds:
            _, upper_counts = thresholds['total_counts']
            if upper_counts < np.inf:
                ax.axvline(upper_counts, color='red', linestyle='--', linewidth=1.5,
                          alpha=0.7, label='Count threshold')
        if 'n_genes_by_counts' in thresholds:
            _, upper_genes = thresholds['n_genes_by_counts']
            if upper_genes < np.inf:
                ax.axhline(upper_genes, color='red', linestyle='--', linewidth=1.5,
                          alpha=0.7, label='Gene threshold')
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_dir / f"{prefix}_scatter.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  ✓ Saved: {prefix}_scatter.png")
    
    # 3. Distribution plots with thresholds
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    plot_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    if 'pct_counts_ribo' in adata.obs.columns:
        plot_metrics.append('pct_counts_ribo')
    
    for idx, metric in enumerate(plot_metrics[:4]):
        if idx >= len(axes):
            break
        
        ax = axes[idx]
        values = adata.obs[metric]
        
        # Histogram
        ax.hist(values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_xlabel(metric.replace('_', ' ').title(), fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f'Distribution: {metric}', fontsize=12)
        
        # Add median line
        median = np.median(values)
        ax.axvline(median, color='green', linestyle='-', linewidth=2,
                  label=f'Median: {median:.1f}')
        
        # Add threshold lines
        if thresholds and metric in thresholds:
            lower, upper = thresholds[metric]
            if upper < np.inf:
                ax.axvline(upper, color='red', linestyle='--', linewidth=2,
                          label=f'Upper: {upper:.1f}')
            if lower > 0 and not metric.startswith('pct_'):
                ax.axvline(lower, color='red', linestyle='--', linewidth=2,
                          label=f'Lower: {lower:.1f}')
        
        ax.legend(fontsize=9)
    
    # Hide unused subplots
    for idx in range(len(plot_metrics), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / f"{prefix}_distributions.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  ✓ Saved: {prefix}_distributions.png")


def apply_qc_filters(
    adata: sc.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_genes: Optional[int] = None,
    max_counts: Optional[int] = None,
    max_mito: Optional[float] = None,
    max_ribo: Optional[float] = None,
    max_hb: Optional[float] = None,
    thresholds: Optional[Dict] = None
) -> Tuple[sc.AnnData, Dict]:
    """Apply QC filters to remove low-quality cells and genes.
    
    Args:
        adata: Input AnnData object
        min_genes: Minimum genes per cell
        min_cells: Minimum cells per gene
        max_genes: Maximum genes per cell (or use thresholds)
        max_counts: Maximum counts per cell (or use thresholds)
        max_mito: Maximum mitochondrial % (or use thresholds)
        max_ribo: Maximum ribosomal %
        max_hb: Maximum hemoglobin %
        thresholds: MAD-based thresholds (overrides manual values)
        
    Returns:
        Tuple of (filtered_adata, filtering_stats)
    """
    logger.info("Applying QC filters...")
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Use MAD thresholds if provided
    if thresholds:
        logger.info("  Using MAD-based thresholds")
        if 'n_genes_by_counts' in thresholds and max_genes is None:
            _, max_genes = thresholds['n_genes_by_counts']
            max_genes = int(max_genes) if max_genes < np.inf else None
        if 'total_counts' in thresholds and max_counts is None:
            _, max_counts = thresholds['total_counts']
            max_counts = int(max_counts) if max_counts < np.inf else None
        if 'pct_counts_mt' in thresholds and max_mito is None:
            _, max_mito = thresholds['pct_counts_mt']
    
    # Track filtering steps
    filtering_stats = {
        'cells_before': n_cells_before,
        'genes_before': n_genes_before
    }
    
    # Filter cells - minimum genes
    cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    cells_removed = cells_before - adata.n_obs
    logger.info(f"  Removed {cells_removed} cells with < {min_genes} genes")
    filtering_stats['cells_min_genes'] = cells_removed
    
    # Filter cells - maximum genes
    if max_genes:
        cells_before = adata.n_obs
        adata = adata[adata.obs.n_genes_by_counts < max_genes, :].copy()
        cells_removed = cells_before - adata.n_obs
        logger.info(f"  Removed {cells_removed} cells with > {max_genes} genes")
        filtering_stats['cells_max_genes'] = cells_removed
    
    # Filter cells - maximum counts
    if max_counts:
        cells_before = adata.n_obs
        adata = adata[adata.obs.total_counts < max_counts, :].copy()
        cells_removed = cells_before - adata.n_obs
        logger.info(f"  Removed {cells_removed} cells with > {max_counts} counts")
        filtering_stats['cells_max_counts'] = cells_removed
    
    # Filter cells - mitochondrial %
    if max_mito:
        cells_before = adata.n_obs
        adata = adata[adata.obs.pct_counts_mt < max_mito, :].copy()
        cells_removed = cells_before - adata.n_obs
        logger.info(f"  Removed {cells_removed} cells with > {max_mito:.1f}% MT")
        filtering_stats['cells_max_mito'] = cells_removed
    
    # Filter cells - ribosomal %
    if max_ribo and 'pct_counts_ribo' in adata.obs.columns:
        cells_before = adata.n_obs
        adata = adata[adata.obs.pct_counts_ribo < max_ribo, :].copy()
        cells_removed = cells_before - adata.n_obs
        logger.info(f"  Removed {cells_removed} cells with > {max_ribo:.1f}% ribosomal")
        filtering_stats['cells_max_ribo'] = cells_removed
    
    # Filter cells - hemoglobin %
    if max_hb and 'pct_counts_hb' in adata.obs.columns:
        cells_before = adata.n_obs
        adata = adata[adata.obs.pct_counts_hb < max_hb, :].copy()
        cells_removed = cells_before - adata.n_obs
        logger.info(f"  Removed {cells_removed} cells with > {max_hb:.1f}% hemoglobin")
        filtering_stats['cells_max_hb'] = cells_removed
    
    # Filter genes
    genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    genes_removed = genes_before - adata.n_vars
    logger.info(f"  Removed {genes_removed} genes expressed in < {min_cells} cells")
    filtering_stats['genes_min_cells'] = genes_removed
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    filtering_stats['cells_after'] = n_cells_after
    filtering_stats['genes_after'] = n_genes_after
    filtering_stats['cells_removed_total'] = n_cells_before - n_cells_after
    filtering_stats['genes_removed_total'] = n_genes_before - n_genes_after
    filtering_stats['cells_retained_pct'] = (n_cells_after / n_cells_before) * 100
    filtering_stats['genes_retained_pct'] = (n_genes_after / n_genes_before) * 100
    
    logger.info(f"  Final: {n_cells_after} cells ({filtering_stats['cells_retained_pct']:.1f}%) × "
                f"{n_genes_after} genes ({filtering_stats['genes_retained_pct']:.1f}%)")
    
    return adata, filtering_stats


def save_qc_report(
    filtering_stats: Dict,
    thresholds: Optional[Dict],
    output_dir: Path,
    adata_before: sc.AnnData,
    adata_after: Optional[sc.AnnData] = None
) -> None:
    """Save comprehensive QC report.
    
    Args:
        filtering_stats: Statistics from filtering
        thresholds: Applied thresholds
        output_dir: Output directory
        adata_before: AnnData before filtering
        adata_after: AnnData after filtering (if filtered)
    """
    logger.info("Generating QC report...")
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'filtering_statistics': filtering_stats,
        'thresholds_applied': thresholds if thresholds else 'manual',
        'qc_metrics_before': {
            'n_genes_by_counts': {
                'median': float(adata_before.obs['n_genes_by_counts'].median()),
                'mean': float(adata_before.obs['n_genes_by_counts'].mean()),
                'std': float(adata_before.obs['n_genes_by_counts'].std())
            },
            'total_counts': {
                'median': float(adata_before.obs['total_counts'].median()),
                'mean': float(adata_before.obs['total_counts'].mean()),
                'std': float(adata_before.obs['total_counts'].std())
            },
            'pct_counts_mt': {
                'median': float(adata_before.obs['pct_counts_mt'].median()),
                'mean': float(adata_before.obs['pct_counts_mt'].mean()),
                'std': float(adata_before.obs['pct_counts_mt'].std())
            }
        }
    }
    
    if adata_after is not None:
        report['qc_metrics_after'] = {
            'n_genes_by_counts': {
                'median': float(adata_after.obs['n_genes_by_counts'].median()),
                'mean': float(adata_after.obs['n_genes_by_counts'].mean()),
                'std': float(adata_after.obs['n_genes_by_counts'].std())
            },
            'total_counts': {
                'median': float(adata_after.obs['total_counts'].median()),
                'mean': float(adata_after.obs['total_counts'].mean()),
                'std': float(adata_after.obs['total_counts'].std())
            },
            'pct_counts_mt': {
                'median': float(adata_after.obs['pct_counts_mt'].median()),
                'mean': float(adata_after.obs['pct_counts_mt'].mean()),
                'std': float(adata_after.obs['pct_counts_mt'].std())
            }
        }
    
    # Save JSON report
    report_path = output_dir / "qc_report.json"
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    logger.info(f"  ✓ Saved: qc_report.json")
    
    # Save metrics CSV
    metrics_df = adata_before.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].describe()
    metrics_path = output_dir / "qc_metrics_summary.csv"
    metrics_df.to_csv(metrics_path)
    logger.info(f"  ✓ Saved: qc_metrics_summary.csv")


def main():
    """Main execution function."""
    args = parse_args()
    
    # Setup output directory
    output_dir = Path(args.qc_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("scAnnex Enhanced Quality Control")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"QC Directory: {output_dir}")
    
    try:
        # Load data
        logger.info("\nLoading data...")
        adata = sc.read_h5ad(args.input)
        logger.info(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Keep copy for before/after comparison
        adata_before = adata.copy()
        
        # Calculate QC metrics
        adata = calculate_qc_metrics(adata)
        
        # Calculate MAD thresholds if requested
        thresholds = None
        if args.use_mad_thresholds:
            thresholds = calculate_mad_thresholds(adata, nmads=args.mad_threshold)
        
        # Generate "before" plots
        logger.info("\nGenerating pre-filtering QC plots...")
        plot_qc_metrics(adata, output_dir, prefix="qc_before", thresholds=thresholds)
        
        # Apply filtering (unless skipped)
        if args.skip_filtering:
            logger.info("\n⚠ Filtering skipped (--skip-filtering enabled)")
            adata_filtered = adata
            filtering_stats = {'skipped': True}
        else:
            logger.info("\nApplying QC filters...")
            adata_filtered, filtering_stats = apply_qc_filters(
                adata,
                min_genes=args.min_genes,
                min_cells=args.min_cells,
                max_genes=args.max_genes,
                max_counts=args.max_counts,
                max_mito=args.max_mito,
                max_ribo=args.max_ribo,
                max_hb=args.max_hb,
                thresholds=thresholds
            )
            
            # Generate "after" plots
            logger.info("\nGenerating post-filtering QC plots...")
            plot_qc_metrics(adata_filtered, output_dir, prefix="qc_after", thresholds=None)
        
        # Save QC report
        save_qc_report(
            filtering_stats,
            thresholds,
            output_dir,
            adata_before,
            adata_filtered if not args.skip_filtering else None
        )
        
        # Save filtered data
        logger.info(f"\nSaving filtered data to: {args.output}")
        adata_filtered.write_h5ad(args.output, compression='gzip')
        
        logger.info("\n" + "=" * 70)
        logger.info("✓ Quality Control Complete")
        logger.info("=" * 70)
        
        if not args.skip_filtering:
            logger.info(f"Results:")
            logger.info(f"  Cells: {filtering_stats['cells_before']} → {filtering_stats['cells_after']} "
                       f"({filtering_stats['cells_retained_pct']:.1f}% retained)")
            logger.info(f"  Genes: {filtering_stats['genes_before']} → {filtering_stats['genes_after']} "
                       f"({filtering_stats['genes_retained_pct']:.1f}% retained)")
        
        logger.info(f"\nOutputs:")
        logger.info(f"  Filtered data: {args.output}")
        logger.info(f"  QC plots: {output_dir}/")
        logger.info(f"  QC report: {output_dir}/qc_report.json")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
