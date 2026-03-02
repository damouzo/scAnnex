#!/usr/bin/env python3

"""
=============================================================================
Differential Gene Expression Analysis with Contrasts Support
=============================================================================
Performs pairwise differential expression analysis using scanpy.

Features:
- Contrasts.csv support for multiple comparisons
- Simple mode (auto-generate contrasts from metadata)
- Multiple statistical methods (Wilcoxon, t-test)
- Per-contrast filtering by cell subsets
- Volcano plots and summary tables

Requirements:
- scanpy >= 1.9.0
- pandas >= 1.5.0
- matplotlib >= 3.5.0
- seaborn >= 0.12.0

Usage:
    # With contrasts file
    python differential_expression.py \\
        --input integrated.h5ad \\
        --output-dir dge_results \\
        --contrasts-file contrasts.csv \\
        --method wilcoxon
    
    # Simple mode (auto-generate contrasts)
    python differential_expression.py \\
        --input integrated.h5ad \\
        --output-dir dge_results \\
        --groupby condition \\
        --reference control \\
        --method wilcoxon
=============================================================================
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Differential gene expression analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Contrasts mode
  %(prog)s --input integrated.h5ad --output-dir dge_results \\
    --contrasts-file contrasts.csv --method wilcoxon
  
  # Simple mode
  %(prog)s --input integrated.h5ad --output-dir dge_results \\
    --groupby condition --reference control --method t-test
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
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for results"
    )
    
    # Contrasts mode
    parser.add_argument(
        "--contrasts-file",
        type=str,
        default=None,
        help="Path to contrasts.csv file (optional)"
    )
    
    # Simple mode (if no contrasts file)
    parser.add_argument(
        "--groupby",
        type=str,
        default='condition',
        help="Column in obs to group by (default: condition)"
    )
    parser.add_argument(
        "--reference",
        type=str,
        default='control',
        help="Reference group for comparison (default: control)"
    )
    
    # Statistical parameters
    parser.add_argument(
        "--method",
        type=str,
        default="wilcoxon",
        choices=["wilcoxon", "t-test", "t-test_overestim_var"],
        help="DE test method (default: wilcoxon)"
    )
    parser.add_argument(
        "--min-pct",
        type=float,
        default=0.1,
        help="Minimum fraction of cells expressing gene (default: 0.1)"
    )
    parser.add_argument(
        "--logfc-threshold",
        type=float,
        default=0.25,
        help="Minimum log fold-change threshold (default: 0.25)"
    )
    parser.add_argument(
        "--pval-cutoff",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for reporting (default: 0.05)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Number of top genes to report per contrast (default: 100)"
    )
    
    # Output options
    parser.add_argument(
        "--save-plots",
        action="store_true",
        default=True,
        help="Generate volcano and MA plots (default: True)"
    )
    parser.add_argument(
        "--no-save-plots",
        action="store_false",
        dest="save_plots",
        help="Skip plot generation"
    )
    
    return parser.parse_args()


def load_contrasts_file(contrasts_file: str) -> pd.DataFrame:
    """Load contrasts from CSV file."""
    logger.info(f"Loading contrasts from: {contrasts_file}")
    
    if not Path(contrasts_file).exists():
        logger.error(f"Contrasts file not found: {contrasts_file}")
        sys.exit(1)
    
    contrasts_df = pd.read_csv(contrasts_file)
    
    # Validate required columns
    required_cols = ['contrast_id', 'variable', 'group1', 'group2']
    missing_cols = [col for col in required_cols if col not in contrasts_df.columns]
    if missing_cols:
        logger.error(f"Missing required columns in contrasts file: {missing_cols}")
        sys.exit(1)
    
    logger.info(f"Loaded {len(contrasts_df)} contrasts")
    return contrasts_df


def auto_generate_contrasts(adata, groupby: str, reference: str) -> pd.DataFrame:
    """Auto-generate contrasts from metadata."""
    logger.info(f"Auto-generating contrasts from column: {groupby}")
    
    if groupby not in adata.obs.columns:
        logger.error(f"Column '{groupby}' not found in adata.obs")
        sys.exit(1)
    
    unique_groups = adata.obs[groupby].unique()
    unique_groups = [g for g in unique_groups if str(g) != 'nan']
    
    if reference not in unique_groups:
        logger.error(
            f"Reference group '{reference}' not found in column '{groupby}'. "
            f"Available groups: {unique_groups}"
        )
        sys.exit(1)
    
    # Generate contrasts: each group vs reference
    contrasts = []
    for group in unique_groups:
        if group != reference:
            contrasts.append({
                'contrast_id': f"{group}_vs_{reference}",
                'variable': groupby,
                'group1': group,
                'group2': reference,
                'filter_column': '',
                'filter_value': ''
            })
    
    contrasts_df = pd.DataFrame(contrasts)
    logger.info(f"Generated {len(contrasts_df)} contrasts: {contrasts_df['contrast_id'].tolist()}")
    
    return contrasts_df


def run_single_contrast(
    adata,
    contrast_id: str,
    variable: str,
    group1: str,
    group2: str,
    filter_column: Optional[str] = None,
    filter_value: Optional[str] = None,
    method: str = "wilcoxon"
) -> Tuple[pd.DataFrame, int, int]:
    """
    Run differential expression for a single contrast.
    
    Returns:
        Tuple of (results_df, n_group1_cells, n_group2_cells)
    """
    logger.info(f"Running contrast: {contrast_id}")
    logger.info(f"  Comparing: {variable} = {group1} vs {group2}")
    
    # Subset data if filter specified
    adata_subset = adata
    if filter_column and pd.notna(filter_column) and filter_value and pd.notna(filter_value):
        logger.info(f"  Filtering: {filter_column} = {filter_value}")
        mask = adata.obs[filter_column].astype(str) == str(filter_value)
        adata_subset = adata[mask].copy()
        logger.info(f"  Filtered to {adata_subset.n_obs} cells")
    
    # Check if variable and groups exist
    if variable not in adata_subset.obs.columns:
        logger.error(f"Variable '{variable}' not found in adata.obs")
        return pd.DataFrame(), 0, 0
    
    # Count cells in each group
    mask_group1 = adata_subset.obs[variable].astype(str) == str(group1)
    mask_group2 = adata_subset.obs[variable].astype(str) == str(group2)
    n_group1 = mask_group1.sum()
    n_group2 = mask_group2.sum()
    
    logger.info(f"  Group1 ({group1}): {n_group1} cells")
    logger.info(f"  Group2 ({group2}): {n_group2} cells")
    
    if n_group1 < 3 or n_group2 < 3:
        logger.warning(f"Skipping contrast {contrast_id}: insufficient cells (need ≥3 per group)")
        return pd.DataFrame(), n_group1, n_group2
    
    # Subset to only these two groups
    mask_both = mask_group1 | mask_group2
    adata_contrast = adata_subset[mask_both].copy()
    
    # Run DE test
    try:
        sc.tl.rank_genes_groups(
            adata_contrast,
            groupby=variable,
            groups=[str(group1)],
            reference=str(group2),
            method=method,
            use_raw=False,
            key_added='rank_genes_groups_temp'
        )
    except Exception as e:
        logger.error(f"DE test failed for contrast {contrast_id}: {e}")
        return pd.DataFrame(), n_group1, n_group2
    
    # Extract results
    result = adata_contrast.uns['rank_genes_groups_temp']
    group_key = str(group1)
    
    if group_key not in result['names'].dtype.names:
        logger.error(f"Group {group1} not found in results")
        return pd.DataFrame(), n_group1, n_group2
    
    # Build results DataFrame
    de_results = []
    genes = result['names'][group_key]
    scores = result['scores'][group_key]
    pvals = result['pvals'][group_key]
    pvals_adj = result['pvals_adj'][group_key]
    logfoldchanges = result['logfoldchanges'][group_key]
    
    for i in range(len(genes)):
        de_results.append({
            'contrast_id': contrast_id,
            'gene': genes[i],
            'log2fc': logfoldchanges[i],
            'pval': pvals[i],
            'pval_adj': pvals_adj[i],
            'score': scores[i],
            'group1': group1,
            'group2': group2,
            'n_group1': n_group1,
            'n_group2': n_group2
        })
    
    de_df = pd.DataFrame(de_results)
    logger.info(f"  Identified {len(de_df)} genes")
    
    return de_df, n_group1, n_group2


def plot_volcano(
    de_df: pd.DataFrame,
    contrast_id: str,
    output_dir: str,
    pval_cutoff: float = 0.05,
    logfc_threshold: float = 0.25
):
    """Generate volcano plot for a contrast."""
    if de_df.empty:
        logger.warning(f"Skipping volcano plot for {contrast_id}: no results")
        return
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Calculate -log10(pval)
    de_df['neg_log10_pval'] = -np.log10(de_df['pval_adj'] + 1e-300)
    
    # Classify genes
    de_df['sig'] = 'Not Sig'
    de_df.loc[
        (de_df['pval_adj'] < pval_cutoff) & (de_df['log2fc'] > logfc_threshold),
        'sig'
    ] = f'Up in {de_df["group1"].iloc[0]}'
    de_df.loc[
        (de_df['pval_adj'] < pval_cutoff) & (de_df['log2fc'] < -logfc_threshold),
        'sig'
    ] = f'Down in {de_df["group1"].iloc[0]}'
    
    # Colors
    colors = {
        'Not Sig': 'lightgray',
        f'Up in {de_df["group1"].iloc[0]}': 'red',
        f'Down in {de_df["group1"].iloc[0]}': 'blue'
    }
    
    # Plot
    for sig_type, color in colors.items():
        subset = de_df[de_df['sig'] == sig_type]
        ax.scatter(
            subset['log2fc'],
            subset['neg_log10_pval'],
            c=color,
            alpha=0.6,
            s=10,
            label=f"{sig_type} ({len(subset)})"
        )
    
    # Add threshold lines
    ax.axhline(-np.log10(pval_cutoff), color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(logfc_threshold, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(-logfc_threshold, color='gray', linestyle='--', linewidth=0.8)
    
    # Labels
    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
    ax.set_title(f'Volcano Plot: {contrast_id}', fontsize=14, fontweight='bold')
    ax.legend(loc='best', frameon=True, fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{contrast_id}_volcano.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/{contrast_id}_volcano.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved volcano plot: {contrast_id}_volcano.pdf")


def save_contrast_summary(
    all_results: List[pd.DataFrame],
    output_dir: str,
    top_n: int = 100,
    pval_cutoff: float = 0.05,
    logfc_threshold: float = 0.25
):
    """Save summary tables for all contrasts."""
    if not all_results:
        logger.warning("No results to save")
        return
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Save full results
    full_output = f"{output_dir}/all_contrasts_full_results.csv"
    combined_df.to_csv(full_output, index=False)
    logger.info(f"Saved full results: {full_output}")
    
    # Save significant genes only
    sig_df = combined_df[
        (combined_df['pval_adj'] < pval_cutoff) &
        (np.abs(combined_df['log2fc']) > logfc_threshold)
    ].copy()
    sig_output = f"{output_dir}/all_contrasts_significant_genes.csv"
    sig_df.to_csv(sig_output, index=False)
    logger.info(f"Saved significant genes: {sig_output} ({len(sig_df)} genes)")
    
    # Save top N genes per contrast
    top_df = combined_df.groupby('contrast_id').apply(
        lambda x: x.nsmallest(top_n, 'pval_adj')
    ).reset_index(drop=True)
    top_output = f"{output_dir}/all_contrasts_top{top_n}_genes.csv"
    top_df.to_csv(top_output, index=False)
    logger.info(f"Saved top {top_n} genes per contrast: {top_output}")
    
    # Per-contrast files
    for contrast_id in combined_df['contrast_id'].unique():
        contrast_df = combined_df[combined_df['contrast_id'] == contrast_id]
        contrast_output = f"{output_dir}/{contrast_id}_results.csv"
        contrast_df.to_csv(contrast_output, index=False)
        
        # Top genes for this contrast
        contrast_top = contrast_df.nsmallest(top_n, 'pval_adj')
        contrast_top_output = f"{output_dir}/{contrast_id}_top{top_n}.csv"
        contrast_top.to_csv(contrast_top_output, index=False)
    
    logger.info(f"Saved per-contrast result files")


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load H5AD
    logger.info(f"Loading H5AD file: {args.input}")
    try:
        adata = sc.read_h5ad(args.input)
        logger.info(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    except Exception as e:
        logger.error(f"Failed to load H5AD file: {e}")
        sys.exit(1)
    
    # Load or generate contrasts
    if args.contrasts_file:
        contrasts_df = load_contrasts_file(args.contrasts_file)
    else:
        contrasts_df = auto_generate_contrasts(adata, args.groupby, args.reference)
    
    # Run DE for each contrast
    all_results = []
    for idx, row in contrasts_df.iterrows():
        contrast_id = row['contrast_id']
        variable = row['variable']
        group1 = row['group1']
        group2 = row['group2']
        filter_column = row.get('filter_column', None)
        filter_value = row.get('filter_value', None)
        
        # Run contrast
        de_df, n_group1, n_group2 = run_single_contrast(
            adata,
            contrast_id,
            variable,
            group1,
            group2,
            filter_column,
            filter_value,
            args.method
        )
        
        if not de_df.empty:
            all_results.append(de_df)
            
            # Generate volcano plot
            if args.save_plots:
                plot_volcano(
                    de_df,
                    contrast_id,
                    args.output_dir,
                    args.pval_cutoff,
                    args.logfc_threshold
                )
    
    # Save summary tables
    save_contrast_summary(
        all_results,
        args.output_dir,
        args.top_n,
        args.pval_cutoff,
        args.logfc_threshold
    )
    
    logger.info(f"✓ Differential expression analysis complete")
    logger.info(f"  Processed {len(all_results)}/{len(contrasts_df)} contrasts successfully")
    logger.info(f"  Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
