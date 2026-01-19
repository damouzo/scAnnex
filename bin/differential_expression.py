#!/usr/bin/env python3
"""
Differential expression analysis
"""

import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description="Differential expression analysis")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--groupby", type=str, required=True,
                       help="Column in obs to group by (e.g., 'leiden', 'cell_type')")
    parser.add_argument("--method", type=str, default="wilcoxon",
                       choices=["wilcoxon", "t-test", "logreg"],
                       help="DE test method")
    parser.add_argument("--reference", type=str, default="rest",
                       help="Reference group for comparison")
    return parser.parse_args()


def run_differential_expression(adata, groupby, method="wilcoxon", reference="rest"):
    """Run differential expression analysis"""
    print(f"Running differential expression ({method})...")
    print(f"Grouping by: {groupby}")
    
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        reference=reference,
        use_raw=False
    )
    
    return adata


def extract_de_results(adata, groupby):
    """Extract DE results into DataFrame"""
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    de_list = []
    for group in groups:
        genes = result['names'][group]
        scores = result['scores'][group]
        pvals = result['pvals'][group]
        pvals_adj = result['pvals_adj'][group]
        logfoldchanges = result['logfoldchanges'][group]
        
        for i in range(len(genes)):
            de_list.append({
                'group': group,
                'gene': genes[i],
                'score': scores[i],
                'pval': pvals[i],
                'pval_adj': pvals_adj[i],
                'log2fc': logfoldchanges[i]
            })
    
    de_df = pd.DataFrame(de_list)
    return de_df


def plot_de_results(adata, groupby, output_dir):
    """Generate DE plots"""
    # Dotplot of top DE genes
    fig = sc.pl.rank_genes_groups_dotplot(
        adata, 
        n_genes=10, 
        show=False,
        dendrogram=True
    )
    plt.savefig(f"{output_dir}/de_dotplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Heatmap of top DE genes
    fig = sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=10,
        show=False,
        dendrogram=True,
        swap_axes=True
    )
    plt.savefig(f"{output_dir}/de_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Violin plot for top genes
    result = adata.uns['rank_genes_groups']
    top_genes = [result['names'][g][0] for g in result['names'].dtype.names[:5]]
    
    fig, axes = plt.subplots(1, len(top_genes), figsize=(4*len(top_genes), 4))
    sc.pl.violin(adata, keys=top_genes, groupby=groupby, ax=axes, show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/de_violin.png", dpi=300, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Check groupby column exists
    if args.groupby not in adata.obs.columns:
        raise ValueError(f"Column '{args.groupby}' not found in adata.obs")
    
    # Run DE analysis
    adata = run_differential_expression(adata, args.groupby, args.method, args.reference)
    
    # Extract results
    de_df = extract_de_results(adata, args.groupby)
    
    # Save results
    import os
    os.makedirs(args.output_dir, exist_ok=True)
    de_df.to_csv(f"{args.output_dir}/de_results.csv", index=False)
    
    # Save top genes per group
    top_genes = de_df.groupby('group').apply(
        lambda x: x.nsmallest(50, 'pval_adj')
    ).reset_index(drop=True)
    top_genes.to_csv(f"{args.output_dir}/top_de_genes.csv", index=False)
    
    # Generate plots
    print("Generating plots...")
    plot_de_results(adata, args.groupby, args.output_dir)
    
    print(f"✓ Differential expression analysis complete")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
