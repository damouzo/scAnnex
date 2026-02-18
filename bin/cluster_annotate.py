#!/usr/bin/env python3
"""
Clustering and cell type annotation
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
    parser = argparse.ArgumentParser(description="Cluster and annotate cell types")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output", type=str, required=True, help="Output H5AD file")
    parser.add_argument("--resolution", type=float, default=1.0,
                       help="Clustering resolution")
    parser.add_argument("--markers", type=str, default=None,
                       help="Path to marker genes CSV for annotation")
    return parser.parse_args()


def compute_neighbors_umap(adata):
    """Compute neighbors and UMAP"""
    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    
    print("Computing UMAP...")
    sc.tl.umap(adata)
    
    return adata


def perform_clustering(adata, resolution=1.0):
    """Perform Leiden clustering"""
    print(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)
    
    n_clusters = adata.obs['leiden'].nunique()
    print(f"Identified {n_clusters} clusters")
    
    return adata


def find_marker_genes(adata):
    """Find marker genes for each cluster"""
    print("Finding marker genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Get top markers
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    markers_dict = {}
    for group in groups:
        markers_dict[group] = [
            result['names'][group][i] 
            for i in range(min(10, len(result['names'][group])))
        ]
    
    return markers_dict


def annotate_clusters(adata, marker_file):
    """Annotate clusters using provided marker genes"""
    print(f"Loading marker genes from: {marker_file}")
    markers_df = pd.read_csv(marker_file)
    
    # Simple annotation based on top expressed markers
    # This is a placeholder - more sophisticated annotation can be added
    adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
    
    # If marker file has cluster -> cell_type mapping
    if 'cluster' in markers_df.columns and 'cell_type' in markers_df.columns:
        cluster_map = dict(zip(markers_df['cluster'].astype(str), 
                             markers_df['cell_type']))
        adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_map).fillna('Unknown')
    
    return adata


def plot_results(adata):
    """Generate clustering and annotation plots"""
    # UMAP colored by cluster
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='leiden', legend_loc='on data', ax=ax, show=False)
    plt.savefig("umap_clusters.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP colored by cell type (if annotated)
    if 'cell_type' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        sc.pl.umap(adata, color='cell_type', ax=ax, show=False)
        plt.savefig("umap_celltypes.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Dotplot of top markers
    if 'rank_genes_groups' in adata.uns:
        fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
        plt.savefig("marker_dotplot.png", dpi=300, bbox_inches='tight')
        plt.close()


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Compute neighbors and UMAP
    adata = compute_neighbors_umap(adata)
    
    # Clustering
    adata = perform_clustering(adata, args.resolution)
    
    # Find marker genes
    markers = find_marker_genes(adata)
    
    # Save markers
    markers_df = pd.DataFrame([
        {'cluster': k, 'markers': ', '.join(v)} 
        for k, v in markers.items()
    ])
    markers_df.to_csv("cluster_markers.csv", index=False)
    
    # Annotate if markers provided
    if args.markers:
        adata = annotate_clusters(adata, args.markers)
    
    # Generate plots
    print("Generating plots...")
    plot_results(adata)
    
    # Save results
    print(f"Saving results to: {args.output}")
    adata.write_h5ad(args.output)
    
    print("✓ Clustering and annotation complete")


if __name__ == "__main__":
    main()
