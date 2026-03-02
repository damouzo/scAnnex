#!/usr/bin/env python3
"""
Prepare merged H5AD file for DGE testing (skip individual QC steps)
Combines the 4 demo samples into a single AnnData object
"""

import sys
from pathlib import Path

try:
    import scanpy as sc
    import pandas as pd
    import anndata as ad
    import numpy as np
    
    ad.settings.allow_write_nullable_strings = True
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    sys.exit(1)

# Load all 4 samples
samples = [
    'data_demo/MultiSample/H5AD/control_batch1.h5ad',
    'data_demo/MultiSample/H5AD/control_batch2.h5ad',
    'data_demo/MultiSample/H5AD/treated_batch1.h5ad',
    'data_demo/MultiSample/H5AD/treated_batch2.h5ad',
]

print("Loading samples...")
adatas = []
for sample_path in samples:
    adata = sc.read_h5ad(sample_path)
    print(f"  {sample_path}: {adata.n_obs} cells")
    
    # Basic QC filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Calculate QC metrics manually to avoid pandas BooleanArray issue
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Manual QC calculation to avoid BooleanArray bug
    mt_genes_idx = np.where(adata.var['mt'].to_numpy().astype(bool))[0]
    
    # Calculate total counts
    from scipy.sparse import issparse
    if issparse(adata.X):
        adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        adata.obs['n_genes_by_counts'] = np.array((adata.X > 0).sum(axis=1)).flatten()
        if len(mt_genes_idx) > 0:
            adata.obs['total_counts_mt'] = np.array(adata.X[:, mt_genes_idx].sum(axis=1)).flatten()
        else:
            adata.obs['total_counts_mt'] = 0
    else:
        adata.obs['total_counts'] = adata.X.sum(axis=1)
        adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1)
        if len(mt_genes_idx) > 0:
            adata.obs['total_counts_mt'] = adata.X[:, mt_genes_idx].sum(axis=1)
        else:
            adata.obs['total_counts_mt'] = 0
    
    adata.obs['pct_counts_mt'] = (adata.obs['total_counts_mt'] / adata.obs['total_counts']) * 100
    
    # Filter by mitochondrial content
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    
    adatas.append(adata)
    print(f"    After QC: {adata.n_obs} cells")

# Merge all samples
print("\nMerging samples...")
adata_merged = ad.concat(adatas, join='outer', merge='same')
adata_merged.obs_names_make_unique()

print(f"Merged dataset: {adata_merged.n_obs} cells × {adata_merged.n_vars} genes")

# Normalize and process
print("\nNormalizing and processing...")
sc.pp.normalize_total(adata_merged, target_sum=1e4)
sc.pp.log1p(adata_merged)
sc.pp.highly_variable_genes(adata_merged, n_top_genes=2000)

# Run PCA
print("Running PCA...")
sc.pp.pca(adata_merged, n_comps=50)

# Run Harmony integration
print("Running Harmony integration...")
try:
    import harmonypy as hm
    ho = hm.run_harmony(adata_merged.obsm['X_pca'], adata_merged.obs, 'batch', max_iter_harmony=20)
    adata_merged.obsm['X_pca_harmony'] = ho.Z_corr.T
    print("  Harmony completed")
except ImportError:
    print("  WARNING: Harmony not available, skipping batch correction")
    adata_merged.obsm['X_pca_harmony'] = adata_merged.obsm['X_pca'].copy()

# Compute neighbors and UMAP on integrated space
print("Computing neighbors and UMAP...")
sc.pp.neighbors(adata_merged, use_rep='X_pca_harmony', n_neighbors=15)
sc.tl.umap(adata_merged)

# Clustering (skip if igraph not available - not needed for DGE testing)
print("Running clustering...")
try:
    sc.tl.leiden(adata_merged, resolution=0.5, key_added='leiden_0.5')
    sc.tl.leiden(adata_merged, resolution=1.0, key_added='leiden_1.0')
    print("  Leiden clustering completed")
except (ImportError, ModuleNotFoundError):
    print("  Clustering skipped (igraph not available - not needed for DGE)")

print("\nSample distribution:")
print(adata_merged.obs.groupby(['batch', 'condition']).size())

# Save
output_file = 'data_demo/MultiSample/merged_processed.h5ad'
print(f"\nSaving to {output_file}...")
adata_merged.write_h5ad(output_file)

print("✓ Done!")
print(f"\nTo test DGE, run:")
print(f"  python bin/differential_expression.py \\")
print(f"    --input {output_file} \\")
print(f"    --output test_dge_results.h5ad \\")
print(f"    --output-dir test_dge_results \\")
print(f"    --contrasts data_demo/MultiSample/contrasts_example.csv \\")
print(f"    --save-plots")
