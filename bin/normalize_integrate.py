#!/usr/bin/env python3
"""
Normalization and batch integration
"""

import argparse
import scanpy as sc
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize and integrate data")
    parser.add_argument("--input", type=str, required=True, help="Input H5AD file")
    parser.add_argument("--output", type=str, required=True, help="Output H5AD file")
    parser.add_argument("--method", type=str, default="log", choices=["log", "scran"],
                       help="Normalization method")
    parser.add_argument("--target-sum", type=float, default=1e4,
                       help="Target sum for normalization")
    parser.add_argument("--batch-key", type=str, default=None,
                       help="Column in obs for batch correction")
    parser.add_argument("--integration-method", type=str, default="harmony",
                       choices=["harmony", "scanorama", "bbknn"],
                       help="Batch integration method")
    return parser.parse_args()


def normalize_data(adata, method="log", target_sum=1e4):
    """Normalize count data"""
    print(f"Normalizing with method: {method}")
    
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()
    
    if method == "log":
        # Standard log normalization
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    elif method == "scran":
        # Scran normalization (requires scran package)
        try:
            import scanpy.external as sce
            sce.pp.scran_normalization(adata)
        except ImportError:
            print("Warning: scran not available, falling back to log normalization")
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
    
    return adata


def integrate_batches(adata, batch_key, method="harmony"):
    """Perform batch correction"""
    print(f"Running batch integration with: {method}")
    
    # First perform PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)
    
    if method == "harmony":
        try:
            import scanpy.external as sce
            sce.pp.harmony_integrate(adata, batch_key)
        except ImportError:
            print("Warning: harmonypy not available, skipping integration")
            
    elif method == "scanorama":
        try:
            import scanpy.external as sce
            sce.pp.scanorama_integrate(adata, batch_key)
        except ImportError:
            print("Warning: scanorama not available, skipping integration")
            
    elif method == "bbknn":
        try:
            import bbknn
            bbknn.bbknn(adata, batch_key=batch_key)
        except ImportError:
            print("Warning: bbknn not available, skipping integration")
    
    return adata


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Normalize
    adata = normalize_data(adata, args.method, args.target_sum)
    
    # Batch integration (optional)
    if args.batch_key:
        if args.batch_key not in adata.obs.columns:
            print(f"Warning: batch key '{args.batch_key}' not found in obs")
        else:
            adata = integrate_batches(adata, args.batch_key, args.integration_method)
    else:
        # Still compute HVGs and PCA for downstream analysis
        print("Computing highly variable genes and PCA...")
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=50)
    
    # Save normalized data
    print(f"Saving normalized data to: {args.output}")
    adata.write_h5ad(args.output)
    
    print("✓ Normalization and integration complete")


if __name__ == "__main__":
    main()
