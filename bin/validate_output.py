#!/usr/bin/env python3
"""
Validate scAnnex unified H5AD output

This script validates that UNIFY_INPUT produces correctly structured
H5AD files following scAnnex conventions (InitProject.md Section 9).

Usage: python validate_output.py <path_to_unified.h5ad>
"""

import sys
from pathlib import Path
from typing import List, Tuple

try:
    import scanpy as sc
    import pandas as pd
except ImportError as e:
    print(f"ERROR: {e}")
    print("Install with: pip install scanpy pandas")
    sys.exit(1)


def validate_h5ad(h5ad_path: Path) -> Tuple[bool, List[str]]:
    """
    Validate H5AD structure against scAnnex conventions.
    
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    print(f"Validating: {h5ad_path}")
    print("=" * 70)
    
    # Load file
    try:
        adata = sc.read_h5ad(h5ad_path)
        print(f"✓ File loaded successfully")
        print(f"  Dimensions: {adata.n_obs} cells × {adata.n_vars} genes")
    except Exception as e:
        issues.append(f"Failed to load H5AD: {e}")
        return False, issues
    
    # Check .obs required columns
    print("\n1. Checking .obs metadata...")
    required_obs = ['sample_id', 'batch', 'condition']
    for col in required_obs:
        if col in adata.obs.columns:
            unique_vals = adata.obs[col].unique()
            print(f"  ✓ {col}: {unique_vals[0] if len(unique_vals) == 1 else f'{len(unique_vals)} unique values'}")
        else:
            issues.append(f"Missing required .obs column: {col}")
            print(f"  ✗ Missing: {col}")
    
    # Check for cell_id index
    if adata.obs.index.name:
        print(f"  ✓ Cell index name: {adata.obs.index.name}")
    else:
        print(f"  ⚠ Cell index has no name (recommended: 'cell_id')")
    
    # Check .var structure
    print("\n2. Checking .var metadata...")
    if adata.var.index.name:
        print(f"  ✓ Gene index name: {adata.var.index.name}")
    else:
        print(f"  ⚠ Gene index has no name (recommended: 'gene_id' or 'gene_symbols')")
    
    # Check for duplicate genes
    n_duplicates = adata.var.index.duplicated().sum()
    if n_duplicates > 0:
        issues.append(f"Found {n_duplicates} duplicate gene names")
        print(f"  ✗ {n_duplicates} duplicate gene names")
    else:
        print(f"  ✓ No duplicate gene names")
    
    # Check .layers
    print("\n3. Checking .layers...")
    if 'counts' in adata.layers:
        print(f"  ✓ 'counts' layer present")
        # Verify it's raw counts (integers or close to it)
        counts_sum = adata.layers['counts'].sum()
        x_sum = adata.X.sum()
        print(f"    - counts sum: {counts_sum:.0f}")
        print(f"    - X sum: {x_sum:.0f}")
    else:
        issues.append("Missing 'counts' layer for raw counts")
        print(f"  ✗ Missing 'counts' layer")
    
    if len(adata.layers.keys()) > 0:
        print(f"  Available layers: {list(adata.layers.keys())}")
    
    # Check .uns metadata
    print("\n4. Checking .uns metadata...")
    if 'scannex' in adata.uns:
        print(f"  ✓ scannex metadata present")
        if 'conversion' in adata.uns['scannex']:
            conv = adata.uns['scannex']['conversion']
            print(f"    - input_type: {conv.get('input_type', 'N/A')}")
            print(f"    - sample_id: {conv.get('sample_id', 'N/A')}")
            print(f"    - timestamp: {conv.get('timestamp', 'N/A')}")
        if 'versions' in adata.uns['scannex']:
            vers = adata.uns['scannex']['versions']
            print(f"    - scanpy: {vers.get('scanpy', 'N/A')}")
            print(f"    - anndata: {vers.get('anndata', 'N/A')}")
    else:
        print(f"  ⚠ No scannex metadata (recommended)")
    
    # Check .obsm
    print("\n5. Checking .obsm (multi-dimensional annotations)...")
    if len(adata.obsm.keys()) > 0:
        print(f"  Available: {list(adata.obsm.keys())}")
    else:
        print(f"  ⚠ Empty (will be populated during dimensionality reduction)")
    
    # Memory usage estimate
    print("\n6. File information...")
    file_size = h5ad_path.stat().st_size / (1024 ** 2)
    print(f"  File size: {file_size:.2f} MB")
    
    # Summary
    print("\n" + "=" * 70)
    if len(issues) == 0:
        print("✓ VALIDATION PASSED - File meets scAnnex conventions")
        return True, []
    else:
        print(f"✗ VALIDATION FAILED - {len(issues)} issue(s) found:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
        return False, issues


def main():
    if len(sys.argv) != 2:
        print("Usage: python validate_output.py <path_to_unified.h5ad>")
        sys.exit(1)
    
    h5ad_path = Path(sys.argv[1])
    
    if not h5ad_path.exists():
        print(f"ERROR: File not found: {h5ad_path}")
        sys.exit(1)
    
    is_valid, issues = validate_h5ad(h5ad_path)
    
    sys.exit(0 if is_valid else 1)


if __name__ == "__main__":
    main()
