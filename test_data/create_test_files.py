#!/usr/bin/env python3
"""
Create test files for scAnnex UNIFY_INPUT testing

This script:
1. Loads the PBMC 1k MTX dataset
2. Creates a minimal H5AD test file (100 cells)
3. Validates the files can be loaded
"""

import sys
from pathlib import Path

try:
    import scanpy as sc
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install: pip install scanpy numpy pandas")
    sys.exit(1)

# Set paths
mtx_dir = Path(__file__).parent / "mtx" / "filtered_feature_bc_matrix"
h5ad_dir = Path(__file__).parent / "h5ad"

print("=" * 60)
print("Creating scAnnex Test Files")
print("=" * 60)

# Load MTX data
print(f"\n1. Loading MTX data from: {mtx_dir}")
if not mtx_dir.exists():
    print(f"ERROR: MTX directory not found: {mtx_dir}")
    sys.exit(1)

try:
    adata_full = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', make_unique=True)
    print(f"   ✓ Loaded: {adata_full.n_obs} cells × {adata_full.n_vars} genes")
except Exception as e:
    print(f"ERROR loading MTX: {e}")
    sys.exit(1)

# Create minimal H5AD (100 cells subset)
print("\n2. Creating minimal H5AD test file (100 cells)...")
n_cells = min(100, adata_full.n_obs)
cell_indices = np.random.choice(adata_full.n_obs, size=n_cells, replace=False)
adata_small = adata_full[cell_indices, :].copy()

# Add some basic metadata
adata_small.obs['sample_id'] = 'test_sample'
adata_small.obs['batch'] = 'batch1'
adata_small.obs['condition'] = 'control'
adata_small.obs['n_counts'] = adata_small.X.sum(axis=1).A1
adata_small.obs['n_genes'] = (adata_small.X > 0).sum(axis=1).A1

# Store raw counts in layers
adata_small.layers['counts'] = adata_small.X.copy()

# Add some metadata to .uns
adata_small.uns['test_data'] = {
    'source': 'PBMC 1k v3 subset',
    'n_cells_original': adata_full.n_obs,
    'subset_size': n_cells,
    'created_for': 'scAnnex testing'
}

# Save H5AD
h5ad_output = h5ad_dir / "pbmc_100cells.h5ad"
h5ad_output.parent.mkdir(parents=True, exist_ok=True)

try:
    adata_small.write_h5ad(h5ad_output, compression='gzip')
    file_size = h5ad_output.stat().st_size / 1024
    print(f"   ✓ Created: {h5ad_output}")
    print(f"   ✓ Size: {file_size:.2f} KB")
    print(f"   ✓ Dimensions: {adata_small.n_obs} cells × {adata_small.n_vars} genes")
except Exception as e:
    print(f"ERROR saving H5AD: {e}")
    sys.exit(1)

# Verify we can read it back
print("\n3. Verifying H5AD can be reloaded...")
try:
    adata_test = sc.read_h5ad(h5ad_output)
    print(f"   ✓ Successfully reloaded: {adata_test.n_obs} cells × {adata_test.n_vars} genes")
    print(f"   ✓ .obs columns: {list(adata_test.obs.columns)}")
    print(f"   ✓ .layers: {list(adata_test.layers.keys())}")
    print(f"   ✓ .uns keys: {list(adata_test.uns.keys())}")
except Exception as e:
    print(f"ERROR reloading H5AD: {e}")
    sys.exit(1)

print("\n" + "=" * 60)
print("✓ Test file creation complete!")
print("=" * 60)
print(f"\nCreated files:")
print(f"  - {h5ad_output}")
print(f"  - {mtx_dir} (existing)")
print(f"\nNote: RDS file creation requires R with Seurat installed.")
print("      You can create it separately using the create_test_rds.R script.")
