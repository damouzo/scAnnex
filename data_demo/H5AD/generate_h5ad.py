#!/usr/bin/env python3
"""
Generate H5AD demo file from MTX data for scAnnex
Creates a 1k cell H5AD file with minimal preprocessing
"""

import sys
from pathlib import Path

try:
    import scanpy as sc
    import pandas as pd
    import anndata as ad
    
    # Enable writing of nullable strings (required for anndata >= 0.11)
    ad.settings.allow_write_nullable_strings = True
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Please install: pip install scanpy pandas")
    sys.exit(1)

# Paths
mtx_dir = Path(__file__).parent / "10xMTX" / "filtered_feature_bc_matrix"
output_file = Path(__file__).parent / "H5AD" / "pbmc_1k.h5ad"

print("Generating H5AD demo file...")
print(f"Input: {mtx_dir}")

# Load MTX data
if not mtx_dir.exists():
    print(f"ERROR: MTX directory not found: {mtx_dir}")
    sys.exit(1)

try:
    adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', make_unique=True)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
except Exception as e:
    print(f"ERROR loading MTX: {e}")
    sys.exit(1)

# Add basic metadata
adata.obs['sample_id'] = 'PBMC_1k'
adata.obs['batch'] = 'batch1'
adata.obs['condition'] = 'control'

# Store raw counts in layers
adata.layers['counts'] = adata.X.copy()

# Add metadata to .uns
adata.uns['demo_data'] = {
    'source': 'PBMC 1k v3 - 10x Genomics',
    'n_cells': adata.n_obs,
    'created_for': 'scAnnex demo'
}

# Save H5AD
output_file.parent.mkdir(parents=True, exist_ok=True)
try:
    adata.write_h5ad(output_file, compression='gzip')
    file_size = output_file.stat().st_size / (1024 * 1024)
    print(f"✓ Created: {output_file}")
    print(f"✓ Size: {file_size:.2f} MB")
    print(f"✓ Dimensions: {adata.n_obs} cells × {adata.n_vars} genes")
except Exception as e:
    print(f"ERROR saving H5AD: {e}")
    sys.exit(1)

print("\n✓ H5AD demo file created successfully!")
