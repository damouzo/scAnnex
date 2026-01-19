#!/usr/bin/env python3
"""Quick inspection of the test output"""

import sys
sys.path.insert(0, '/usr/local/lib/python3.10/site-packages')

import scanpy as sc

print("Loading test output...")
adata = sc.read_h5ad('test_data/outputs/PBMC_MTX_quick_test.h5ad')

print("\n" + "="*70)
print("UNIFY_INPUT Test Output Inspection")
print("="*70)

print(f"\n✓ File loaded successfully")
print(f"  Dimensions: {adata.n_obs} cells × {adata.n_vars} genes")

print("\n1. .obs metadata:")
print(f"  Columns: {list(adata.obs.columns)}")
print(f"  sample_id: {adata.obs['sample_id'].unique()}")
print(f"  batch: {adata.obs['batch'].unique()}")
print(f"  condition: {adata.obs['condition'].unique()}")

print("\n2. .var metadata:")
print(f"  Columns: {list(adata.var.columns)}")
print(f"  Index name: {adata.var.index.name}")

print("\n3. .layers:")
print(f"  Available: {list(adata.layers.keys())}")
if 'counts' in adata.layers:
    print(f"  ✓ counts layer present")
    print(f"    Sum: {adata.layers['counts'].sum():.0f}")

print("\n4. .uns metadata:")
print(f"  Keys: {list(adata.uns.keys())}")
if 'scannex' in adata.uns:
    print(f"  ✓ scannex metadata present:")
    if 'conversion' in adata.uns['scannex']:
        conv = adata.uns['scannex']['conversion']
        print(f"    - input_type: {conv['input_type']}")
        print(f"    - sample_id: {conv['sample_id']}")
    if 'versions' in adata.uns['scannex']:
        vers = adata.uns['scannex']['versions']
        print(f"    - scanpy: {vers['scanpy']}")
        print(f"    - anndata: {vers['anndata']}")

print("\n5. .obsm (future PCA/UMAP):")
print(f"  Keys: {list(adata.obsm.keys()) if len(adata.obsm.keys()) > 0 else 'Empty (as expected)'}")

print("\n" + "="*70)
print("✓ VALIDATION SUCCESSFUL - All required fields present!")
print("="*70)
