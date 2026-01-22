#!/usr/bin/env python3
"""
Helper script to convert H5AD metadata to R-compatible format
This avoids reticulate's py_to_r() issues with pandas DataFrames
"""

import anndata as ad
import json
import sys

def convert_obs_to_json(h5ad_path, output_path):
    """Convert adata.obs to JSON file that R can read"""
    adata = ad.read_h5ad(h5ad_path)
    
    # Convert to dict of lists
    obs_dict = adata.obs.to_dict('list')
    
    # Convert any non-JSON-serializable types
    for col in obs_dict:
        obs_dict[col] = [str(x) if not isinstance(x, (int, float, bool, str, type(None))) else x 
                         for x in obs_dict[col]]
    
    # Add cell IDs
    obs_dict['cell_id'] = adata.obs_names.tolist()
    
    # Write to JSON
    with open(output_path, 'w') as f:
        json.dump(obs_dict, f)
    
    print(f"Converted {adata.n_obs} cells x {len(obs_dict)-1} columns to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_obs_to_json.py <input.h5ad> <output.json>")
        sys.exit(1)
    
    convert_obs_to_json(sys.argv[1], sys.argv[2])
