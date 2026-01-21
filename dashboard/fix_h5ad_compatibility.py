#!/usr/bin/env python3
"""
Fix H5AD file compatibility for older anndata versions
Removes unsupported 'null' encoding from /uns/log1p/base
"""

import h5py
import sys
import shutil
from pathlib import Path

def fix_h5ad_compatibility(h5ad_path):
    """Fix encoding issues in H5AD file for anndata 0.11.x compatibility"""
    
    h5ad_path = Path(h5ad_path)
    backup_path = h5ad_path.with_suffix('.h5ad.backup')
    
    print(f"Fixing H5AD file: {h5ad_path}")
    
    # Create backup
    print(f"Creating backup: {backup_path}")
    shutil.copy2(h5ad_path, backup_path)
    
    try:
        with h5py.File(h5ad_path, 'r+') as f:
            # Check if the problematic encoding exists
            if 'uns/log1p/base' in f:
                print("Found /uns/log1p/base with 'null' encoding")
                
                # Read the value
                base_value = f['uns/log1p/base'][()]
                print(f"  Current value: {base_value}")
                
                # Delete the dataset
                del f['uns/log1p/base']
                print("  Deleted problematic dataset")
                
                # Recreate as a simple scalar without custom encoding
                f['uns/log1p'].create_dataset('base', data=base_value, dtype='float64')
                print(f"  Recreated as simple float64: {base_value}")
                
                # Remove the encoding attributes
                if 'encoding-type' in f['uns/log1p/base'].attrs:
                    del f['uns/log1p/base'].attrs['encoding-type']
                if 'encoding-version' in f['uns/log1p/base'].attrs:
                    del f['uns/log1p/base'].attrs['encoding-version']
                
                print("✓ Fixed successfully!")
                print(f"\nBackup saved at: {backup_path}")
                print("You can delete the backup if everything works correctly.")
                
            else:
                print("No problematic encoding found. File should be compatible.")
                
    except Exception as e:
        print(f"\n✗ Error fixing file: {e}")
        print(f"Restoring from backup...")
        shutil.copy2(backup_path, h5ad_path)
        raise
    
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fix_h5ad_compatibility.py <path/to/file.h5ad>")
        sys.exit(1)
    
    h5ad_file = sys.argv[1]
    
    if not Path(h5ad_file).exists():
        print(f"Error: File not found: {h5ad_file}")
        sys.exit(1)
    
    fix_h5ad_compatibility(h5ad_file)
