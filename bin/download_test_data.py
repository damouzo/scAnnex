#!/usr/bin/env python3
"""
Download Test Data for scAnnex SLC Pipeline

Downloads 2 small PBMC datasets from 10x Genomics for testing
the complete SLC (Simple, Lovable, Complete) pipeline.

Downloads:
- PBMC 1k (v3 chemistry) - batch1
- PBMC 1k (v3 chemistry) - batch2 (same dataset, simulated as batch)

Author: scAnnex SLC Development Team
"""

import argparse
import gzip
import shutil
import sys
from pathlib import Path
from urllib.request import urlretrieve

def download_file(url: str, dest: Path, desc: str = "file"):
    """Download a file with progress."""
    print(f"  Downloading {desc}...")
    print(f"    URL: {url}")
    
    def progress(block_num, block_size, total_size):
        """Show download progress."""
        downloaded = block_num * block_size
        if total_size > 0:
            percent = downloaded * 100 / total_size
            sys.stdout.write(f"\r    Progress: {percent:.1f}% ({downloaded}/{total_size} bytes)")
            sys.stdout.flush()
    
    try:
        urlretrieve(url, dest, reporthook=progress)
        print(f"\n    ✓ Saved to: {dest}")
    except Exception as e:
        print(f"\n    ✗ ERROR: {e}")
        sys.exit(1)


def setup_test_data(output_dir: Path):
    """Download and setup test data."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("scAnnex SLC Test Data Download")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}\n")
    
    # 10x Genomics PBMC 1k v3 dataset
    base_url = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3"
    
    # Download H5 format (simplest for testing)
    h5_file = output_dir / "pbmc_1k_v3_filtered_feature_bc_matrix.h5"
    
    if h5_file.exists():
        print(f"✓ Dataset already exists: {h5_file}")
    else:
        print("Downloading PBMC 1k v3 dataset (H5 format)...")
        download_file(
            f"{base_url}/pbmc_1k_v3_filtered_feature_bc_matrix.h5",
            h5_file,
            "PBMC 1k v3 H5 file"
        )
    
    # Create batch1 and batch2 directories (symlink for simplicity)
    batch1_dir = output_dir / "batch1"
    batch2_dir = output_dir / "batch2"
    
    batch1_dir.mkdir(exist_ok=True)
    batch2_dir.mkdir(exist_ok=True)
    
    batch1_h5 = batch1_dir / "pbmc_1k_batch1.h5"
    batch2_h5 = batch2_dir / "pbmc_1k_batch2.h5"
    
    # Copy files (not symlink to avoid issues)
    if not batch1_h5.exists():
        print(f"\nCreating batch1 dataset...")
        shutil.copy(h5_file, batch1_h5)
        print(f"  ✓ Created: {batch1_h5}")
    
    if not batch2_h5.exists():
        print(f"Creating batch2 dataset...")
        shutil.copy(h5_file, batch2_h5)
        print(f"  ✓ Created: {batch2_h5}")
    
    # Create samplesheet
    samplesheet = output_dir / "samplesheet.csv"
    
    print(f"\nCreating samplesheet: {samplesheet}")
    with open(samplesheet, 'w') as f:
        f.write("sample_id,file_type,file_path,batch,condition\n")
        f.write(f"PBMC_batch1,h5,{batch1_h5.absolute()},batch1,control\n")
        f.write(f"PBMC_batch2,h5,{batch2_h5.absolute()},batch2,control\n")
    
    print(f"  ✓ Created: {samplesheet}")
    
    # Summary
    print("\n" + "=" * 70)
    print("✓ Test Data Setup Complete")
    print("=" * 70)
    print(f"\nDatasets:")
    print(f"  - Batch 1: {batch1_h5}")
    print(f"  - Batch 2: {batch2_h5}")
    print(f"\nSamplesheet:")
    print(f"  {samplesheet}")
    print(f"\nTo run the SLC pipeline:")
    print(f"  nextflow run main.nf --input {samplesheet} --outdir results/")
    print("=" * 70)


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description="Download test data for scAnnex SLC pipeline"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="test_data",
        help="Output directory for test data (default: test_data)"
    )
    
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    
    setup_test_data(output_dir)


if __name__ == "__main__":
    main()
