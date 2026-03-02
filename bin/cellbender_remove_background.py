#!/usr/bin/env python3

"""
=============================================================================
CellBender Ambient RNA Removal
=============================================================================
Wrapper for CellBender remove-background to eliminate ambient RNA contamination
from 10x Genomics raw count matrices.

Features:
- Auto-detection of expected cells from filtered matrix
- Configurable false positive rate
- GPU support (optional, falls back to CPU)
- H5AD output for seamless integration with scAnnex

Requirements:
- cellbender >= 0.3.0
- scanpy >= 1.9.0
- pytorch >= 1.8.0

Note: CellBender works best with RAW 10x counts (before any filtering).
      Only run on MTX format inputs, not pre-processed H5AD/RDS.

Usage:
    python cellbender_remove_background.py \\
        --input raw_feature_bc_matrix/ \\
        --output cellbender_output.h5ad \\
        --expected-cells auto \\
        --total-droplets 25000 \\
        --fpr 0.01
=============================================================================
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="CellBender ambient RNA removal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect expected cells
  %(prog)s --input raw_mtx/ --output output.h5ad --expected-cells auto
  
  # Specify expected cells
  %(prog)s --input raw_mtx/ --output output.h5ad \\
    --expected-cells 5000 --total-droplets 25000
  
  # Use GPU
  %(prog)s --input raw_mtx/ --output output.h5ad --cuda
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input directory (10x MTX format) or H5 file"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output H5AD file"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="cellbender_tmp",
        help="Temporary output directory for CellBender (default: cellbender_tmp)"
    )
    
    # CellBender parameters
    parser.add_argument(
        "--expected-cells",
        type=str,
        default="auto",
        help="Expected number of cells (auto or integer, default: auto)"
    )
    parser.add_argument(
        "--total-droplets",
        type=int,
        default=25000,
        help="Total droplets including empty (default: 25000)"
    )
    parser.add_argument(
        "--fpr",
        type=float,
        default=0.01,
        help="False positive rate for calling cells (default: 0.01)"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=150,
        help="Training epochs (default: 150)"
    )
    
    # Hardware options
    parser.add_argument(
        "--cuda",
        action="store_true",
        help="Use GPU acceleration (requires CUDA)"
    )
    
    return parser.parse_args()


def detect_expected_cells(input_path: str) -> int:
    """
    Auto-detect expected cells from filtered matrix.
    
    Assumes filtered matrix is in sibling directory or tries to estimate
    from raw counts distribution.
    """
    logger.info("Auto-detecting expected number of cells...")
    
    # Try to find filtered matrix in parent directory
    input_path_obj = Path(input_path)
    parent_dir = input_path_obj.parent
    
    # Common 10x directory names
    filtered_dirs = [
        parent_dir / "filtered_feature_bc_matrix",
        parent_dir / "filtered_gene_bc_matrices",
        parent_dir.parent / "filtered_feature_bc_matrix",
    ]
    
    for filt_dir in filtered_dirs:
        if filt_dir.exists():
            try:
                logger.info(f"Found filtered matrix: {filt_dir}")
                adata_filt = sc.read_10x_mtx(str(filt_dir))
                n_cells = adata_filt.n_obs
                logger.info(f"Auto-detected {n_cells} cells from filtered matrix")
                return n_cells
            except Exception as e:
                logger.warning(f"Failed to read filtered matrix: {e}")
    
    # Fallback: estimate from raw matrix
    logger.info("Filtered matrix not found, estimating from raw counts...")
    try:
        if Path(input_path).is_dir():
            adata_raw = sc.read_10x_mtx(input_path)
        else:
            adata_raw = sc.read_10x_h5(input_path)
        
        # Simple heuristic: cells with >1000 UMI
        umi_counts = adata_raw.X.sum(axis=1).A1
        estimated_cells = (umi_counts > 1000).sum()
        
        logger.info(f"Estimated {estimated_cells} cells (>1000 UMI threshold)")
        return int(estimated_cells)
    
    except Exception as e:
        logger.warning(f"Failed to estimate cells: {e}")
        logger.warning("Using default: 5000 cells")
        return 5000


def run_cellbender(
    input_path: str,
    output_dir: str,
    expected_cells: int,
    total_droplets: int,
    fpr: float,
    epochs: int,
    use_cuda: bool
) -> str:
    """
    Run CellBender remove-background.
    
    Returns:
        Path to CellBender output H5 file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    output_h5 = os.path.join(output_dir, "cellbender_output.h5")
    
    # Build CellBender command
    cmd = [
        "cellbender",
        "remove-background",
        "--input", input_path,
        "--output", output_h5,
        "--expected-cells", str(expected_cells),
        "--total-droplets-included", str(total_droplets),
        "--fpr", str(fpr),
        "--epochs", str(epochs)
    ]
    
    if use_cuda:
        cmd.append("--cuda")
    
    logger.info("Running CellBender remove-background...")
    logger.info(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        logger.info("CellBender completed successfully")
        logger.debug(result.stdout)
    
    except subprocess.CalledProcessError as e:
        logger.error(f"CellBender failed with exit code {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        sys.exit(1)
    
    except FileNotFoundError:
        logger.error("CellBender not found. Please install: pip install cellbender")
        sys.exit(1)
    
    # CellBender appends _filtered.h5 to output name
    cellbender_output = output_h5.replace(".h5", "_filtered.h5")
    
    if not os.path.exists(cellbender_output):
        logger.error(f"CellBender output not found: {cellbender_output}")
        sys.exit(1)
    
    return cellbender_output


def convert_to_h5ad(cellbender_h5: str, output_h5ad: str):
    """Convert CellBender H5 output to H5AD format."""
    logger.info(f"Converting CellBender output to H5AD: {output_h5ad}")
    
    try:
        # Read CellBender output
        adata = sc.read_10x_h5(cellbender_h5)
        
        # Add metadata
        adata.uns['cellbender'] = {
            'version': 'cellbender >= 0.3.0',
            'ambient_rna_removed': True,
            'source_file': cellbender_h5
        }
        
        # Save as H5AD
        adata.write_h5ad(output_h5ad, compression='gzip')
        
        logger.info(f"Saved H5AD: {output_h5ad}")
        logger.info(f"  Cells: {adata.n_obs}")
        logger.info(f"  Genes: {adata.n_vars}")
    
    except Exception as e:
        logger.error(f"Failed to convert to H5AD: {e}")
        sys.exit(1)


def main():
    args = parse_args()
    
    # Validate input
    if not os.path.exists(args.input):
        logger.error(f"Input not found: {args.input}")
        sys.exit(1)
    
    # Auto-detect expected cells if needed
    if args.expected_cells.lower() == 'auto':
        expected_cells = detect_expected_cells(args.input)
    else:
        try:
            expected_cells = int(args.expected_cells)
        except ValueError:
            logger.error(f"Invalid expected-cells value: {args.expected_cells}")
            sys.exit(1)
    
    logger.info("=" * 70)
    logger.info("CellBender Ambient RNA Removal")
    logger.info("=" * 70)
    logger.info(f"Input: {args.input}")
    logger.info(f"Expected cells: {expected_cells}")
    logger.info(f"Total droplets: {args.total_droplets}")
    logger.info(f"FPR threshold: {args.fpr}")
    logger.info(f"Epochs: {args.epochs}")
    logger.info(f"CUDA: {args.cuda}")
    logger.info("=" * 70)
    
    # Run CellBender
    cellbender_h5 = run_cellbender(
        args.input,
        args.output_dir,
        expected_cells,
        args.total_droplets,
        args.fpr,
        args.epochs,
        args.cuda
    )
    
    # Convert to H5AD
    convert_to_h5ad(cellbender_h5, args.output)
    
    logger.info("✓ CellBender ambient RNA removal complete")
    logger.info(f"  Output: {args.output}")


if __name__ == "__main__":
    main()
