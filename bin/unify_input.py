#!/usr/bin/env python3
"""
Unify input formats (H5AD, RDS, MTX) to standardized AnnData H5AD

This script converts various single-cell RNA-seq input formats into a 
standardized AnnData object following the scAnnex conventions defined in 
InitProject.md Section 9.

Features:
- H5AD: Direct loading with validation
- RDS: Two-step conversion via SeuratDisk (RDS → h5seurat → h5ad)
- MTX: 10x Genomics format loading
- Metadata integration: sample_id, batch, condition
- AnnData standardization: .obs, .var, .uns, .layers conventions
- Comprehensive error handling and logging

Author: scAnnex Development Team
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Convert various scRNA-seq formats to unified H5AD",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert H5AD
  %(prog)s --input sample.h5ad --input-type h5ad --output unified.h5ad \\
    --sample-id PBMC_1 --batch batch1 --condition control

  # Convert RDS (requires R with Seurat/SeuratDisk)
  %(prog)s --input seurat.rds --input-type rds --output unified.h5ad \\
    --sample-id PBMC_2 --batch batch1 --condition treated

  # Convert 10x MTX
  %(prog)s --input /path/to/filtered_feature_bc_matrix/ --input-type mtx \\
    --output unified.h5ad --sample-id PBMC_3 --batch batch2 --condition control
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--input", 
        type=str, 
        required=True,
        help="Input file or directory (H5AD, RDS, or MTX directory)"
    )
    parser.add_argument(
        "--input-type",
        type=str,
        required=True,
        choices=["h5ad", "rds", "mtx"],
        help="Type of input data"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output H5AD file path"
    )
    
    # Metadata arguments
    parser.add_argument(
        "--sample-id",
        type=str,
        required=True,
        help="Unique sample identifier"
    )
    parser.add_argument(
        "--batch",
        type=str,
        default=None,
        help="Batch identifier (for integration)"
    )
    parser.add_argument(
        "--condition",
        type=str,
        default=None,
        help="Experimental condition (e.g., control, treated)"
    )
    parser.add_argument(
        "--metadata-json",
        type=str,
        default=None,
        help="Additional metadata as JSON string"
    )
    
    # Optional parameters
    parser.add_argument(
        "--rscript-path",
        type=str,
        default="Rscript",
        help="Path to Rscript executable (default: Rscript)"
    )
    parser.add_argument(
        "--make-unique",
        action="store_true",
        help="Make gene names unique (append suffix to duplicates)"
    )
    
    return parser.parse_args()


def log_versions():
    """Log version information for reproducibility."""
    logger.info("=== Software Versions ===")
    logger.info(f"Python: {sys.version.split()[0]}")
    logger.info(f"scanpy: {sc.__version__}")
    logger.info(f"anndata: {ad.__version__}")
    logger.info(f"numpy: {np.__version__}")
    logger.info(f"pandas: {pd.__version__}")
    logger.info("========================")


def validate_input_path(input_path: str, input_type: str) -> Path:
    """Validate that input path exists and is correct type.
    
    Args:
        input_path: Path to input file or directory
        input_type: Type of input (h5ad, rds, mtx)
        
    Returns:
        Path: Validated Path object
        
    Raises:
        FileNotFoundError: If input does not exist
        ValueError: If input type doesn't match path type
    """
    path = Path(input_path)
    
    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")
    
    if input_type == "mtx":
        if not path.is_dir():
            raise ValueError(f"MTX input must be a directory, got file: {input_path}")
        # Check for required files
        required_files = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
        missing = [f for f in required_files if not (path / f).exists()]
        if missing:
            # Try without .gz extension
            required_files_alt = ["matrix.mtx", "features.tsv", "barcodes.tsv"]
            missing_alt = [f for f in required_files_alt if not (path / f).exists()]
            if missing_alt:
                raise ValueError(f"MTX directory missing required files: {missing}")
    else:
        if not path.is_file():
            raise ValueError(f"{input_type.upper()} input must be a file, got directory: {input_path}")
    
    return path


def load_h5ad(input_path: Path) -> ad.AnnData:
    """Load H5AD file directly.
    
    Args:
        input_path: Path to H5AD file
        
    Returns:
        AnnData: Loaded AnnData object
    """
    logger.info(f"Loading H5AD file: {input_path}")
    try:
        adata = ad.read_h5ad(input_path)
        logger.info(f"Successfully loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.error(f"Failed to load H5AD file: {e}")
        raise


def load_rds(input_path: Path, rscript_path: str = "Rscript") -> ad.AnnData:
    """Load RDS (Seurat) file via direct CSV export.
    
    This method converts RDS to CSV/MTX format using R, then reads it with scanpy.
    This approach is more compatible with conda (no SeuratDisk dependency).
    
    Args:
        input_path: Path to RDS file
        rscript_path: Path to Rscript executable
        
    Returns:
        AnnData: Converted AnnData object
        
    Raises:
        RuntimeError: If R conversion fails
    """
    logger.info(f"Loading RDS file via direct conversion: {input_path}")
    
    # Create temporary directory for CSV export
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_dir = Path(tmpdir)
        
        # Step 1: Convert RDS to CSV/MTX using R script
        logger.info("Step 1/2: Exporting Seurat object to CSV/MTX format...")
        
        # Get path to R conversion script
        script_dir = Path(__file__).parent
        r_script = script_dir / "convert_rds_to_h5ad_direct.R"
        
        if not r_script.exists():
            raise FileNotFoundError(
                f"R conversion script not found: {r_script}\n"
                "Expected location: bin/convert_rds_to_h5ad_direct.R"
            )
        
        # Execute R script
        cmd = [rscript_path, str(r_script), str(input_path), str(csv_dir)]
        logger.debug(f"Executing: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        # Log R script output
        if result.stdout:
            for line in result.stdout.splitlines():
                logger.info(f"[R] {line}")
        
        if result.returncode != 0:
            error_msg = f"R conversion failed with exit code {result.returncode}"
            if result.stderr:
                error_msg += f"\nStderr: {result.stderr}"
            raise RuntimeError(error_msg)
        
        # Step 2: Read CSV/MTX data with scanpy
        logger.info("Step 2/2: Reading exported data with scanpy...")
        
        try:
            # Read count matrix
            from scipy.io import mmread
            counts = mmread(csv_dir / "matrix.mtx").T.tocsr()
            
            # Read gene names
            genes = pd.read_csv(csv_dir / "genes.tsv", sep="\t", header=None, names=["gene"])
            
            # Read cell barcodes
            barcodes = pd.read_csv(csv_dir / "barcodes.tsv", sep="\t", header=None, names=["barcode"])
            
            # Read metadata
            metadata = pd.read_csv(csv_dir / "metadata.csv", index_col=0)
            
            # Create AnnData object
            adata = ad.AnnData(
                X=counts,
                obs=metadata,
                var=genes.set_index("gene")
            )
            
            # Set obs_names and var_names
            adata.obs_names = barcodes["barcode"].values
            
            logger.info(f"Successfully converted: {adata.n_obs} cells × {adata.n_vars} genes")
            return adata
            
        except Exception as e:
            logger.error(f"Failed to read exported CSV data: {e}")
            raise RuntimeError(
                f"Failed to convert CSV/MTX to AnnData: {e}\n"
                "Check that the R script exported data correctly."
            )


def load_mtx(input_path: Path, make_unique: bool = False) -> ad.AnnData:
    """Load 10X MTX format and convert to AnnData.
    
    Args:
        input_path: Path to directory containing MTX files
        make_unique: Make gene names unique if duplicates exist
        
    Returns:
        AnnData: Loaded AnnData object
    """
    logger.info(f"Loading 10X MTX directory: {input_path}")
    try:
        adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', make_unique=make_unique)
        logger.info(f"Successfully loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.error(f"Failed to load MTX directory: {e}")
        raise


def integrate_metadata(
    adata: ad.AnnData,
    sample_id: str,
    batch: Optional[str] = None,
    condition: Optional[str] = None,
    additional_metadata: Optional[Dict] = None
) -> ad.AnnData:
    """Integrate sample metadata into AnnData.obs following scAnnex conventions.
    
    Per InitProject.md Section 9, .obs must include:
    - sample_id: Unique sample identifier
    - batch: Batch identifier (for integration)
    - condition: Experimental condition
    
    Args:
        adata: Input AnnData object
        sample_id: Sample identifier
        batch: Batch identifier
        condition: Experimental condition
        additional_metadata: Additional metadata dictionary
        
    Returns:
        AnnData: Modified AnnData with metadata
    """
    logger.info("Integrating sample metadata into .obs...")
    
    # Add required metadata columns
    adata.obs['sample_id'] = sample_id
    
    if batch is not None:
        adata.obs['batch'] = batch
    else:
        adata.obs['batch'] = 'batch_unknown'
        logger.warning("No batch specified, using 'batch_unknown'")
    
    if condition is not None:
        adata.obs['condition'] = condition
    else:
        adata.obs['condition'] = 'condition_unknown'
        logger.warning("No condition specified, using 'condition_unknown'")
    
    # Add additional metadata if provided
    if additional_metadata:
        for key, value in additional_metadata.items():
            if key not in ['sample_id', 'batch', 'condition']:
                adata.obs[key] = value
                logger.info(f"Added metadata column: {key}")
    
    logger.info(f"Metadata columns in .obs: {list(adata.obs.columns)}")
    return adata


def standardize_anndata(
    adata: ad.AnnData,
    input_type: str,
    input_path: str,
    sample_id: str
) -> ad.AnnData:
    """Standardize AnnData structure following InitProject.md Section 9.
    
    Conventions:
    - .obs: Cell-level metadata (sample_id, batch, condition, QC metrics)
    - .var: Gene-level metadata (gene names as index)
    - .uns: Unstructured metadata (parameters, conversion info)
    - .layers: Alternative representations ('counts' for raw)
    - .obsm: Multi-dimensional annotations (empty initially, for PCA/UMAP)
    
    Args:
        adata: Input AnnData object
        input_type: Type of input file
        input_path: Original input path
        sample_id: Sample identifier
        
    Returns:
        AnnData: Standardized AnnData object
    """
    logger.info("Standardizing AnnData structure...")
    
    # Ensure gene names are in index
    if adata.var.index.name is None:
        adata.var.index.name = 'gene_id'
    
    # Ensure cell barcodes have a name
    if adata.obs.index.name is None:
        adata.obs.index.name = 'cell_id'
    
    # Preserve raw counts in .layers['counts'] if not already present
    if 'counts' not in adata.layers:
        logger.info("Preserving raw counts in .layers['counts']")
        adata.layers['counts'] = adata.X.copy()
    
    # Initialize .obsm if empty (for future PCA/UMAP)
    if len(adata.obsm.keys()) == 0:
        logger.debug(".obsm is empty (will be populated during dimensionality reduction)")
    
    # Store conversion metadata in .uns
    adata.uns['scannex'] = {
        'conversion': {
            'input_type': input_type,
            'input_path': str(input_path),
            'sample_id': sample_id,
            'timestamp': datetime.now().isoformat(),
            'tool': 'scAnnex/unify_input',
            'version': '1.0.0'
        }
    }
    
    # Add software versions
    adata.uns['scannex']['versions'] = {
        'python': sys.version.split()[0],
        'scanpy': sc.__version__,
        'anndata': ad.__version__,
        'numpy': np.__version__,
        'pandas': pd.__version__
    }
    
    logger.info("AnnData structure standardized:")
    logger.info(f"  .obs: {adata.obs.shape[1]} columns")
    logger.info(f"  .var: {adata.var.shape[1]} columns")
    logger.info(f"  .uns: {list(adata.uns.keys())}")
    logger.info(f"  .layers: {list(adata.layers.keys())}")
    logger.info(f"  .obsm: {list(adata.obsm.keys())}")
    
    return adata


def validate_anndata(adata: ad.AnnData) -> None:
    """Validate AnnData object meets scAnnex requirements.
    
    Args:
        adata: AnnData object to validate
        
    Raises:
        ValueError: If validation fails
    """
    logger.info("Validating AnnData structure...")
    
    # Check required .obs columns
    required_obs = ['sample_id', 'batch', 'condition']
    missing_obs = [col for col in required_obs if col not in adata.obs.columns]
    if missing_obs:
        raise ValueError(f"Missing required .obs columns: {missing_obs}")
    
    # Check for raw counts in .layers
    if 'counts' not in adata.layers:
        logger.warning("No 'counts' layer found - raw counts may not be preserved")
    
    # Check for valid dimensions
    if adata.n_obs == 0:
        raise ValueError("AnnData object has 0 cells")
    if adata.n_vars == 0:
        raise ValueError("AnnData object has 0 genes")
    
    # Check for gene name duplicates
    if adata.var.index.duplicated().any():
        n_duplicates = adata.var.index.duplicated().sum()
        logger.warning(
            f"Found {n_duplicates} duplicate gene names. "
            "Consider using --make-unique flag."
        )
    
    logger.info("✓ AnnData validation passed")


def main():
    """Main execution function."""
    args = parse_args()
    
    # Log versions
    log_versions()
    
    try:
        # Validate input path
        input_path = validate_input_path(args.input, args.input_type)
        
        # Load data based on input type
        logger.info(f"Processing {args.input_type.upper()} input...")
        
        if args.input_type == "h5ad":
            adata = load_h5ad(input_path)
        elif args.input_type == "rds":
            adata = load_rds(input_path, args.rscript_path)
        elif args.input_type == "mtx":
            adata = load_mtx(input_path, args.make_unique)
        else:
            raise ValueError(f"Unknown input type: {args.input_type}")
        
        # Parse additional metadata if provided
        additional_metadata = None
        if args.metadata_json:
            try:
                additional_metadata = json.loads(args.metadata_json)
            except json.JSONDecodeError as e:
                logger.error(f"Failed to parse metadata JSON: {e}")
                raise
        
        # Integrate metadata
        adata = integrate_metadata(
            adata,
            sample_id=args.sample_id,
            batch=args.batch,
            condition=args.condition,
            additional_metadata=additional_metadata
        )
        
        # Standardize AnnData structure
        adata = standardize_anndata(
            adata,
            input_type=args.input_type,
            input_path=args.input,
            sample_id=args.sample_id
        )
        
        # Validate output
        validate_anndata(adata)
        
        # Save unified H5AD
        logger.info(f"Saving to: {args.output}")
        adata.write_h5ad(args.output, compression='gzip')
        
        # Log final statistics
        output_size = Path(args.output).stat().st_size / (1024 ** 2)
        logger.info(f"✓ Successfully created unified H5AD")
        logger.info(f"  Output file: {args.output}")
        logger.info(f"  File size: {output_size:.2f} MB")
        logger.info(f"  Dimensions: {adata.n_obs} cells × {adata.n_vars} genes")
        logger.info(f"  Sample ID: {args.sample_id}")
        logger.info("✓ Input unification complete")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
