#!/usr/bin/env python3

"""
=============================================================================
Contrasts CSV Validation for Differential Expression Analysis
=============================================================================
Validates contrasts.csv file against AnnData object before running DGE.

Features:
- Schema validation (required/optional columns)
- Column existence checks in adata.obs
- Group value validation
- Filter column and value validation
- Unique contrast_id enforcement

Usage:
    python validate_contrasts.py \\
        --contrasts contrasts.csv \\
        --h5ad integrated.h5ad \\
        --output validated_contrasts.csv
=============================================================================
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
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
        description="Validate contrasts.csv for differential expression analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate contrasts file
  %(prog)s --contrasts contrasts.csv --h5ad integrated.h5ad
  
  # Validate and save cleaned version
  %(prog)s --contrasts contrasts.csv --h5ad integrated.h5ad \\
    --output validated_contrasts.csv
        """
    )
    
    parser.add_argument(
        "--contrasts",
        type=str,
        required=True,
        help="Path to contrasts.csv file"
    )
    parser.add_argument(
        "--h5ad",
        type=str,
        required=True,
        help="Path to H5AD file for validation"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save validated contrasts.csv (optional)"
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Strict mode: exit on warnings (default: warnings allowed)"
    )
    
    return parser.parse_args()


def load_contrasts(contrasts_file: str) -> pd.DataFrame:
    """Load and perform basic validation of contrasts CSV."""
    logger.info(f"Loading contrasts file: {contrasts_file}")
    
    if not Path(contrasts_file).exists():
        logger.error(f"Contrasts file not found: {contrasts_file}")
        sys.exit(1)
    
    try:
        contrasts_df = pd.read_csv(contrasts_file)
    except Exception as e:
        logger.error(f"Failed to read contrasts file: {e}")
        sys.exit(1)
    
    logger.info(f"Loaded {len(contrasts_df)} contrasts")
    return contrasts_df


def validate_schema(contrasts_df: pd.DataFrame) -> Tuple[bool, List[str]]:
    """Validate contrasts DataFrame schema."""
    errors = []
    
    # Required columns
    required_columns = ['contrast_id', 'variable', 'group1', 'group2']
    for col in required_columns:
        if col not in contrasts_df.columns:
            errors.append(f"Missing required column: '{col}'")
    
    # Optional columns (if present, must be paired)
    if 'filter_column' in contrasts_df.columns or 'filter_value' in contrasts_df.columns:
        if 'filter_column' not in contrasts_df.columns:
            errors.append("'filter_value' present but 'filter_column' is missing")
        if 'filter_value' not in contrasts_df.columns:
            errors.append("'filter_column' present but 'filter_value' is missing")
    
    # Check for duplicate contrast_ids
    if 'contrast_id' in contrasts_df.columns:
        duplicates = contrasts_df['contrast_id'].duplicated()
        if duplicates.any():
            dup_ids = contrasts_df.loc[duplicates, 'contrast_id'].tolist()
            errors.append(f"Duplicate contrast_id values found: {dup_ids}")
    
    # Check for empty values in required columns
    for col in required_columns:
        if col in contrasts_df.columns:
            null_count = contrasts_df[col].isnull().sum()
            if null_count > 0:
                errors.append(f"Column '{col}' has {null_count} empty/null values")
    
    return len(errors) == 0, errors


def validate_against_adata(contrasts_df: pd.DataFrame, adata) -> Tuple[bool, List[str], List[str]]:
    """Validate contrasts against AnnData object."""
    errors = []
    warnings = []
    
    for idx, row in contrasts_df.iterrows():
        contrast_id = row['contrast_id']
        variable = row['variable']
        group1 = row['group1']
        group2 = row['group2']
        
        # Check if variable column exists in adata.obs
        if variable not in adata.obs.columns:
            errors.append(
                f"Contrast '{contrast_id}': Column '{variable}' not found in adata.obs. "
                f"Available columns: {', '.join(adata.obs.columns[:10])}..."
            )
            continue
        
        # Check if group1 and group2 exist in variable column
        unique_values = adata.obs[variable].unique()
        unique_values_str = [str(v) for v in unique_values]
        
        if str(group1) not in unique_values_str:
            errors.append(
                f"Contrast '{contrast_id}': group1 value '{group1}' not found in column '{variable}'. "
                f"Available values: {', '.join(unique_values_str)}"
            )
        
        if str(group2) not in unique_values_str:
            errors.append(
                f"Contrast '{contrast_id}': group2 value '{group2}' not found in column '{variable}'. "
                f"Available values: {', '.join(unique_values_str)}"
            )
        
        # Validate filter columns if present
        if 'filter_column' in row and pd.notna(row['filter_column']):
            filter_col = row['filter_column']
            filter_val = row['filter_value']
            
            if filter_col not in adata.obs.columns:
                errors.append(
                    f"Contrast '{contrast_id}': filter_column '{filter_col}' not found in adata.obs"
                )
            elif pd.notna(filter_val):
                filter_unique = adata.obs[filter_col].unique()
                filter_unique_str = [str(v) for v in filter_unique]
                
                if str(filter_val) not in filter_unique_str:
                    errors.append(
                        f"Contrast '{contrast_id}': filter_value '{filter_val}' not found in column '{filter_col}'. "
                        f"Available values: {', '.join(filter_unique_str)}"
                    )
        
        # Warning: Check cell counts (at least 30 cells per group recommended)
        if variable in adata.obs.columns:
            subset_adata = adata
            
            # Apply filter if specified
            if 'filter_column' in row and pd.notna(row['filter_column']) and pd.notna(row['filter_value']):
                filter_col = row['filter_column']
                filter_val = str(row['filter_value'])
                subset_adata = adata[adata.obs[filter_col].astype(str) == filter_val]
            
            # Count cells in each group
            if str(group1) in unique_values_str:
                n_group1 = (subset_adata.obs[variable].astype(str) == str(group1)).sum()
                if n_group1 < 30:
                    warnings.append(
                        f"Contrast '{contrast_id}': group1 '{group1}' has only {n_group1} cells "
                        f"(recommended: ≥30 for robust statistics)"
                    )
            
            if str(group2) in unique_values_str:
                n_group2 = (subset_adata.obs[variable].astype(str) == str(group2)).sum()
                if n_group2 < 30:
                    warnings.append(
                        f"Contrast '{contrast_id}': group2 '{group2}' has only {n_group2} cells "
                        f"(recommended: ≥30 for robust statistics)"
                    )
    
    return len(errors) == 0, errors, warnings


def main():
    args = parse_args()
    
    # Load contrasts file
    contrasts_df = load_contrasts(args.contrasts)
    
    # Validate schema
    logger.info("Validating contrasts schema...")
    schema_valid, schema_errors = validate_schema(contrasts_df)
    
    if not schema_valid:
        logger.error("Schema validation failed:")
        for error in schema_errors:
            logger.error(f"  - {error}")
        sys.exit(1)
    
    logger.info("✓ Schema validation passed")
    
    # Load H5AD file
    logger.info(f"Loading H5AD file: {args.h5ad}")
    try:
        adata = sc.read_h5ad(args.h5ad)
        logger.info(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    except Exception as e:
        logger.error(f"Failed to load H5AD file: {e}")
        sys.exit(1)
    
    # Validate against adata
    logger.info("Validating contrasts against AnnData object...")
    adata_valid, adata_errors, adata_warnings = validate_against_adata(contrasts_df, adata)
    
    if adata_errors:
        logger.error("Validation errors found:")
        for error in adata_errors:
            logger.error(f"  - {error}")
        sys.exit(1)
    
    if adata_warnings:
        logger.warning("Validation warnings:")
        for warning in adata_warnings:
            logger.warning(f"  - {warning}")
        
        if args.strict:
            logger.error("Strict mode enabled: exiting due to warnings")
            sys.exit(1)
    
    logger.info("✓ Contrasts validation passed")
    
    # Save validated contrasts if output specified
    if args.output:
        logger.info(f"Saving validated contrasts to: {args.output}")
        contrasts_df.to_csv(args.output, index=False)
        logger.info("✓ Validated contrasts saved")
    
    logger.info(f"✓ All {len(contrasts_df)} contrasts are valid and ready for DGE analysis")


if __name__ == "__main__":
    main()
