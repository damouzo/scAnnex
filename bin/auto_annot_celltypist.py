#!/usr/bin/env python3

"""
=============================================================================
CellTypist Automated Cell-Type Annotation
=============================================================================
Uses CellTypist models to annotate cell types based on gene expression

Features:
- Pre-trained model support (Immune_All_Low, Immune_All_High, etc.)
- Custom model loading
- Majority voting option for stable predictions
- Standardized CSV output for integration

Requirements:
- celltypist >= 1.6.0
- scanpy >= 1.9.0
- pandas >= 1.5.0

Usage:
    python auto_annot_celltypist.py \\
        --input integrated.h5ad \\
        --output celltypist_annotations.csv \\
        --model Immune_All_Low.pkl \\
        --majority-voting
=============================================================================
"""

import argparse
import sys
import warnings
from pathlib import Path

import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import anndata as ad

# Enable writing of nullable strings (required for anndata >= 0.11)
# Only set if available (older versions don't have this attribute)
if hasattr(ad, 'settings') and hasattr(ad.settings, 'allow_write_nullable_strings'):
    ad.settings.allow_write_nullable_strings = True

warnings.filterwarnings('ignore')


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="CellTypist automated cell-type annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use pre-trained Immune model
  %(prog)s --input integrated.h5ad --output annotations.csv \\
    --model Immune_All_Low.pkl

  # With majority voting for stability
  %(prog)s --input integrated.h5ad --output annotations.csv \\
    --model Immune_All_High.pkl --majority-voting

  # Download and use specific model
  %(prog)s --input integrated.h5ad --output annotations.csv \\
    --model Healthy_Adult_Heart.pkl --download-model
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input H5AD file (integrated data)"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV file (cell_id, label, score, tool)"
    )
    parser.add_argument(
        "--output-h5ad",
        type=str,
        default=None,
        help="Output H5AD file with annotations added to .obs (optional)"
    )
    
    # Model selection
    parser.add_argument(
        "--model",
        type=str,
        default="Immune_All_Low.pkl",
        help="CellTypist model name or path (default: Immune_All_Low.pkl)"
    )
    parser.add_argument(
        "--download-model",
        action="store_true",
        help="Download model if not available locally"
    )
    parser.add_argument(
        "--list-models",
        action="store_true",
        help="List all available CellTypist models and exit"
    )
    
    # Annotation parameters
    parser.add_argument(
        "--majority-voting",
        action="store_true",
        help="Enable majority voting for stable predictions"
    )
    parser.add_argument(
        "--over-clustering-resolution",
        type=float,
        default=0.5,
        help="Resolution for over-clustering (used with majority voting, default: 0.5)"
    )
    parser.add_argument(
        "--min-prop",
        type=float,
        default=0.5,
        help="Minimum proportion for majority voting (default: 0.5)"
    )
    
    # Processing options
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="Use GPU acceleration (if available)"
    )
    parser.add_argument(
        "--backed",
        action="store_true",
        help="Use backed mode for large datasets"
    )
    
    return parser.parse_args()


def list_available_models():
    """List all available CellTypist models."""
    print("\n" + "="*60)
    print("Available CellTypist Models")
    print("="*60 + "\n")
    
    try:
        available_models = models.models_description()
        for model_name, description in available_models.items():
            print(f"{model_name}")
            print(f"  Description: {description}")
            print()
    except Exception as e:
        print(f"Error retrieving models: {e}")
        sys.exit(1)


def load_data(h5ad_path, backed=False):
    """Load H5AD file."""
    print(f"Loading H5AD file: {h5ad_path}")
    
    if backed:
        adata = sc.read_h5ad(h5ad_path, backed='r')
    else:
        adata = sc.read_h5ad(h5ad_path)
    
    print(f"  Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"  Observations: {list(adata.obs.columns)}")
    
    # CellTypist expects log1p normalized data in .X
    # Check if normalized data is stored in layers (Scanpy best practice)
    if 'log1p_norm' in adata.layers:
        print(f"  Using log1p normalized data from .layers['log1p_norm'] for CellTypist")
        adata.X = adata.layers['log1p_norm'].copy()
    elif 'normalized' in adata.layers:
        print(f"  Using normalized data from .layers['normalized'] for CellTypist")
        adata.X = adata.layers['normalized'].copy()
    else:
        print(f"  WARNING: No log1p normalized data found in layers, using .X as-is")
    
    return adata


def load_model(model_name, download=False):
    """Load or download CellTypist model."""
    print(f"\nLoading CellTypist model: {model_name}")
    
    try:
        # Check if model exists locally
        if Path(model_name).exists():
            print(f"  Loading from local path: {model_name}")
            model = models.Model.load(model_name)
        else:
            # Try to load from CellTypist repository
            if download:
                print(f"  Downloading model from CellTypist repository...")
                models.download_models(model=model_name.replace('.pkl', ''))
            
            model = models.Model.load(model_name)
        
        print(f"  Model loaded successfully")
        print(f"  Cell types: {len(model.cell_types)} types")
        
        return model
        
    except Exception as e:
        print(f"ERROR: Failed to load model: {e}")
        print(f"\nTip: Use --list-models to see available models")
        print(f"     Use --download-model to download the model")
        sys.exit(1)


def run_celltypist(adata, model, majority_voting=False, 
                   over_clustering_resolution=0.5, min_prop=0.5):
    """Run CellTypist annotation."""
    print("\nRunning CellTypist annotation...")
    
    try:
        # For majority voting, use existing leiden clustering or compute it
        over_clustering_key = None
        if majority_voting:
            # Check if leiden clustering exists at the specified resolution
            leiden_key = f'leiden_{over_clustering_resolution}'
            if leiden_key not in adata.obs.columns:
                # Compute leiden clustering if it doesn't exist
                print(f"  Computing leiden clustering at resolution {over_clustering_resolution} for majority voting...")
                import scanpy as sc
                sc.pp.neighbors(adata, use_rep='X_pca')
                sc.tl.leiden(adata, resolution=over_clustering_resolution, key_added=leiden_key)
            over_clustering_key = leiden_key
            print(f"  Using {leiden_key} for majority voting")
        
        # Run prediction
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=majority_voting,
            over_clustering=over_clustering_key,
            min_prop=min_prop if majority_voting else 0
        )
        
        print(f"  Annotation completed")
        print(f"  Majority voting: {'enabled' if majority_voting else 'disabled'}")
        
        # Get prediction results
        if majority_voting:
            labels = predictions.predicted_labels.majority_voting
            print(f"  Using majority voting labels")
        else:
            labels = predictions.predicted_labels.predicted_labels
            print(f"  Using direct predictions")
        
        # Get confidence scores
        conf_scores = predictions.probability_matrix.max(axis=1)
        
        print(f"  Unique cell types: {labels.nunique()}")
        print(f"  Mean confidence score: {conf_scores.mean():.3f}")
        
        return labels, conf_scores
        
    except Exception as e:
        print(f"ERROR: CellTypist annotation failed: {e}")
        sys.exit(1)


def save_results(cell_ids, labels, scores, output_path):
    """Save annotations to standardized CSV format."""
    print(f"\nSaving results to: {output_path}")
    
    # Create DataFrame
    results_df = pd.DataFrame({
        'cell_id': cell_ids,
        'label': labels.values,
        'score': scores.values,
        'tool': 'celltypist'
    })
    
    # Save to CSV
    results_df.to_csv(output_path, index=False)
    
    print(f"  Saved {len(results_df)} annotations")
    print(f"\nCell type distribution:")
    print(results_df['label'].value_counts().head(10))
    
    return results_df


def save_h5ad_with_annotations(adata, labels, scores, output_h5ad_path):
    """Save H5AD file with CellTypist annotations added to .obs."""
    print(f"\nAdding CellTypist annotations to H5AD...")
    
    # Add annotations to .obs
    adata.obs['auto_annot_celltypist'] = labels.values
    adata.obs['auto_annot_celltypist_score'] = scores.values
    
    print(f"  Added columns to .obs:")
    print(f"    - auto_annot_celltypist: cell type labels")
    print(f"    - auto_annot_celltypist_score: confidence scores")
    
    # Save annotated H5AD
    print(f"\nSaving annotated H5AD to: {output_h5ad_path}")
    adata.write_h5ad(output_h5ad_path)
    
    print(f"  Successfully saved annotated H5AD")


def main():
    """Main execution function."""
    args = parse_args()
    
    # List models and exit if requested
    if args.list_models:
        list_available_models()
        sys.exit(0)
    
    print("="*60)
    print("CellTypist Automated Annotation")
    print("="*60)
    
    # Step 1: Load data
    # Force backed=False if we need to save to H5AD (backed mode is read-only)
    use_backed = args.backed and not args.output_h5ad
    if args.output_h5ad and args.backed:
        print("  Note: Disabling backed mode because --output-h5ad requires write access")
    
    adata = load_data(args.input, backed=use_backed)
    
    # Step 2: Load model
    model = load_model(args.model, download=args.download_model)
    
    # Step 3: Run CellTypist
    labels, scores = run_celltypist(
        adata,
        model,
        majority_voting=args.majority_voting,
        over_clustering_resolution=args.over_clustering_resolution,
        min_prop=args.min_prop
    )
    
    # Step 4: Save results
    results_df = save_results(
        cell_ids=adata.obs_names,
        labels=labels,
        scores=scores,
        output_path=args.output
    )
    
    # Step 5: Save annotated H5AD if requested
    if args.output_h5ad:
        save_h5ad_with_annotations(
            adata=adata,
            labels=labels,
            scores=scores,
            output_h5ad_path=args.output_h5ad
        )
    
    print("\n" + "="*60)
    print("CellTypist annotation completed successfully!")
    print("="*60)


if __name__ == "__main__":
    main()
