#!/usr/bin/env python3
"""
Build pseudo-bulk count matrices per sample x cell type from raw counts.

Aggregates raw, un-normalized counts by summing over cells within each
(sample, cell_type) pair. Output values are integer and are NOT normalized
or log-transformed: DESeq2 performs its own median-of-ratios normalization.

This avoids pseudo-replication by using the sample x cell type aggregate as
the unit of biological replication (never individual cells).

Usage:
    python pseudobulk_build.py \\
        --input auto_annotated_global.h5ad \\
        --counts-out pseudobulk_counts.tsv \\
        --metadata-out pseudobulk_metadata.csv \\
        --summary-out pseudobulk_build_summary.json
"""

import argparse
import json
import logging
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build pseudo-bulk count matrices per sample x cell type"
    )
    parser.add_argument(
        "--input", required=True, help="Auto-annotated H5AD (needs .layers['counts'])"
    )
    parser.add_argument(
        "--counts-out",
        required=True,
        help="Output genes x aggregates TSV (integer counts)",
    )
    parser.add_argument(
        "--metadata-out", required=True, help="Output aggregate metadata CSV"
    )
    parser.add_argument(
        "--summary-out", required=True, help="Output build summary JSON"
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=10,
        help="Minimum number of cells per aggregate; smaller aggregates are dropped",
    )
    parser.add_argument(
        "--groupby",
        default="cell_type",
        help="Annotation column to group aggregates by",
    )
    parser.add_argument(
        "--sample-col",
        default="sample_id",
        help="Column identifying the biological sample",
    )
    parser.add_argument(
        "--condition-col",
        default="condition",
        help="Column identifying the experimental condition",
    )
    return parser.parse_args()


def require_counts_layer(adata: ad.AnnData) -> np.ndarray:
    """Return center count matrix (integer) preserving gene order."""
    if "counts" not in adata.layers:
        raise RuntimeError(
            "No '.layers['counts']' found in input. Cannot build pseudo-bulk from raw counts."
        )
    matrix = adata.layers["counts"]
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix)


def build_pseudobulk(
    adata,
    min_cells=10,
    groupby="cell_type",
    sample_col="sample_id",
    condition_col="condition",
):
    """Aggregate raw counts per sample x cell type.

    Returns (counts_df, meta_df, summary). Values are integer raw sums and are
    never normalized/log-transformed (DESeq2 performs its own normalization).
    """
    for col in (sample_col, condition_col, groupby):
        if col not in adata.obs.columns:
            raise ValueError(f"Column '{col}' not found in adata.obs")

    counts = require_counts_layer(adata)

    if not np.isfinite(counts).all():
        raise ValueError(
            "Counts matrix contains non-finite values; cannot build pseudo-bulk"
        )
    if not np.all(counts == np.round(counts)):
        raise ValueError(
            "Counts matrix is not integer-valued. Pseudo-bulk DESeq2 requires raw integer counts; "
            "check that '.layers['counts']' holds non-normalized data."
        )

    sample_ids = adata.obs[sample_col].astype(str).values
    conditions = adata.obs[condition_col].astype(str).values
    cell_types = adata.obs[groupby].astype(str).values

    valid = (
        ~pd.isna(adata.obs[sample_col].values)
        & ~pd.isna(adata.obs[condition_col].values)
        & ~pd.isna(adata.obs[groupby].values)
        & (adata.obs[groupby].astype(str).str.strip() != "")
    )
    n_dropped = int((~valid).sum())
    if n_dropped:
        logger.warning(
            f"Dropping {n_dropped} cells with missing sample/condition/cell_type"
        )

    gene_names = adata.var_names.astype(str).values
    aggregate_keys: List[Tuple[str, str]] = []
    aggregate_counts: List[np.ndarray] = []
    aggregate_cell_count: Dict[Tuple, int] = {}
    aggregate_condition: Dict[Tuple, str] = {}

    keys = set(zip(sample_ids[valid], cell_types[valid]))
    for sample_id, cell_type in sorted(keys):
        mask = valid & (sample_ids == sample_id) & (cell_types == cell_type)
        n_cells = int(mask.sum())
        summed = np.asarray(counts[mask, :].sum(axis=0)).astype(np.int64).ravel()
        agg_key = (sample_id, cell_type)
        aggregate_cell_count[agg_key] = n_cells
        # condition from first (all share it)
        aggregate_condition[agg_key] = conditions[np.where(mask)[0][0]]
        if n_cells >= min_cells:
            aggregate_counts.append(summed)
            aggregate_keys.append(agg_key)

    # Replicate group size = number of distinct samples per condition within each cell type
    replicate_count: Dict[Tuple, int] = {}
    for sample_id, cell_type in aggregate_keys:
        cond = aggregate_condition[(sample_id, cell_type)]
        n_samples = len(
            {
                s
                for (s, ct) in aggregate_keys
                if ct == cell_type and aggregate_condition[(s, ct)] == cond
            }
        )
        replicate_count[(sample_id, cell_type)] = n_samples

    if not aggregate_keys:
        raise ValueError(
            "No aggregates passed the minimum cell filter. No pseudo-bulk could be built."
        )

    counts_matrix = np.vstack(aggregate_counts)
    agg_labels = [
        f"{sample_id}::{cell_type}" for (sample_id, cell_type) in aggregate_keys
    ]

    counts_df = pd.DataFrame(counts_matrix.T, index=gene_names, columns=agg_labels)
    counts_df = counts_df[(counts_df > 0).any(axis=1)]
    counts_df.index.name = "gene"

    meta_rows = []
    for sample_id, cell_type in aggregate_keys:
        cond = aggregate_condition[(sample_id, cell_type)]
        raw = aggregate_cell_count[(sample_id, cell_type)]
        meta_rows.append(
            {
                "sample_id": sample_id,
                "cell_type": cell_type,
                "condition": cond,
                "n_cells_raw": raw,
                "cells_after_min_cells": raw,
                "excluded": False,
                "replicate_group_size_condition": replicate_count[
                    (sample_id, cell_type)
                ],
            }
        )
    meta_df = pd.DataFrame(meta_rows)

    summary = {
        "success": True,
        "n_cells_input": int(adata.n_obs),
        "n_genes": int(counts_df.shape[0]),
        "n_aggregates": len(aggregate_keys),
        "min_cells": min_cells,
        "groupby": groupby,
        "n_cells_dropped_no_metadata": n_dropped,
        "n_cell_types": int(meta_df["cell_type"].nunique()),
        "n_samples": int(meta_df["sample_id"].nunique()),
        "n_conditions": int(meta_df["condition"].nunique()),
        "cell_types": sorted(meta_df["cell_type"].unique().tolist()),
    }
    return counts_df, meta_df, summary


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    counts_df, meta_df, summary = build_pseudobulk(
        adata,
        min_cells=args.min_cells,
        groupby=args.groupby,
        sample_col=args.sample_col,
        condition_col=args.condition_col,
    )
    counts_df.to_csv(args.counts_out, sep="\t")
    meta_df.to_csv(args.metadata_out, index=False)
    with open(args.summary_out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, default=str)

    logger.info(
        f"Built {summary['n_aggregates']} aggregates from {summary['n_samples']} samples "
        f"across {summary['n_cell_types']} cell types"
    )
    logger.info(f"  Counts: {args.counts_out}")
    logger.info(f"  Metadata: {args.metadata_out}")


if __name__ == "__main__":
    main()
