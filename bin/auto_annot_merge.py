#!/usr/bin/env python3

import argparse
import json
import re

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
    ad.settings.allow_write_nullable_strings = True


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge automatic annotation CSV outputs into one H5AD"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--celltypist", required=True)
    parser.add_argument("--celltypist-status", required=True)
    parser.add_argument("--sctype", required=True)
    parser.add_argument("--sctype-status", required=True)
    parser.add_argument("--azimuth", required=True)
    parser.add_argument("--azimuth-status", required=True)
    parser.add_argument("--singler", required=True)
    parser.add_argument("--singler-status", required=True)
    parser.add_argument(
        "--cell-type-column",
        default=None,
        help="Canonical cell type column. Default: auto-detect SingleR pruned column",
    )
    return parser.parse_args()


TOOL_PREFIX = {
    "celltypist": "auto_annot_celltypist_",
    "sctype": "auto_annot_sctype",
    "azimuth": "auto_annot_azimuth_",
    "singler": "auto_annot_singler_",
}

SUFFIXES = ("_score", "_delta_next", "_pruned", "_l1", "_l2", "_cluster")


def _label_columns(cols, prefix):
    """Return the main label columns for a tool prefix (excluding score/sibling columns)."""
    found = []
    for col in cols:
        if col == prefix:  # sctype main column has no trailing suffix
            found.append(col)
        elif col.startswith(prefix) and not col.endswith(SUFFIXES):
            found.append(col)
    return found


def _score_column_for(label_col):
    if label_col == "auto_annot_sctype":
        return "auto_annot_sctype_score"
    return f"{label_col}_score"


def _normalize_label(val):
    """Normalize a cell type label for comparison (lowercase, strip non-alphanumerics)."""
    if val is pd.NA or val is None:
        return None
    if isinstance(val, float) and np.isnan(val):
        return None
    s = str(val).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s if s else None


def read_status(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except Exception:
        return {"success": False, "message": f"Failed to parse status file: {path}"}


def read_annotations(path):
    df = pd.read_csv(path)
    if "cell_id" not in df.columns:
        return pd.DataFrame()
    valid_mask = df["cell_id"].notna()
    df = df.loc[valid_mask].copy()
    df["cell_id"] = df["cell_id"].astype(str).str.strip()
    df = df[df["cell_id"] != ""]
    return df.set_index("cell_id")


def merge_df(adata, csv_path):
    incoming = read_annotations(csv_path)
    if incoming.empty:
        return {"columns": [], "n_rows": 0, "n_common_cells": 0}

    obs_names = pd.Index(adata.obs_names.astype(str))
    common_cells = obs_names.intersection(incoming.index)
    if len(common_cells) == 0:
        return {"columns": [], "n_rows": int(incoming.shape[0]), "n_common_cells": 0}

    cols_added = []
    for col in incoming.columns:
        source = incoming.loc[common_cells, col]

        if pd.api.types.is_numeric_dtype(source):
            target = pd.Series(np.nan, index=adata.obs_names, dtype="float64")
            target.loc[common_cells] = pd.to_numeric(source, errors="coerce").values
        else:
            target = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
            target.loc[common_cells] = source.astype("string").values

        adata.obs[col] = target
        cols_added.append(col)
    return {
        "columns": cols_added,
        "n_rows": int(incoming.shape[0]),
        "n_common_cells": int(len(common_cells)),
    }


def resolve_cell_type_tool(config):
    """Pick the SingleR label column (prefer '_pruned') used as canonical source."""
    cols = list(config.get("singler", {}).get("columns", []) or [])
    label_cols = _label_columns(cols, TOOL_PREFIX["singler"])
    if not label_cols:
        return None, None
    # Prefer the pruned label column (pruning applied) for a stable, QC-filtered label
    for label_col in label_cols:
        pruned_col = f"{label_col}_pruned"
        if pruned_col in cols:
            return label_col, pruned_col
    return label_cols[0], label_cols[0]


def add_canonical_cell_type(adata, merge_report, configured):
    """Add canonical cell_type, cell_type_agreement and cell_type_confidence columns.

    cell_type default = SingleR pruned label (column stable regardless of reference).
    cell_type_agreement = number of tools (0-4) with a non-NA label equal to cell_type.
    cell_type_confidence = SingleR score (for the default source) or NA otherwise.
    """
    single_label_col, single_pruned_col = resolve_cell_type_tool(merge_report)

    if configured:
        cell_type_col = configured
        if cell_type_col not in adata.obs.columns:
            raise ValueError(
                f"--cell-type-column '{cell_type_col}' not found in merged annotations. "
                f"Available: {list(adata.obs.columns)}"
            )
        confidence_col = (
            _score_column_for(cell_type_col)
            if cell_type_col.startswith(TOOL_PREFIX["singler"])
            else None
        )
        source_tool = next(
            (
                t
                for t, p in TOOL_PREFIX.items()
                if cell_type_col == p or cell_type_col.startswith(p)
            ),
            None,
        )
    elif single_label_col is not None:
        cell_type_col = single_pruned_col
        confidence_col = _score_column_for(single_label_col)
        source_tool = "singler"
    else:
        cell_type_col = None
        confidence_col = None
        source_tool = None

    if cell_type_col is not None:
        adata.obs["cell_type"] = adata.obs[cell_type_col]

    # Per-tool main label columns
    tool_labels = {}
    for tool, prefix in TOOL_PREFIX.items():
        if tool == "singler" and source_tool == "singler":
            # Compare against the exact canonical value (same as default source)
            ref_col = cell_type_col if cell_type_col else (single_label_col or "")
        else:
            cols = merge_report.get(tool, {}).get("columns", []) or []
            labels = _label_columns(cols, prefix)
            ref_col = labels[0] if labels else None
        if ref_col and ref_col in adata.obs.columns:
            tool_labels[tool] = ref_col

    canonical = (
        adata.obs["cell_type"]
        if cell_type_col
        else pd.Series(pd.NA, index=adata.obs.index)
    )

    normalized_canonical = canonical.map(_normalize_label)
    agreement = pd.Series(0, index=adata.obs.index, dtype=int)
    annotated_count = pd.Series(0, index=adata.obs.index, dtype=int)
    for tool, col in tool_labels.items():
        has_label = adata.obs[col].notna()
        normalized = adata.obs[col].map(_normalize_label)
        matches = (
            has_label
            & (normalized == normalized_canonical)
            & normalized_canonical.notna()
        )
        agreement = agreement + matches.astype(int)
        annotated_count = annotated_count + has_label.astype(int)

    adata.obs["cell_type_agreement"] = agreement
    adata.obs["n_tools_annotated"] = annotated_count

    if confidence_col and confidence_col in adata.obs.columns:
        adata.obs["cell_type_confidence"] = adata.obs[confidence_col]

    info = {
        "cell_type_column": cell_type_col,
        "cell_type_source_tool": source_tool,
        "cell_type_confidence_column": confidence_col,
        "tools_compared": {t: c for t, c in tool_labels.items()},
    }
    if cell_type_col is not None:
        annotated = adata.obs["cell_type"].notna().sum()
        total = adata.n_obs
        mean_agreement = (
            agreement[adata.obs["cell_type"].notna()].mean() if annotated else 0.0
        )
        info["percent_annotated"] = round(annotated / total * 100, 2) if total else 0.0
        info["percent_agreement"] = round(
            float(mean_agreement if mean_agreement == mean_agreement else 0.0), 4
        )
        info["n_cells_no_cell_type"] = int(total - annotated)
    else:
        info["percent_annotated"] = 0.0
        info["percent_agreement"] = 0.0
        info["warning"] = "No SingleR annotation found; cell_type not created"
    return info


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)

    status = {
        "celltypist": read_status(args.celltypist_status),
        "sctype": read_status(args.sctype_status),
        "azimuth": read_status(args.azimuth_status),
        "singler": read_status(args.singler_status),
    }

    merge_report = {
        "celltypist": merge_df(adata, args.celltypist),
        "sctype": merge_df(adata, args.sctype),
        "azimuth": merge_df(adata, args.azimuth),
        "singler": merge_df(adata, args.singler),
    }

    columns_added = {k: v["columns"] for k, v in merge_report.items()}

    cell_type_info = add_canonical_cell_type(adata, merge_report, args.cell_type_column)

    adata.uns["auto_annotation_summary"] = {
        "status": {k: bool(v.get("success", False)) for k, v in status.items()},
        "status_json": json.dumps(status, ensure_ascii=True),
        "columns_added": columns_added,
        "merge_report": merge_report,
        "cell_type": cell_type_info,
    }
    adata.write_h5ad(args.output)

    summary = {
        "success": True,
        "status": status,
        "columns_added": columns_added,
        "merge_report": merge_report,
        "cell_type": cell_type_info,
    }
    with open(args.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, default=str)


if __name__ == "__main__":
    main()
