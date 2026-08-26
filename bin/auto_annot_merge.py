#!/usr/bin/env python3

import argparse
import json
import os
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
        help="Canonical cell type column. If set, forces a single source (no auto-detect).",
    )
    parser.add_argument(
        "--custom-cell-type-file",
        default=None,
        help=(
            "Optional CSV with a barcode/cell_id column and a 'cell_type' column. "
            "When provided, overrides the automatic source for the matching cells "
            "(e.g. exported from the dashboard Annotation Station)."
        ),
    )
    parser.add_argument(
        "--min-cell-type-annotated-pct",
        type=float,
        default=0.5,
        help=(
            "Minimum fraction of cells that must carry a cell_type for the run to "
            "succeed. If exceeded, the pipeline fails explicitly. Default: 0.5."
        ),
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


def _tool_for_column(col):
    """Return the tool name owning a merged annotation column."""
    if col is None:
        return None
    for tool, prefix in TOOL_PREFIX.items():
        if col == prefix or col.startswith(prefix):
            return tool
    return None


def load_custom_types(path):
    """Load a (barcode/cell_id, cell_type) CSV into a cell_id->str mapping.

    Only rows with a non-empty id are kept. An explicit error is raised when either
    required column is missing so a misnamed dashboard export fails loudly.
    """
    df = pd.read_csv(path)
    id_col = next((c for c in ("cell_id", "barcode", "cellbarcode", "barcode_id")
                   if c in df.columns), None)
    label_col = next((c for c in ("cell_type", "label", "annotation")
                      if c in df.columns), None)
    if id_col is None or label_col is None:
        raise ValueError(
            f"--custom-cell-type-file must contain an id column "
            f"(cell_id/barcode) and a 'cell_type' column. Found: {list(df.columns)}"
        )
    df = df.loc[df[id_col].notna()].copy()
    df[id_col] = df[id_col].astype(str).str.strip()
    df = df[df[id_col] != ""]
    df[label_col] = df[label_col].astype(str).str.strip()
    return df.set_index(id_col)[label_col]


def add_canonical_cell_type(adata, merge_report, configured, min_pct, custom_file=None):
    """Build the canonical ``cell_type`` column and agreement/confidence metadata.

    Resolution (no automatic tool switching):
      1. ``configured`` column, when set (any source).
      2. Default: SingleR pruned label, or its base label when pruned holds no data.
      ``custom_file``, when provided, overrides ``cell_type`` for the matching cells.

    If fewer than ``min_pct`` of cells end up annotated, the pipeline FAILS with an
    explicit, actionable message. It never silently picks a different annotator.
    """
    base_col = None
    source_tool = None

    if configured:
        if configured not in adata.obs.columns:
            raise ValueError(
                f"--cell-type-column '{configured}' not found in merged annotations. "
                f"Available: {list(adata.obs.columns)}"
            )
        base_col = configured
        source_tool = _tool_for_column(configured)
    else:
        cols = merge_report.get("singler", {}).get("columns", []) or []
        label_cols = _label_columns(cols, TOOL_PREFIX["singler"])
        if label_cols:
            base_col = label_cols[0]
            pruned = f"{base_col}_pruned"
            if pruned in cols and adata.obs[pruned].notna().any():
                base_col = pruned
            source_tool = "singler"

    if base_col is None:
        raise ValueError(
            "No default cell type source found (SingleR produced no label columns). "
            "Specify --cell-type-column or --custom-cell-type-file."
        )

    adata.obs["cell_type"] = adata.obs[base_col].astype("string").copy()

    n_custom = 0
    if custom_file:
        custom = load_custom_types(custom_file)
        if len(custom) > 0:
            inter = adata.obs.index.astype(str).intersection(custom.index)
            if len(inter) == 0:
                raise ValueError(
                    f"--custom-cell-type-file '{custom_file}' matched no cells. "
                    "Check that its id column matches the H5AD obs_names."
                )
            adata.obs.loc[inter, "cell_type"] = custom.loc[inter].values
            n_custom = int(len(inter))
        # An empty file (0 rows) means no override was supplied.

    annotated = int(adata.obs["cell_type"].notna().sum())
    total = adata.n_obs
    final_pct = annotated / total if total else 0.0

    if final_pct < min_pct:
        remedy = (
            f"Fix the '{base_col}' annotation ({source_tool}) or override with "
            f"--cell-type-column (another source) or --custom-cell-type-file."
        )
        if custom_file:
            remedy = (
                f"Custom file matched {n_custom} cells and the base source "
                f"'{base_col}' is missing labels for the rest. {remedy}"
            )
        raise ValueError(
            f"Only {round(final_pct * 100, 2)}% of {total} cells have a cell_type "
            f"(< min_cell_type_annotated_pct={min_pct}). {remedy}"
        )

    chose_custom = n_custom > 0
    chosen_col = "custom_cell_type_file" if chose_custom else base_col
    chosen_tool = "custom" if chose_custom else source_tool

    # Confidence from the source score when the source produced one.
    confidence_col = None
    score_base = base_col[:-len("_pruned")] if base_col.endswith("_pruned") else base_col
    candidate = _score_column_for(score_base)
    if candidate in adata.obs.columns:
        confidence_col = candidate
        adata.obs["cell_type_confidence"] = adata.obs[candidate]

    # Per-tool main label columns used to compute agreement against the canonical.
    tool_labels = {}
    for tool, prefix in TOOL_PREFIX.items():
        cols = merge_report.get(tool, {}).get("columns", []) if tool in merge_report else []
        labels = _label_columns(cols, prefix)
        ref_col = labels[0] if labels else None
        if ref_col and ref_col in adata.obs.columns:
            tool_labels[tool] = ref_col

    canonical = adata.obs["cell_type"]
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

    mean_agreement = (
        agreement[canonical.notna()].mean() if annotated else 0.0
    )

    info = {
        "cell_type_column": base_col,
        "cell_type_source_tool": source_tool,
        "cell_type_used": chosen_col,
        "cell_type_used_tool": chosen_tool,
        "cell_type_confidence_column": confidence_col,
        "custom_cell_type_file": os.path.basename(custom_file) if custom_file else None,
        "n_custom_assigned": n_custom,
        "tools_compared": {t: c for t, c in tool_labels.items()},
        "min_cell_type_annotated_pct": min_pct,
        "percent_annotated": round(final_pct * 100, 2) if total else 0.0,
        "percent_agreement": round(
            float(mean_agreement if mean_agreement == mean_agreement else 0.0), 4
        ),
        "n_cells_no_cell_type": int(total - annotated),
    }
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

    cell_type_info = add_canonical_cell_type(
        adata,
        merge_report,
        args.cell_type_column,
        args.min_cell_type_annotated_pct,
        args.custom_cell_type_file,
    )

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
