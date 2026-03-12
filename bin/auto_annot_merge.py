#!/usr/bin/env python3

import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
    ad.settings.allow_write_nullable_strings = True


def parse_args():
    parser = argparse.ArgumentParser(description="Merge automatic annotation CSV outputs into one H5AD")
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
    return parser.parse_args()


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

    adata.uns["auto_annotation_summary"] = {
        "status": {k: bool(v.get("success", False)) for k, v in status.items()},
        "status_json": json.dumps(status, ensure_ascii=True),
        "columns_added": columns_added,
        "merge_report": merge_report,
    }
    adata.write_h5ad(args.output)

    summary = {
        "success": True,
        "status": status,
        "columns_added": columns_added,
        "merge_report": merge_report,
    }
    with open(args.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)


if __name__ == "__main__":
    main()
