#!/usr/bin/env python3

import argparse
import json
import re
import warnings

import anndata as ad
import celltypist
import pandas as pd
import scanpy as sc
from celltypist import models

if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
    ad.settings.allow_write_nullable_strings = True

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(description="CellTypist annotation to metadata CSV")
    parser.add_argument("--input", required=True, help="Input H5AD")
    parser.add_argument("--output", required=True, help="Output CSV with annotation columns")
    parser.add_argument("--status", required=True, help="Status JSON output")
    parser.add_argument(
        "--models",
        default="Immune_All_Low.pkl",
        help="Comma-separated CellTypist models",
    )
    parser.add_argument("--majority-voting", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def slugify_model(model_name):
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", model_name.strip().lower())
    slug = re.sub(r"_+", "_", slug).strip("_")
    return slug


def load_data(h5ad_path):
    adata = sc.read_h5ad(h5ad_path)
    if "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"].copy()
    elif "normalized" in adata.layers:
        adata.X = adata.layers["normalized"].copy()
    return adata


def run_single_model(adata, model_name, majority_voting):
    model_obj = models.Model.load(model_name)
    predictions = celltypist.annotate(
        adata,
        model=model_obj,
        majority_voting=majority_voting,
    )
    if majority_voting and hasattr(predictions.predicted_labels, "majority_voting"):
        labels = predictions.predicted_labels.majority_voting
    else:
        labels = predictions.predicted_labels.predicted_labels
    scores = predictions.probability_matrix.max(axis=1)
    return labels, scores


def main():
    args = parse_args()
    models_list = [m.strip() for m in args.models.split(",") if m.strip()]
    adata = load_data(args.input)

    out_df = pd.DataFrame(index=adata.obs_names)
    out_df.index.name = "cell_id"
    status = {
        "tool": "celltypist",
        "success": True,
        "models": [],
        "errors": [],
    }

    for model_name in models_list:
        model_slug = slugify_model(model_name)
        label_col = f"auto_annot_celltypist_{model_slug}"
        score_col = f"{label_col}_score"
        try:
            labels, scores = run_single_model(adata, model_name, args.majority_voting)
            out_df[label_col] = labels.values
            out_df[score_col] = scores.values
            status["models"].append({"model": model_name, "success": True})
        except Exception as exc:
            status["models"].append({"model": model_name, "success": False})
            status["errors"].append(f"{model_name}: {exc}")
            if not args.continue_on_error:
                status["success"] = False
                break

    if len(out_df.columns) == 0:
        status["success"] = False

    out_df.reset_index().to_csv(args.output, index=False)
    with open(args.status, "w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2)

    if not status["success"] and not args.continue_on_error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
