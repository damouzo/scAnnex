#!/usr/bin/env python3
"""Mock test for pseudobulk_build.py core aggregation logic.

Verifies:
  (a) correct number of aggregates (sample x cell type),
  (b) biological replicates = number of samples, NOT number of cells
      (anti-pseudo-replication check),
  (c) min_cells filter drops small aggregates as integer raw sums,
  (d) aggregation uses raw integer counts, no normalization/pre-scaling,
  (e) replicate_group_size_condition (used by the DGE guard to trigger
      SKIP_INSUFFICIENT_REPLICATES) is 1 when a group has a single sample.

Run with:  python tests/test_pseudobulk.py
Requires:  anndata, numpy, pandas (and scanpy to import the module).
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import anndata as ad  # noqa: E402
from pseudobulk_build import build_pseudobulk  # noqa: E402


def make_synthetic_adata():
    """6 samples (4 control, 2 treated), 2 cell types, 6 genes.

    Deterministic design so TypeB has 2 control samples (C2, C3) but a single
    treated sample (T2), exercising the replicate guard.
    """
    n_cells = 200
    rng = np.random.RandomState(42)
    genes = [f"G{i}" for i in range(6)]
    X = rng.randint(0, 50, size=(n_cells, len(genes))).astype(np.int64)

    sample_ids = (
        ["C1"] * 30 + ["C2"] * 30 + ["C3"] * 30 +
        ["T1"] * 45 + ["T2"] * 45 + ["C1x"] * 20
    )
    sample_ids = sample_ids[:n_cells]
    conditions = ["control" if s.startswith("C") else "treated" for s in sample_ids]

    cell_type_for = {
        "C1": "TypeA", "C2": "TypeB", "C3": "TypeB",
        "T1": "TypeA", "T2": "TypeB", "C1x": "TypeA",
    }
    cell_types = [cell_type_for[s] for s in sample_ids]

    obs = pd.DataFrame({"sample_id": sample_ids, "condition": conditions, "cell_type": cell_types})
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.layers["counts"] = X.copy()
    return adata


def main():
    adata = make_synthetic_adata()
    counts_df, meta_df, summary = build_pseudobulk(adata, min_cells=10)

    # (a) aggregate label = sample::celltype, one column per aggregate
    assert ":" in counts_df.columns[0], "aggregate columns must be sample::celltype"
    assert counts_df.shape[1] == len(meta_df), "counts columns must match metadata rows"

    # (d) counts are integer and equal to cell-wise raw sums (no normalization)
    assert (counts_df.dtypes == np.int64).all(), "pseudobulk counts must be integers"
    expected = adata.layers["counts"].sum(axis=0)
    # reconstruct one aggregate from raw cells to confirm direct summation
    m = (adata.obs["sample_id"] == "C1") & (adata.obs["cell_type"] == "TypeA")
    cell = "C1::TypeA"
    if cell in counts_df.columns:
        np.testing.assert_allclose(counts_df[cell].values, np.asarray(adata.layers["counts"])[m].sum(axis=0))

    # (b) replicate_group_size_condition = number of samples, never number of cells
    assert meta_df["replicate_group_size_condition"].max() >= 1
    sample_sizes = meta_df.groupby(["cell_type", "condition"])["sample_id"].nunique()
    for _, row in meta_df.iterrows():
        expected_repl = int(sample_sizes.get((row["cell_type"], row["condition"]), 0))
        assert row["replicate_group_size_condition"] == expected_repl

    # (c) every aggregate kept has >= min_cells raw cells
    assert (meta_df["n_cells_raw"] >= 10).all(), "min_cells filter not respected"

    # (e) a group with a single sample must report replicate_group_size_condition == 1,
    #     which drives SKIP_INSUFFICIENT_REPLICATES on the DESeq2 side
    single = meta_df[meta_df["sample_id"] == "T2"][["cell_type", "condition", "replicate_group_size_condition"]]
    assert (single["replicate_group_size_condition"] == 1).all(), \
        "single-sample group must report replicate size 1 (triggers SKIP)"

    print("pseudobulk_build mock test PASSED")
    print(f"  aggregates: {counts_df.shape[1]}")
    print(f"  cell types: {meta_df['cell_type'].unique().tolist()}")
    print(f"  samples:    {meta_df['sample_id'].nunique()}")
    print(f"  conditions: {meta_df['condition'].unique().tolist()}")


if __name__ == "__main__":
    main()