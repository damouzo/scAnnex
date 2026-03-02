# Contrasts CSV Schema for Differential Gene Expression

## Overview

The `contrasts.csv` file defines pairwise comparisons for differential expression analysis in scAnnex. This follows nf-core best practices for reproducible analysis specification.

## Schema Definition

### Required Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `contrast_id` | string | Unique identifier for the contrast | `TreatedVsControl` |
| `variable` | string | Column name in `adata.obs` to compare | `condition` |
| `group1` | string | Test group (numerator in fold-change) | `treated` |
| `group2` | string | Reference group (denominator in fold-change) | `control` |

### Optional Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `filter_column` | string | Column to filter cells before comparison | `leiden_res_0.5` |
| `filter_value` | string | Value to filter by (subset of cells) | `3` |

## Example Contrasts

### Example 1: Simple Condition Comparison

```csv
contrast_id,variable,group1,group2,filter_column,filter_value
TreatedVsControl,condition,treated,control,,
```

**Interpretation:** Compare all cells from `treated` vs `control` conditions.

---

### Example 2: Per-Cluster Comparisons

```csv
contrast_id,variable,group1,group2,filter_column,filter_value
TreatedVsControl_Cluster0,condition,treated,control,leiden_res_0.5,0
TreatedVsControl_Cluster1,condition,treated,control,leiden_res_0.5,1
TreatedVsControl_Cluster2,condition,treated,control,leiden_res_0.5,2
```

**Interpretation:** Compare `treated` vs `control` within each cluster separately (pseudo-bulk approach).

---

### Example 3: Cell Type Comparisons

```csv
contrast_id,variable,group1,group2,filter_column,filter_value
Bcells_vs_Tcells,celltypist_majority_voting,B cells,T cells,,
CD8_vs_CD4,celltypist_majority_voting,CD8 T cells,CD4 T cells,,
```

**Interpretation:** Compare annotated cell types directly.

---

### Example 4: Complex Multi-Contrast

```csv
contrast_id,variable,group1,group2,filter_column,filter_value
TreatedVsControl,condition,treated,control,,
HighVsLow,condition,high_dose,low_dose,,
TreatedVsControl_Bcells,condition,treated,control,celltypist_majority_voting,B cells
TreatedVsControl_Tcells,condition,treated,control,celltypist_majority_voting,T cells
Cluster3_vs_Cluster1,leiden_res_0.5,3,1,,
```

**Interpretation:** Multiple comparisons in a single analysis run.

---

## Validation Rules

The pipeline validates the contrasts file with these checks:

1. **Required columns:** `contrast_id`, `variable`, `group1`, `group2` must be present
2. **Unique IDs:** All `contrast_id` values must be unique
3. **Column existence:** All `variable` columns must exist in `adata.obs`
4. **Group existence:** `group1` and `group2` values must exist in the specified `variable` column
5. **Filter validity:** If `filter_column` is specified, `filter_value` must also be provided
6. **Filter column existence:** All `filter_column` values must exist in `adata.obs`
7. **Filter value existence:** All `filter_value` values must exist in the specified `filter_column`

## Output Files

For each contrast, the pipeline generates:

```
results/differential_expression/
├── {contrast_id}_results.csv          # Full DE results table
├── {contrast_id}_top100.csv           # Top 100 DE genes
├── {contrast_id}_volcano.pdf          # Volcano plot
└── {contrast_id}_stats.txt            # Summary statistics
```

## Integration with Nextflow Parameters

If no `contrasts.csv` is provided, the pipeline auto-generates a simple contrast using:

```bash
nextflow run main.nf \
  --run_dge true \
  --dge_groupby condition \
  --dge_reference control
```

This creates: `{group}_vs_{reference}` for all unique values in the `dge_groupby` column.

## Best Practices

1. **Naming:** Use descriptive `contrast_id` names (e.g., `TreatedVsControl` not `c1`)
2. **Replicates:** Ensure sufficient biological replicates per group (≥3 recommended)
3. **Cell counts:** Each group should have ≥30 cells (ideally >100)
4. **Per-cluster:** Use per-cluster comparisons for cell-type-specific effects
5. **Multiple testing:** Adjust p-values when running many contrasts

## See Also

- `data_demo/MultiSample/contrasts.csv` - Example contrasts file
- `bin/validate_contrasts.py` - Validation script
- `bin/differential_expression.py` - DGE implementation
