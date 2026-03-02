# Dashboard Environment Update Instructions

## Issue

The dashboard environment is missing two required R packages for the DGE visualization:
- `r-ggrepel` (for non-overlapping gene labels on volcano plots)
- `r-dplyr` (for data manipulation)

## Solution

You have two options to update the environment:

### Option 1: Update Existing Environment (Faster)

Run this command to install the missing packages:

```bash
conda install -n scannex-dashboard -c conda-forge r-ggrepel r-dplyr -y
```

This takes ~2-3 minutes.

### Option 2: Recreate Environment from Updated YAML (Clean)

If you prefer a clean installation:

```bash
# Remove old environment
conda env remove -n scannex-dashboard

# Create new environment from updated YAML
cd /data/BCI-KRP/projects/scAnnex/dashboard
conda env create -f environment_dashboard.yml -n scannex-dashboard
```

This takes ~5-10 minutes.

## After Updating

Launch the dashboard normally:

```bash
cd /data/BCI-KRP/projects/scAnnex/dashboard
bash launch_dashboard_hpc.sh /data/BCI-KRP/projects/scAnnex/results_dge_test
```

The DGE tab should now work correctly.

## Changes Made

The following files were updated to add DGE visualization:

1. **dashboard/environment_dashboard.yml** - Added r-ggrepel and r-dplyr
2. **dashboard/ui.R** - Added DGE tab with volcano plot and controls
3. **dashboard/server.R** - Added DGE reactive logic and plotting
4. **dashboard/global.R** - Added ggrepel and dplyr library imports
5. **dashboard/test_dge_load.R** - Created test script for validation

## Testing the DGE Tab

Once the dashboard loads:

1. Navigate to "Differential Expression" tab
2. Enter DGE results directory: `results_dge_test/dge/dge_results`
3. Click "Load DGE Results"
4. Select contrast from dropdown (e.g., `treated_vs_control`)
5. Adjust filtering thresholds:
   - P-value threshold: 0.05 (default)
   - Log2 FC threshold: 0.25 (default)
6. View interactive volcano plot
7. Check significant genes table below
8. Test download buttons for plots and CSVs

## Expected Results

With the demo data (treated vs control), you should see:
- ~43 significant genes (p < 0.05, |log2FC| >= 0.25)
- Volcano plot with upregulated (red) and downregulated (blue) genes
- Top genes labeled on the plot
- Downloadable high-resolution plots and CSV files
