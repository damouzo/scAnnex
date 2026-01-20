# Conda Environment Strategy for scAnnex

## Single Unified Environment

We use **ONE** conda environment file (`env/scanpy.yml`) for the entire pipeline to ensure:
- ✅ Consistency across all modules
- ✅ Reproducibility on any machine
- ✅ No version conflicts between processes
- ✅ Simplified maintenance

---

## Environment: `scannex`

**File:** `env/scanpy.yml`

**Python Version:** 3.10.x (stable, widely supported)

**Core Packages:**
- Scanpy 1.9.8 (latest stable)
- AnnData 0.10.5 (compatible with Scanpy 1.9.8)
- CellTypist 1.6.2 (cell type annotation)
- Scrublet 0.2.3 (doublet detection)
- HarmonyPy 0.0.9 (batch correction)

---

## What Was Removed and Why

### ❌ scanpy-scripts (1.9.301)
**Problem:** Forces downgrade to Scanpy 1.9.3
**Solution:** Our custom `bin/` scripts provide all CLI functionality

### ❌ loompy (3.0.7)
**Problem:** Not available in current conda channels
**Solution:** Pipeline uses H5AD (native AnnData format)

### ❌ R packages (rpy2, r-seurat, r-seuratdisk)
**Problem:** Complex dependencies, conda solver conflicts
**Solution:** 
- For RDS conversion: Use separate Docker container
- See `containers/Dockerfile.seurat` for R environment
- Most users work with H5AD/MTX, not RDS

### ❌ plotly (5.18.0)
**Problem:** Large dependency tree, slows environment creation
**Solution:** Only needed in dashboard (separate container)

### ❌ jupyter, ipykernel
**Problem:** Not required for automated pipeline
**Solution:** Install separately if needed: `conda install jupyter -n scannex`

---

## Version Pinning Strategy

| Package | Pin Level | Example | Reason |
|---------|-----------|---------|--------|
| Critical | Major.Minor | `scanpy=1.9.8` | Exact version for reproducibility |
| Core | Major.Minor | `numpy=1.26.4` | API stability within minor versions |
| Tools | Major | `matplotlib-base=3.8` | Allow patch updates |
| Utils | Major | `tqdm=4.66` | Non-critical, flexible |

---

## Creating the Environment

### First Time Setup
```bash
cd /home/damo/scAnnex
mamba env create -f env/scanpy.yml
conda activate scannex
```

### Updating Existing Environment
```bash
mamba env update -n scannex -f env/scanpy.yml --prune
```

### Complete Reset
```bash
mamba env remove -n scannex
mamba env create -f env/scanpy.yml
```

---

## Nextflow Integration

All modules reference the same environment file:

```groovy
process EXAMPLE_MODULE {
    conda "${projectDir}/env/scanpy.yml"
    container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"
    
    // ... rest of process definition
}
```

**Why both conda and container?**
- `conda`: For local development and systems without containers
- `container`: For HPC/Cloud with Singularity/Docker

Nextflow automatically chooses based on active profile:
- `-profile conda` → uses conda environment
- `-profile docker` → uses container
- `-profile singularity` → uses container

---

## Handling RDS Input Files

If you need to process Seurat RDS files:

### Option 1: Separate R Environment (Recommended for local)
```bash
# Create R-specific environment
mamba create -n scannex-r \
    r-base=4.3 \
    r-seurat=5.0 \
    bioconductor-singlecellexperiment \
    rpy2=3.5

# Activate for RDS conversion only
conda activate scannex-r
Rscript bin/convert_rds_to_h5ad.R input.rds output.h5ad

# Switch back to main environment
conda activate scannex
```

### Option 2: Docker Container (Recommended for production)
```bash
# Use pre-built container with R + Python
docker run -v $(pwd):/data scannex/seurat:latest \
    Rscript /opt/scannex/bin/convert_rds_to_h5ad.R \
    /data/input.rds /data/output.h5ad
```

---

## Troubleshooting

### Conda Solver Hangs
```bash
# Use mamba (faster solver)
mamba env create -f env/scanpy.yml

# If still hanging, check channels
conda config --show channels
# Should be: conda-forge, bioconda, defaults
```

### Package Not Found
```bash
# Update conda
conda update -n base conda mamba

# Clear cache
conda clean --all

# Retry
mamba env create -f env/scanpy.yml
```

### Version Conflicts
```bash
# Check which package is causing conflict
mamba create -n test-env --dry-run -f env/scanpy.yml

# Relax version constraint if needed
# e.g., change `scipy=1.12.0` to `scipy=1.12`
```

---

## Testing the Environment

After creation, verify all imports work:

```bash
conda activate scannex

python << 'EOF'
import sys
import scanpy as sc
import anndata as ad
import scrublet
import harmonypy
import celltypist

print("✓ Environment Test Passed")
print(f"  Python: {sys.version.split()[0]}")
print(f"  Scanpy: {sc.__version__}")
print(f"  AnnData: {ad.__version__}")
print(f"  Scrublet: {scrublet.__version__}")
print(f"  CellTypist: {celltypist.__version__}")
EOF
```

Expected output:
```
✓ Environment Test Passed
  Python: 3.10.19
  Scanpy: 1.9.8
  AnnData: 0.10.5.post1
  Scrublet: 0.2.3
  CellTypist: 1.6.2
```

---

## Environment Size

**Disk Space Required:** ~2.5 GB

**Installation Time:**
- With mamba: 5-10 minutes
- With conda: 15-30 minutes

**Tip:** Use mamba for faster environment creation:
```bash
conda install -n base mamba
mamba env create -f env/scanpy.yml
```

---

## CI/CD Integration

For automated testing:

```yaml
# .github/workflows/test.yml
- name: Setup Conda
  uses: conda-incubator/setup-miniconda@v2
  with:
    environment-file: env/scanpy.yml
    activate-environment: scannex
    auto-activate-base: false

- name: Test Pipeline
  run: |
    conda activate scannex
    nextflow run main.nf -profile test,conda
```

---

## Summary

✅ **Single environment file** eliminates conflicts
✅ **No R dependencies** speeds up environment creation
✅ **Version pinned** for reproducibility
✅ **Tested on Python 3.10** for stability
✅ **Works with Nextflow** automatic environment creation
✅ **Documented alternatives** for edge cases (RDS, Jupyter)

**Last Updated:** 2026-01-20
**Tested On:** Ubuntu 18.04 (WSL2), CentOS 7, macOS 14
**Nextflow Versions:** 23.04+ (DSL2)
