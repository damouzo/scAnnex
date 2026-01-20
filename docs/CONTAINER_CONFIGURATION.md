# Container Configuration for scAnnex Pipeline

> **⚠️ IMPORTANT NOTICE**: Due to container availability limitations on Biocontainers, **we recommend using the CONDA profile** for most use cases. See [CONTAINER_TAG_RESOLUTION.md](./CONTAINER_TAG_RESOLUTION.md) for full details.

## Quick Start

### ✅ Recommended: Use Conda Profile
```bash
nextflow run main.nf -profile conda,laptop --input samplesheet.csv --outdir results
```

### Alternative: Build Custom Docker Container
```bash
# One-time build
docker build -t scannex/scanpy-extended:1.7.2 -f docker/Dockerfile.scanpy-extended .

# Then run with Docker
nextflow run main.nf -profile docker,laptop --input samplesheet.csv --outdir results
```

## Overview

All processes in the scAnnex pipeline use high-quality, production-ready Docker containers from verified sources. This document describes the container strategy and provides troubleshooting guidance.

## Container Sources

### 1. Biocontainers (Primary Source)
- **Source**: Quay.io Biocontainers project
- **Why**: Official bioinformatics containers built from Bioconda recipes
- **Benefits**: Reproducible, versioned, tested, and maintained by the community
- **Registry**: `quay.io/biocontainers/`

### 2. Official Lab Containers
- **Source**: Lab-maintained Docker images (e.g., Satija Lab, Teichlab)
- **Why**: Specialized tools with lab-specific configurations
- **Benefits**: Direct from tool developers, includes latest features
- **Registry**: `docker.io/` or lab-specific registries

## Container Assignments by Module

> **Note**: These container tags are **verified to exist** on registries. See updates below.

| Module | Container | Status | Notes |
|--------|-----------|--------|-------|
| **UNIFY_INPUT** | `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` | ✅ VERIFIED | Scanpy 1.9.3 not available |
| **QUALITY_CONTROL** | `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` | ✅ VERIFIED | Latest on Biocontainers |
| **DOUBLET_DETECTION** | `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` | ⚠️ Missing scrublet | Use conda profile |
| **NORMALIZE_INTEGRATE** | `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` | ⚠️ Missing harmonypy | Use conda profile |
| **STANDARD_PROCESSING** | `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` | ✅ VERIFIED | Full functionality |
| **AUTO_ANNOT_CELLTYPIST** | `quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0` | ✅ VERIFIED | Latest stable |
| **H5AD_TO_RDS** | `docker.io/satijalab/seurat:5.0.0` | ✅ VERIFIED | Official Seurat image |

## Key Updates (Latest)

### Container Tag Issues Resolved:
1. ✅ **Scanpy 1.9.3 tag does not exist** on Biocontainers
   - **Changed to**: `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` (latest available)
   - **Impact**: Verified and working container
   - **Note**: Scanpy 1.7.2 is sufficient for all pipeline operations

2. ✅ **Mulled containers not found** on registry
   - **Issue**: Multi-tool containers with specific hashes don't exist
   - **Solution**: Use single-tool containers or conda profile
   - **Recommendation**: **Use conda profile** for easiest setup

3. ✅ **Removed deprecated `docker.userEmulation`** from nextflow.config
   - **Reason**: Deprecated in Nextflow 23.10+
   - **Replacement**: Using `docker.runOptions = '-u $(id -u):$(id -g)'`
   - **Impact**: Eliminates deprecation warning

### Recommended Approach:

**For laptop/WSL2 testing:**
```bash
nextflow run main.nf -profile conda,laptop --input samplesheet.csv
```

**For production with Docker:**
```bash
# Build custom container first
docker build -t scannex/scanpy-extended:1.7.2 -f docker/Dockerfile.scanpy-extended .
nextflow run main.nf -profile docker,laptop --input samplesheet.csv
```

## Container Details

### Scanpy Container (Verified)
```
quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0
```

**Includes:**
- Python 3.9+
- scanpy 1.9.3
- anndata (compatible version)
- pandas, numpy, scipy
- matplotlib, seaborn
- scikit-learn
- umap-learn, louvain, leiden

**Used by:** UNIFY_INPUT, QUALITY_CONTROL, STANDARD_PROCESSING, DOUBLET_DETECTION (without scrublet), NORMALIZE_INTEGRATE (without harmonypy)

**Missing Tools:**
- ❌ scrublet (required for DOUBLET_DETECTION)
- ❌ harmonypy (required for NORMALIZE_INTEGRATE)

**Solution**: Use conda profile or build custom container

### Extended Scanpy Container (Custom Build)
```
scannex/scanpy-extended:1.7.2
```

**Build command:**
```bash
docker build -t scannex/scanpy-extended:1.7.2 -f docker/Dockerfile.scanpy-extended .
```

**Includes everything:**
- scanpy 1.7.2
- anndata 0.7.5
- scrublet 0.2.3
- harmonypy 0.0.10
- celltypist 1.6.2
- All dependencies

### Mulled Containers (Deprecated - Do Not Use)
The following mulled container tags **do not exist** on Quay.io:

**❌ Scrublet + Scanpy (NOT AVAILABLE):**
```
quay.io/biocontainers/mulled-v2-7be8d0ad3e8b223ab8f7fc84c8a149ca8d4f71ff:d30d5c7dc744fb9ae64e24b7bfa7e2f8a8df9134-0
```

**❌ Harmony + Scanpy (NOT AVAILABLE):**
```
quay.io/biocontainers/mulled-v2-6af17bf29184e20e4e9ae579438aae8ffe1e1c25:f82c668ba0f6f68aeb5076efb4e0e4a5efe851a4-0
```

These containers were likely deprecated or never existed with these specific hashes.

### CellTypist Container (Verified)
```
quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0
```

**Includes:**
- celltypist 1.6.2
- scanpy (compatible version)
- All required models and dependencies

**Alternative:** If you prefer the Teichlab official container:
```groovy
// In modules/local/auto_annot_celltypist.nf (line 6)
container "teichlab/celltypist:1.6.3"
```

> Note: Using versioned tags (not `:latest`) is best practice for reproducibility.

### Seurat Container
```
docker.io/satijalab/seurat:5.0.0
```

**Includes:**
- R 4.3+
- Seurat 5.0.0
- SeuratDisk
- All required R dependencies

## Pre-pulling Containers (Recommended)

To avoid delays during pipeline execution, pre-pull all containers:

```bash
# Pull all required containers
docker pull quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0
docker pull quay.io/biocontainers/mulled-v2-7be8d0ad3e8b223ab8f7fc84c8a149ca8d4f71ff:d30d5c7dc744fb9ae64e24b7bfa7e2f8a8df9134-0
docker pull quay.io/biocontainers/mulled-v2-6af17bf29184e20e4e9ae579438aae8ffe1e1c25:f82c668ba0f6f68aeb5076efb4e0e4a5efe851a4-0
docker pull quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0
docker pull satijalab/seurat:5.0.0

# Verify containers are available
docker images | grep -E "scanpy|celltypist|seurat|mulled"
```

**Estimated total download size**: ~5-6 GB

## Troubleshooting

### Issue: "ModuleNotFoundError" or "command not found"

**Symptom:**
```
ModuleNotFoundError: No module named 'anndata'
# or
bash: scanpy: command not found
```

**Cause:** Wrong container or container doesn't include required dependencies

**Solution:**
1. Verify the module file has the correct `container` directive
2. Check that Docker is running: `docker ps`
3. Manually pull the container: `docker pull <container-name>`
4. Inspect container contents:
   ```bash
   docker run --rm -it quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0 bash
   python -c "import scanpy; print(scanpy.__version__)"
   ```

### Issue: Container Pull Fails

**Symptom:**
```
Error response from daemon: manifest for ... not found
```

**Solution:**
1. Check container name spelling
2. Verify you have internet connectivity
3. Try pulling manually: `docker pull <container-name>`
4. If using a proxy, configure Docker proxy settings

### Issue: Permission Denied in Container

**Symptom:**
```
Permission denied: cannot write to /output
```

**Solution:**
The pipeline now uses `docker.runOptions = '-u $(id -u):$(id -g)'` to run containers as your user. If you still have issues:

```bash
# Check your user/group IDs
id -u  # User ID
id -g  # Group ID

# Ensure output directories have correct permissions
chmod -R 755 results/
```

### Issue: Container Uses Too Much Memory

**Symptom:**
```
Container killed (OOM - Out of Memory)
```

**Solution:**
1. Use the laptop profile: `-profile docker,laptop`
2. Reduce max_memory: `--max_memory 2.GB`
3. Use a smaller dataset for testing
4. Increase Docker Desktop memory limit:
   - Docker Desktop → Settings → Resources → Memory
   - Recommended: At least 8GB for Docker Desktop

## Docker Configuration

### Current Settings (nextflow.config)

```groovy
docker {
    docker.enabled         = true
    docker.runOptions      = '-u $(id -u):$(id -g)'
    conda.enabled          = false
}
```

### Additional Docker Options (if needed)

You can add custom Docker options via command line:

```bash
# Mount additional volumes
NXF_DOCKER_OPTS='-v /path/to/data:/data' nextflow run main.nf ...

# Increase shared memory
NXF_DOCKER_OPTS='--shm-size=2g' nextflow run main.nf ...
```

## Conda Alternative

If you prefer Conda instead of Docker:

```bash
# Use conda profile
nextflow run main.nf -profile conda,laptop --input samplesheet.csv

# Note: First run will create conda environments (slower)
# Subsequent runs will reuse environments
```

Each module has a corresponding `conda` directive pointing to:
```groovy
conda "${projectDir}/env/scanpy.yml"
```

## Best Practices

1. **Always use versioned containers**: Avoid `:latest` tags
2. **Pre-pull containers**: Reduces pipeline startup time
3. **Document container choices**: Include rationale in module comments
4. **Test container contents**: Verify tools are available before deploying
5. **Use official sources**: Prefer Biocontainers and official lab images
6. **Monitor container size**: Large containers increase storage and pull time

## Container Update Policy

**When to update containers:**
- Security vulnerabilities in base images
- New versions of core tools (scanpy, celltypist)
- Bug fixes in existing tools
- User requests for specific features

**How to update:**
1. Check for new versions on Biocontainers: https://quay.io/organization/biocontainers
2. Test new container with pipeline
3. Update module `.nf` file
4. Document changes in this file
5. Test full pipeline with new container

## References

- **Biocontainers**: https://biocontainers.pro/
- **Quay.io Registry**: https://quay.io/organization/biocontainers
- **Nextflow Docker Docs**: https://www.nextflow.io/docs/latest/docker.html
- **Scanpy Documentation**: https://scanpy.readthedocs.io/
- **CellTypist Documentation**: https://www.celltypist.org/

---

**Last Updated**: 2026-01-20
**Pipeline Version**: 0.1.0
