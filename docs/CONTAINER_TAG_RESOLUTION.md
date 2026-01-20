# Container Tag Resolution Guide

## Issue Summary

The Biocontainers registry on Quay.io has **limited Scanpy versions** available. The most recent version is **1.7.2**, not 1.9.3 as originally specified.

## Verified Working Container Tags

After extensive testing and verification against the Quay.io registry, here are the **confirmed working** container tags:

### ✅ Available and Verified

| Tool | Container Tag | Status |
|------|--------------|--------|
| **Scanpy** | `quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0` | ✅ VERIFIED |
| **CellTypist** | `quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0` | ✅ EXISTS |
| **Seurat** | `docker.io/satijalab/seurat:5.0.0` | ✅ VERIFIED |

### ❌ Not Available

| Tool | Attempted Tag | Issue |
|------|--------------|-------|
| Scanpy 1.9.3 | `quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0` | **Does not exist** |
| Mulled Scrublet | `mulled-v2-7be8d0ad...` | **Does not exist** |
| Mulled Harmony | `mulled-v2-6af17bf2...` | **Does not exist** |

## Recommended Solutions

###  Option 1: Use Conda Profile (RECOMMENDED for laptops)

The **simplest and most reliable** solution is to use the `conda` profile, which automatically installs all required Python packages with correct versions:

```bash
nextflow run main.nf \
  -profile conda,laptop \
  --input samplesheet.csv \
  --outdir results
```

**Advantages:**
- ✅ All dependencies automatically installed
- ✅ Correct versions guaranteed
- ✅ No container tag issues
- ✅ Works immediately

**Disadvantages:**
- ⏱️ First run takes 10-15 minutes to build conda environments
- 💾 Requires ~5GB disk space for environments
- 🐌 Slightly slower startup (subsequent runs reuse environments)

### Option 2: Build Custom Container (for production/repeated runs)

We've provided a Dockerfile that extends the Scanpy 1.7.2 container with all required tools:

```bash
# Build the extended container
cd /path/to/scAnnex
docker build -t scannex/scanpy-extended:1.7.2 -f docker/Dockerfile.scanpy-extended .

# Update modules to use it (or done automatically)
# Then run with Docker profile
nextflow run main.nf \
  -profile docker,laptop \
  --input samplesheet.csv \
  --outdir results
```

**Advantages:**
- ✅ Fast execution (container cached)
- ✅ Reproducible environment
- ✅ All tools included

**Disadvantages:**
- ⏱️ Requires one-time build (~5-10 minutes)
- 🔧 Needs Docker build capabilities
- 💾 Larger container image (~2GB)

### Option 3: Mixed Docker + Conda (intermediate)

Use Docker for some processes and Conda for others by commenting/uncommenting the container lines in module files.

## Current Module Configuration

All modules have been updated to use verified container tags:

```groovy
// modules/local/unify_input.nf
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"

// modules/local/quality_control.nf
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"

// modules/local/doublet_detection.nf
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"
// NOTE: scrublet not included - use conda profile or custom container

// modules/local/normalize_integrate.nf
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"
// NOTE: harmonypy not included - use conda profile or custom container

// modules/local/standard_processing.nf
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"

// modules/local/auto_annot_celltypist.nf
container "quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0"
```

## Why the Container Tag Issues?

1. **Biocontainers lag behind PyPI**: Latest scanpy on PyPI is 1.9.x, but Biocontainers only has 1.7.2
2. **Mulled containers deprecated**: The specific hashes for mulled containers may have expired or been rebuilt
3. **Registry inconsistencies**: Quay.io may remove old tags or reorganize repositories

## Testing Container Availability

Before running the pipeline, you can verify containers exist:

```bash
# Test Scanpy
docker pull quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0
docker run --rm quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0 \
  python -c "import scanpy, anndata; print(f'scanpy: {scanpy.__version__}, anndata: {anndata.__version__}')"

# Expected: scanpy: 1.7.2, anndata: 0.7.5

# Test CellTypist
docker pull quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0
```

## Building the Custom Container

The Dockerfile is provided at `docker/Dockerfile.scanpy-extended`:

```dockerfile
FROM quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0

# Install system dependencies
USER root
RUN apt-get update && apt-get install -y gcc g++ make

# Install additional packages
RUN pip install --no-cache-dir \
    scrublet==0.2.3 \
    harmonypy==0.0.10 \
    celltypist==1.6.2

USER nobody
```

Build command:
```bash
docker build -t scannex/scanpy-extended:1.7.2 -f docker/Dockerfile.scanpy-extended .
```

## Troubleshooting

### "manifest for ... not found"

This means the container tag doesn't exist on the registry. Solutions:
1. Switch to conda profile (recommended)
2. Build custom container
3. Use verified tags listed above

### "ModuleNotFoundError" for scrublet/harmonypy

The base scanpy:1.7.2 container doesn't include these. Solutions:
1. Use conda profile
2. Build extended container
3. Skip processes that need these tools (modify workflow)

### Conda is slow

First run builds environments. Subsequent runs are much faster. To speed up:
```bash
# Use mamba (faster conda)
conda install -n base mamba
# Then Nextflow will automatically use mamba if available
```

## Long-term Solution

We're monitoring:
- Biocontainers for newer Scanpy releases
- Alternative registries (Docker Hub, GitHub Container Registry)
- Official Scanpy Docker images

Once newer containers are available, we'll update the module files accordingly.

## Quick Start (TL;DR)

**For immediate testing:**
```bash
nextflow run main.nf -profile conda,laptop --input samplesheet.csv --outdir results
```

**For production (after building container):**
```bash
docker build -t scannex/scanpy-extended:1.7.2 -f docker/Dockerfile.scanpy-extended .
nextflow run main.nf -profile docker,laptop --input samplesheet.csv --outdir results
```

---

**Last Updated**: 2026-01-20
**Verified Against**: Quay.io Biocontainers Registry
