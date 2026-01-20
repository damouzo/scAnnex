# Container Strategy for scAnnex SLC Production

## Overview

scAnnex SLC uses a **multi-level container strategy** to balance development speed, reproducibility, and production readiness.

---

## Current Implementation (Development Phase)

### 1. Conda Environments (Development & Testing)

**Purpose:** Fast iteration during development and testing on local machines/WSL2

**Location:** `env/scanpy.yml` and `env/scanpy-minimal.yml`

**Usage:**
```bash
# Full environment (includes R/Seurat for RDS conversion)
mamba env create -f env/scanpy.yml
conda activate scannex-scanpy

# Minimal environment (faster install, Python-only)
mamba env create -f env/scanpy-minimal.yml
conda activate scannex-minimal
```

**When to use:**
- Local development and debugging
- WSL2/laptop environments
- Quick parameter testing
- Interactive analysis with Jupyter

**Limitations:**
- Less reproducible across systems
- Conda solver can change package versions
- Larger disk footprint

---

### 2. Biocontainers (Current Default)

**Purpose:** Validated, community-maintained containers with specific tool versions

**Current Container:**
```groovy
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"
```

**Why Scanpy 1.7.2?**
- Stable, well-tested version
- Available in Biocontainers registry
- Compatible with most datasets

**Note:** We specify Scanpy 1.9.8 in conda environments for newer features, but fall back to 1.7.2 in containers for stability.

**When to use:**
- HPC environments with Singularity/Apptainer
- Production runs requiring reproducibility
- Multi-user shared environments

---

## Production Strategy (Deployment Phase)

### 3. Custom Docker Images

**Purpose:** Full control over dependencies, optimized for scAnnex SLC

**Dockerfile:** `containers/Dockerfile.scanpy`

**Build command:**
```bash
docker build -t scannex/scanpy:1.9.8 -f containers/Dockerfile.scanpy .
```

**Push to registry:**
```bash
# Docker Hub
docker tag scannex/scanpy:1.9.8 damouzo/scannex:1.9.8
docker push damouzo/scannex:1.9.8

# GitHub Container Registry (recommended for private repos)
docker tag scannex/scanpy:1.9.8 ghcr.io/damouzo/scannex:1.9.8
docker push ghcr.io/damouzo/scannex:1.9.8
```

**Update modules to use custom container:**
```groovy
// In modules/local/unify_input.nf
container "ghcr.io/damouzo/scannex:1.9.8"
```

---

### 4. Singularity/Apptainer Images

**Purpose:** HPC-compatible containerization with better security model

**Convert Docker to Singularity:**
```bash
# Option 1: Build from Docker Hub
singularity build scannex_1.9.8.sif docker://damouzo/scannex:1.9.8

# Option 2: Build from local Docker daemon
singularity build scannex_1.9.8.sif docker-daemon://scannex/scanpy:1.9.8
```

**Use in Nextflow config:**
```groovy
profiles {
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        singularity.cacheDir = "${projectDir}/singularity-cache"
        
        process {
            container = "${projectDir}/containers/scannex_1.9.8.sif"
        }
    }
}
```

---

## Version Pinning Strategy

### Rationale

For SLC production, **all dependencies must be pinned** to ensure:
- Reproducibility across runs
- Consistent results for publications
- Traceable analysis provenance

### Implementation

| Component | Development | Production |
|-----------|-------------|------------|
| Python | `python=3.10` | `python=3.10.13` |
| Scanpy | `scanpy=1.9.8` | `scanpy==1.9.8` |
| AnnData | `anndata=0.10.5` | `anndata==0.10.5` |
| NumPy | `numpy=1.26.4` | `numpy==1.26.4` |

**Conda:** Use `=` (allows patch updates within minor version)
**Pip/Docker:** Use `==` (exact version, no updates)

---

## Recommended Workflow by Environment

### Local Development (WSL2/Laptop)
```bash
nextflow run main.nf -profile conda,laptop --input data.csv
```
- Uses: Conda environment (`env/scanpy.yml`)
- Fast: No container overhead
- Flexible: Easy to add packages for debugging

### HPC (SLURM/PBS)
```bash
nextflow run main.nf -profile singularity,cluster --input data.csv
```
- Uses: Singularity `.sif` images
- Reproducible: Identical software stack
- Compatible: Works without root access

### Cloud (AWS/GCP/Azure)
```bash
nextflow run main.nf -profile docker,cloud --input s3://bucket/data.csv
```
- Uses: Docker containers
- Scalable: Easy to spawn multiple instances
- Cost-effective: Pay per use

### CI/CD Testing
```bash
nextflow run main.nf -profile docker,test
```
- Uses: Docker containers
- Fast: Cached layers speed up builds
- Isolated: No system dependencies

---

## Migration Path to Production

### Phase 1: Development (Current)
- ✅ Conda environments for local testing
- ✅ Biocontainers as fallback
- ✅ Version-pinned environment files

### Phase 2: Pre-Production (Next 2-4 weeks)
- [ ] Build custom Docker images
- [ ] Push to container registry (GHCR or Docker Hub)
- [ ] Convert to Singularity for HPC testing
- [ ] Update all modules to use custom containers
- [ ] Test on HPC cluster

### Phase 3: Production (Release)
- [ ] Lock all container tags (no `:latest`)
- [ ] Create immutable release containers (`scannex:v0.1.0`)
- [ ] Store `.sif` files in shared HPC storage
- [ ] Document container provenance in methods
- [ ] Set up automated container rebuilds on dependency updates

---

## Container Registry Strategy

### Recommended: GitHub Container Registry (GHCR)

**Advantages:**
- Free for public repositories
- Native GitHub integration
- Good for academic/open-source projects
- Fine-grained access control

**Setup:**
```bash
# Login to GHCR
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin

# Tag and push
docker tag scannex/scanpy:1.9.8 ghcr.io/damouzo/scannex:1.9.8
docker push ghcr.io/damouzo/scannex:1.9.8
```

### Alternative: Docker Hub

**Advantages:**
- Most widely used
- Good documentation
- Easy discovery

**Disadvantages:**
- Rate limits on free tier
- Requires separate account

---

## Testing Strategy

### Container Validation Checklist

Before deploying a new container to production:

- [ ] All Python packages import successfully
- [ ] Scanpy version matches specification
- [ ] Scripts in `bin/` execute without errors
- [ ] Memory footprint acceptable (<2GB base)
- [ ] Security scan passes (Trivy/Snyk)
- [ ] Documented in CHANGELOG.md

### Automated Testing

```bash
# Test Docker container
docker run --rm scannex/scanpy:1.9.8 python -c "import scanpy; print(scanpy.__version__)"

# Test Singularity container
singularity exec scannex_1.9.8.sif python -c "import scanpy; print(scanpy.__version__)"
```

---

## Troubleshooting

### Conda environment creation fails

**Problem:** Conda solver hangs or fails with conflicts

**Solution:**
```bash
# Use Mamba (faster solver)
mamba env create -f env/scanpy.yml

# Or use minimal environment
mamba env create -f env/scanpy-minimal.yml
```

### Container pull fails on HPC

**Problem:** No internet access on compute nodes

**Solution:**
```bash
# Pre-pull on login node
singularity pull docker://damouzo/scannex:1.9.8

# Or use local .sif file
singularity.cacheDir = "/shared/containers"
```

### Version mismatch between Conda and containers

**Problem:** Results differ between local (Conda) and HPC (Singularity)

**Solution:**
- Always use containers for production runs
- Use Conda only for development/debugging
- Document which environment was used in analysis logs

---

## Summary

| Aspect | Development | Production |
|--------|-------------|------------|
| **Container Type** | Conda | Docker/Singularity |
| **Version Pinning** | Minor versions | Exact versions |
| **Registry** | conda-forge/bioconda | GHCR/Docker Hub |
| **Update Frequency** | Weekly | Only for critical fixes |
| **Testing Level** | Manual | Automated CI/CD |

**Current Status:** Phase 1 (Development)
**Next Milestone:** Build and test custom Docker images

---

**Last Updated:** 2026-01-20
**Pipeline Version:** 0.1.0
