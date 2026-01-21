# Using scAnnex with Singularity/Apptainer

**Recommended for:** HPC clusters, production environments, reproducible science

Singularity/Apptainer is the **recommended container runtime** for scAnnex because:
- ✅ Designed for HPC and shared computing environments
- ✅ No root/sudo required after installation
- ✅ Better security model than Docker
- ✅ Automatic image caching (saves bandwidth)
- ✅ Native MPI support for parallel computing
- ✅ Standard in most research clusters worldwide

---

## Quick Start

```bash
# Run scAnnex with Singularity (recommended)
nextflow run main.nf \
  --input data.h5ad \
  --outdir results \
  -profile singularity

# Or with Apptainer (newer name)
nextflow run main.nf \
  --input data.h5ad \
  --outdir results \
  -profile apptainer
```

**First run will be slow** (downloading/converting containers ~5-10 min)  
**Subsequent runs are fast** (containers cached)

---

## Installation

### Option 1: Pre-installed (HPC Clusters)

Most research clusters have Singularity/Apptainer pre-installed:

```bash
# Check if available
module avail singularity
module avail apptainer

# Load module
module load singularity/3.8.0
# or
module load apptainer/1.2.0

# Verify
singularity --version
```

### Option 2: Install on Ubuntu/Debian (WSL2)

```bash
# Install dependencies
sudo apt update
sudo apt install -y \
  build-essential \
  libseccomp-dev \
  pkg-config \
  squashfs-tools \
  cryptsetup \
  wget

# Download and install Apptainer (recommended - newer)
APPTAINER_VERSION=1.2.5
wget https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/apptainer_${APPTAINER_VERSION}_amd64.deb
sudo dpkg -i apptainer_${APPTAINER_VERSION}_amd64.deb

# Or install SingularityCE (classic)
SINGULARITY_VERSION=3.11.4
wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-ce_${SINGULARITY_VERSION}-focal_amd64.deb
sudo dpkg -i singularity-ce_${SINGULARITY_VERSION}-focal_amd64.deb

# Verify installation
singularity --version
# or
apptainer --version
```

### Option 3: Install from Source (Advanced)

See official documentation:
- Apptainer: https://apptainer.org/docs/admin/main/installation.html
- Singularity: https://sylabs.io/guides/latest/admin-guide/

---

## How It Works

### Automatic Container Management

scAnnex uses Nextflow's built-in Singularity support:

1. **First run:** Nextflow automatically:
   - Downloads Docker images from registries
   - Converts to Singularity `.sif` format
   - Caches in `~/.singularity/cache/` or `~/.apptainer/cache/`

2. **Subsequent runs:**
   - Uses cached `.sif` images
   - No re-download or conversion needed
   - Fast execution

### Container Sources

All scAnnex containers come from public registries:

```groovy
// In modules/local/*.nf
container "quay.io/biocontainers/scanpy:1.10.0--pyhdfd78af_0"
```

Nextflow automatically converts to:
```
~/.singularity/cache/quay.io-biocontainers-scanpy-1.10.0--pyhdfd78af_0.img
```

---

## Cache Management

### Default Cache Locations

```bash
# Singularity
~/.singularity/cache/

# Apptainer  
~/.apptainer/cache/
```

### Custom Cache Directory

```bash
# Set custom location (useful for shared filesystems)
export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache

# Or in nextflow.config
singularity.cacheDir = '/scratch/shared/singularity_cache'
```

### Check Cache Size

```bash
du -sh ~/.singularity/cache/
# Expected: 2-5 GB for all scAnnex containers
```

### Clean Cache

```bash
# Remove all cached images
rm -rf ~/.singularity/cache/*

# Next pipeline run will re-download containers
```

---

## Troubleshooting

### Issue: "Failed to pull image"

**Cause:** Network issues or registry unavailable

**Solution:**
```bash
# Pre-pull containers manually
singularity pull docker://quay.io/biocontainers/scanpy:1.10.0--pyhdfd78af_0

# Or use a different cache directory with more space
export SINGULARITY_CACHEDIR=/tmp/singularity_cache
```

### Issue: "No space left on device"

**Cause:** Cache directory full (often /tmp)

**Solution:**
```bash
# Check available space
df -h ~/.singularity/cache

# Move cache to larger filesystem
export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR
```

### Issue: "Operation not permitted"

**Cause:** Running on NFS-mounted home directory

**Solution:**
```bash
# Use local scratch space
export SINGULARITY_TMPDIR=/tmp/$USER/singularity_tmp
mkdir -p $SINGULARITY_TMPDIR
```

### Issue: Slow first run

**Expected behavior:** First run downloads ~2-5 GB of containers  
**Time:** 5-15 minutes depending on network

**Solution:** Be patient, subsequent runs are fast!

---

## HPC-Specific Configuration

### SLURM Cluster Example

```groovy
// conf/hpc_slurm.config
process {
    executor = 'slurm'
    queue    = 'normal'
    
    withLabel: 'process_low' {
        cpus   = 4
        memory = 16.GB
        time   = 2.h
    }
    withLabel: 'process_high' {
        cpus   = 16
        memory = 64.GB
        time   = 8.h
    }
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/scratch/shared/singularity_cache'
}
```

**Usage:**
```bash
nextflow run main.nf \
  --input data.h5ad \
  --outdir results \
  -profile singularity \
  -c conf/hpc_slurm.config
```

---

## Performance Tips

### 1. Use Shared Cache on HPC

```bash
# Admin creates shared cache
sudo mkdir -p /scratch/shared/singularity_cache
sudo chmod 1777 /scratch/shared/singularity_cache

# Users configure
export SINGULARITY_CACHEDIR=/scratch/shared/singularity_cache
```

**Benefit:** Download containers once, all users benefit

### 2. Pre-Pull Containers Before Large Runs

```bash
# Pull all scAnnex containers before analysis
bash scripts/prepull_containers.sh

# Then run pipeline (no download delays)
nextflow run main.nf -profile singularity ...
```

### 3. Enable Image Compression

```groovy
// nextflow.config
singularity {
    enabled = true
    autoMounts = true
    runOptions = '--no-home --cleanenv'
}
```

---

## Comparison: Singularity vs Docker vs Conda

| Feature | Singularity | Docker | Conda |
|---------|-------------|--------|-------|
| **HPC Support** | ✅ Excellent | ❌ Poor | ✅ Good |
| **Root Required** | ❌ No | ⚠️ Yes (setup) | ❌ No |
| **Reproducibility** | ✅ Perfect | ✅ Perfect | ⚠️ Good |
| **Speed** | ✅ Fast | ✅ Fast | ❌ Slower |
| **Disk Usage** | 2-5 GB | 3-6 GB | 5-10 GB |
| **First Run** | 10 min | 5 min | 20-30 min |
| **Subsequent** | Fast | Fast | Fast |
| **Portability** | ✅ Excellent | ⚠️ Limited | ✅ Good |

**Recommendation:** Use Singularity for all production work

---

## Frequently Asked Questions

### Q: Singularity vs Apptainer - which to use?

**A:** They're the same! Apptainer is the new name (since 2021).  
Use whichever is installed on your system.

### Q: Can I use pre-built Singularity images?

**A:** Yes! scAnnex containers are automatically converted from Docker images.  
No manual building required.

### Q: Do I need to download containers manually?

**A:** No! Nextflow handles everything automatically.  
Just run the pipeline with `-profile singularity`

### Q: Can I share containers between projects?

**A:** Yes! Cache directory is shared by default.  
All Nextflow pipelines using Singularity benefit from the same cache.

### Q: What if my HPC doesn't have Singularity?

**A:** Request installation from your HPC admin (it's standard).  
Or use `-profile conda` (slower but works everywhere)

---

## Further Reading

- **Apptainer Docs:** https://apptainer.org/docs/
- **Singularity Docs:** https://sylabs.io/docs/
- **Nextflow Containers:** https://nextflow.io/docs/latest/container.html
- **Seqera Best Practices:** https://seqera.io/blog/

---

**Last Updated:** 2026-01-21  
**scAnnex Version:** 1.0 post-fixes
