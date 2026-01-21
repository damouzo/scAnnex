# scAnnex Execution Profiles

This document describes the available execution profiles for running scAnnex and provides guidance on choosing the right profile for your environment.

## Quick Start

Choose a profile when running the pipeline:

```bash
nextflow run main.nf --input <your_data> -profile <profile_name>
```

## Available Profiles

### 1. Conda Profile (Recommended for Testing)

**Usage:** `-profile conda`

**Best for:**
- Testing and development
- Environments without container support
- Systems where you have conda installed

**Pros:**
- Works on any system with conda/mamba
- No special permissions required
- Easy to customize environments

**Cons:**
- Slower than containers (environment creation overhead)
- Larger disk space usage
- Dependency conflicts possible

**Example:**
```bash
nextflow run main.nf \
    --input test_data/outputs/PBMC_MTX_quick_test.h5ad \
    -profile conda
```

### 2. Wave Profile (Recommended for Modern Workflows)

**Usage:** `-profile wave`

**Best for:**
- Modern cloud-native workflows
- CI/CD pipelines
- Users who want reproducibility without pre-built containers

**Pros:**
- Builds containers on-demand from conda specs
- Maximum reproducibility
- No need to maintain container registries
- Combines conda flexibility with container isolation

**Cons:**
- Requires internet connection
- First run is slower (container building)
- Requires Seqera Wave service access

**Example:**
```bash
nextflow run main.nf \
    --input your_data.h5ad \
    -profile wave
```

**Requirements:**
- A container runtime (Docker, Singularity, or Podman)
- Internet connection (to build containers)

**Important:** Wave builds containers on-demand from conda specifications, then executes them. This means:
- You **must** have a container runtime installed (Docker, Singularity, etc.)
- First run is slower (builds containers)
- Subsequent runs are fast (containers are cached)

**WSL2 Note:** Wave requires Docker to be installed. If you get errors like "Failed to create Conda environment", it means Wave is falling back to local conda because no container runtime is available. Either:
- Install Docker Desktop for Windows with WSL2 integration, OR
- Use `-profile conda` instead (works without Docker)

**Wave caching:** Containers are built once and cached. First run ~10-15 min, subsequent runs ~5-7 min.

### 3. Singularity/Apptainer Profile (Recommended for HPC)

**Usage:** `-profile singularity` or `-profile apptainer`

**Best for:**
- HPC clusters and shared computing environments
- Production pipelines
- Environments where Docker is not available

**Pros:**
- Standard on HPC systems
- Good performance
- No root privileges required
- Better security than Docker

**Cons:**
- Requires Singularity/Apptainer installation
- Container images must exist in registry
- May have issues with overlay filesystems in WSL2

**Example:**
```bash
nextflow run main.nf \
    --input your_data.h5ad \
    -profile singularity
```

**Cache location:** `~/.singularity/cache` or `~/.apptainer/cache`

### 4. Docker Profile (For Local Development)

**Usage:** `-profile docker`

**Best for:**
- Local development on macOS or Linux
- Environments with Docker installed
- Users familiar with Docker ecosystem

**Pros:**
- Fast execution
- Good isolation
- Popular and well-documented

**Cons:**
- Requires Docker daemon
- Needs root or docker group membership
- Not available on most HPC systems
- May not work in WSL2 without Docker Desktop

**Example:**
```bash
nextflow run main.nf \
    --input your_data.h5ad \
    -profile docker
```

## Comparison Table

| Profile | Speed | Setup Difficulty | HPC Compatible | Reproducibility | Disk Usage |
|---------|-------|------------------|----------------|-----------------|------------|
| Conda | Slow | Easy | Yes | Good | High |
| Wave | Medium-Fast | Easy | Yes | Excellent | Medium |
| Singularity | Fast | Medium | Yes | Excellent | Medium |
| Docker | Fast | Easy | Limited | Excellent | Medium |

## Default Behavior

If no profile is specified, Nextflow will attempt to use conda environments by default. We recommend explicitly specifying a profile for better control and reproducibility.

## Profile-Specific Configuration

### Conda Profile
- Uses channels: `conda-forge`, `bioconda`, `defaults`
- Creates environments in `work/conda/`
- Can be customized via `params.conda_channels`

### Wave Profile
- Automatically builds containers from conda specifications
- Uses Seqera Wave service (https://seqera.io/wave/)
- Caches built containers for reuse
- Strategy: `conda,container` (tries conda first, then containers)

### Singularity/Apptainer Profiles
- Auto-mounts home directory
- Cache directory: `~/.singularity/cache` or `~/.apptainer/cache`
- Pull timeout: 20 minutes (configurable)
- Uses docker:// protocol to pull from Docker registries

### Docker Profile
- Runs as current user: `-u $(id -u):$(id -g)`
- Mounts required directories automatically
- Uses Docker daemon socket

## Troubleshooting

### Conda Profile Issues

**Problem:** Conda environment creation is slow
```bash
# Solution: Use mamba instead
conda install -c conda-forge mamba
# Nextflow will automatically use mamba if available
```

**Problem:** Dependency conflicts
```bash
# Solution: Clean conda cache
conda clean --all
rm -rf work/conda/
```

### Wave Profile Issues

**Problem:** Wave container building times out
```bash
# Solution: Check internet connection and Wave service status
# Increase timeout in nextflow.config if needed
```

**Problem:** Authentication error
```bash
# Solution: Wave is a free service, but may require registration for heavy use
# Visit: https://seqera.io/wave/
```

### Singularity Profile Issues

**Problem:** Image pull fails (404 error)
```bash
# This means the container tag doesn't exist
# Solution 1: Use Wave profile instead (builds from conda specs)
nextflow run main.nf --input data.h5ad -profile wave

# Solution 2: Use Conda profile
nextflow run main.nf --input data.h5ad -profile conda
```

**Problem:** Permission denied / overlay filesystem issues (WSL2)
```bash
# WSL2 may have issues with Singularity overlays
# Solution: Use conda or Wave profile instead
```

### Docker Profile Issues

**Problem:** Docker daemon not running
```bash
# Solution: Start Docker daemon
sudo systemctl start docker  # Linux
# or launch Docker Desktop (macOS/Windows)
```

**Problem:** Permission denied
```bash
# Solution: Add user to docker group
sudo usermod -aG docker $USER
# Log out and back in
```

### Dashboard Compatibility Issues

**Problem:** Dashboard fails to load H5AD files with `IORegistryError`
```
Error: anndata._io.specs.registry.IORegistryError: No read method registered for IOSpec
```

**Cause:** Version mismatch between pipeline and dashboard anndata versions.

**Solution:** Recreate dashboard environment with updated dependencies
```bash
cd dashboard
conda env remove -n scannex-dashboard -y
conda env create -f environment_dashboard.yml
```

The dashboard now requires:
- **Python 3.11+** (was 3.10)
- **anndata >=0.12.0** (was >=0.10.0)

This ensures compatibility with H5AD files generated by all execution profiles (Wave, Conda, Singularity, Docker).

## Recommendations by Environment

### Local Development (Linux/macOS)
1. **Docker** - Fast and convenient
2. **Conda** - If Docker not available
3. **Wave** - For testing production setup

### HPC Cluster
1. **Singularity/Apptainer** - Standard on most clusters
2. **Wave** - Modern alternative with excellent reproducibility
3. **Conda** - Fallback if containers unavailable

### Cloud (AWS, GCP, Azure)
1. **Wave** - Best for cloud-native workflows
2. **Docker** - If using Docker-based compute
3. **Singularity** - For traditional HPC-style cloud

### CI/CD Pipelines
1. **Wave** - Reproducible and cacheable
2. **Docker** - Fast for quick tests
3. **Conda** - For environment testing

### WSL2 (Windows)
1. **Conda** - Most reliable
2. **Wave** - Good alternative
3. **Docker** - Requires Docker Desktop
4. **Singularity** - May have overlay filesystem issues

## Multiple Profiles

You can combine profiles:

```bash
# Use test config with conda
nextflow run main.nf -profile test,conda

# Use laptop config with wave
nextflow run main.nf -profile laptop,wave
```

## Custom Profile

Create a custom profile in your local `nextflow.config`:

```groovy
profiles {
    myprofile {
        process {
            executor = 'slurm'
            queue = 'compute'
            memory = '32 GB'
            cpus = 8
        }
        singularity.enabled = true
        singularity.cacheDir = '/scratch/singularity_cache'
    }
}
```

Then use: `nextflow run main.nf -profile myprofile`

## Environment Variables

Relevant environment variables:

```bash
# Conda
export CONDA_PKGS_DIRS=/path/to/cache  # Conda package cache

# Singularity
export SINGULARITY_CACHEDIR=/path/to/cache  # Image cache
export SINGULARITY_TMPDIR=/tmp  # Temp directory

# Docker
export DOCKER_HOST=unix:///var/run/docker.sock  # Docker daemon

# Wave
export TOWER_ACCESS_TOKEN=your_token  # For authenticated Wave access
```

## Performance Tips

1. **Use SSD storage** for `work/` directory
2. **Set appropriate cache directories** for containers
3. **Use `-resume`** to restart failed runs
4. **Clean up with** `nextflow clean -f` between major runs
5. **Pre-pull containers** for singularity:
   ```bash
   singularity pull docker://quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0
   ```

## More Information

- Nextflow profiles: https://www.nextflow.io/docs/latest/config.html#config-profiles
- Conda environments: https://docs.conda.io/
- Wave containers: https://seqera.io/wave/
- Singularity: https://sylabs.io/docs/
- Docker: https://docs.docker.com/

## Support

If you encounter issues with any profile:
1. Check this documentation
2. Check `.nextflow.log` for detailed error messages
3. Run with `-with-trace -with-report` for diagnostics
4. Open an issue on GitHub with error details and profile used
