# Changelog

All notable changes to scAnnex will be documented in this file.

## [Unreleased]

### Changed
- **Container system optimization**: Cleaned up legacy profiles and improved container management
  - Removed experimental/unsupported profiles: `wave`, `apptainer`, `podman`, `docker`
  - Maintained production profiles: `singularity`, `conda`, `slurm`, `apocrita`
  - GSEA container migrated to DockerHub: `docker://damouzo/scannex:gsea-1.0`
  - Singularity now auto-downloads GSEA container from DockerHub (no manual build required)
  - Moved conda package cache to scratch in Apocrita (`CONDA_PKGS_DIRS` in scratch to avoid filling `/data/home`)
- **Dashboard environment setup optimization (HPC)**:
  - Hybrid auto-creation with mamba (10-20x faster than conda)
  - **Configurable location**: New param `--dashboard_conda_dir` to override default location
  - Environment variable `SCANNEX_DASHBOARD_CONDA_DIR` for both HPC and local launchers
  - Smart defaults: HPC scratch (`/gpfs/scratch/$USER/conda_envs/`) or home (`~/.conda/envs/`)
  - **Age warning**: Alerts when environment is >60 days old (before 65-day auto-cleanup)
  - 10-minute timeout for environment creation with fallback to manual instructions
  - Automatic `module load miniforge` detection on HPC clusters
  - Improved error handling: detects missing `timeout` command, permission issues
  - Non-critical failures skip dashboard without stopping pipeline
  - Reduced setup time from 30-60 minutes to 2-5 minutes
  - Fixed permissions issue with system-wide miniforge installation
  - Applied configurable location to both `launch_dashboard.sh` (local) and `launch_dashboard_hpc.sh` (HPC)
- Updated all modules to use modern Wave/Seqera containers for Singularity profile
- Replaced old biocontainers (scanpy 1.7.2, anndata 0.8.x) with Wave containers (scanpy 1.12, anndata 0.12.6)
- UNIFY_INPUT: Now uses `pip_scanpy:46ad0720691ef95a` (scanpy 1.12)
- QUALITY_CONTROL: Now uses `pip_scanpy:46ad0720691ef95a` (scanpy 1.12)
- DOUBLET_DETECTION: Now uses `scrublet_scanpy:31ec5c9d8d3579e3` (scanpy 1.12 + scrublet 0.2.3)
- STANDARD_PROCESSING: Now uses `leidenalg_python-igraph_scanpy:b3d23ac8b00c1980` (scanpy 1.12 + leidenalg 0.11 + igraph 1.0)
- NORMALIZE_INTEGRATE: Now uses `harmonypy_scanpy:a8efcfdf23c8acc8` (scanpy 1.12 + harmonypy 0.2.0)
- AUTO_ANNOT_CELLTYPIST: Now uses `celltypist_scanpy:20c2e982b26fecc1` (scanpy 1.12 + celltypist 1.7.1)
- AUTO_ANNOT_SUMMARIZE: Updated conda directive to include scipy for consistency with Wave container

### Added
- Script for building and pushing GSEA container to DockerHub (`scripts/build_and_push_gsea_container.sh`)
- Script for testing GSEA container download from DockerHub (`scripts/test_gsea_container.sh`)
- Comprehensive container documentation (`Info_containers_Status.md`) with nf-core migration guide

### Removed
- Legacy Dockerfile (`containers/Dockerfile.scanpy`) - no longer used
- Unsupported execution profiles from `nextflow.config`

### Fixed
- Fixed AUTO_ANNOT_CELLTYPIST anndata incompatibility (nullable-string format issue)
- Old celltypist biocontainer couldn't read H5AD files from modern anndata, now resolved with Wave container
- All Python scripts now handle anndata compatibility (nullable strings) defensively
- Fixed container compatibility issues with modern anndata versions

### Planned for v1.1.0
- Cell-cell communication analysis (CellPhoneDB, NicheNet)
- Gene regulatory network inference
- Spatial transcriptomics support
- Integration with Seurat v5
- Support for multimodal data (CITE-seq, ATAC-seq)
- Automated parameter optimization
- Cloud execution profiles (AWS, GCP, Azure)

### Planned for v1.2.0
- Real-time analysis monitoring dashboard
- Advanced trajectory methods (Monocle3, Slingshot)
- Gene set enrichment analysis
- Automated figure generation for publications
- Benchmarking module for method comparison
