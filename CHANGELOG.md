# Changelog

All notable changes to scAnnex will be documented in this file.

## [Unreleased]

### Changed
- Updated all modules to use modern Wave/Seqera containers for Singularity profile
- Replaced old biocontainers (scanpy 1.7.2, anndata 0.8.x) with Wave containers (scanpy 1.12, anndata 0.12.6)
- UNIFY_INPUT: Now uses `pip_scanpy:46ad0720691ef95a` (scanpy 1.12)
- QUALITY_CONTROL: Now uses `pip_scanpy:46ad0720691ef95a` (scanpy 1.12)
- DOUBLET_DETECTION: Now uses `scrublet_scanpy:31ec5c9d8d3579e3` (scanpy 1.12 + scrublet 0.2.3)
- STANDARD_PROCESSING: Now uses `leidenalg_python-igraph_scanpy:b3d23ac8b00c1980` (scanpy 1.12 + leidenalg 0.11 + igraph 1.0)
- NORMALIZE_INTEGRATE: Now uses `harmonypy_scanpy:a8efcfdf23c8acc8` (scanpy 1.12 + harmonypy 0.2.0)
- AUTO_ANNOT_CELLTYPIST: Now uses `celltypist_scanpy:20c2e982b26fecc1` (scanpy 1.12 + celltypist 1.7.1)

### Fixed
- Fixed AUTO_ANNOT_CELLTYPIST anndata incompatibility (nullable-string format issue)
- Old celltypist biocontainer couldn't read H5AD files from modern anndata, now resolved with Wave container
- All Python scripts now handle anndata compatibility (nullable strings) defensively
- Fixed container compatibility issues with modern anndata versions

### Planned for v1.1.0
- Interactive Shiny dashboard for manual curation
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
