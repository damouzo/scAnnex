# Changelog

All notable changes to scAnnex will be documented in this file.

## [1.0.0] - 2025-01-XX

### Added
- Initial release of scAnnex pipeline
- Multi-format input support (H5AD, RDS, 10X MTX)
- Quality control with automated filtering
- Doublet detection using Scrublet
- Normalization with multiple methods (log, scran)
- Batch integration (Harmony, Scanorama, BBKNN)
- Dimensionality reduction (PCA, UMAP)
- Leiden clustering with configurable resolution
- Automated cell type annotation
- Differential expression analysis (Wilcoxon, t-test, logistic regression)
- Trajectory inference (PAGA, DPT)
- HTML report generation
- Docker, Singularity, and Conda profiles
- Comprehensive test profile with example data
- Full DSL2 Nextflow implementation
- Modular process organization
- Extensive parameter validation
- Resource optimization with process labels

### Features
- ✅ Strict Nextflow DSL2 syntax (v25.10+ compatible)
- ✅ All files pass `nextflow lint` without errors
- ✅ Modular design with reusable processes
- ✅ Comprehensive error handling
- ✅ Automated retry strategies
- ✅ Progress tracking and logging
- ✅ MultiQC-style HTML reports

### Pipeline Steps
1. **UNIFY_INPUT**: Convert various formats to standardized H5AD
2. **QUALITY_CONTROL**: Filter cells and genes based on QC metrics
3. **DOUBLET_DETECTION**: Identify and annotate doublets
4. **NORMALIZE_INTEGRATE**: Normalize and optionally integrate batches
5. **DIMENSIONALITY_REDUCTION**: Compute PCA and UMAP
6. **CLUSTER_ANNOTATE**: Perform clustering and annotation
7. **DIFFERENTIAL_EXPRESSION**: Identify marker genes
8. **TRAJECTORY_ANALYSIS** (optional): Infer developmental trajectories
9. **GENERATE_REPORT**: Create comprehensive HTML report

### Dependencies
- Nextflow ≥23.04.0
- Python packages: scanpy, anndata, scrublet, harmonypy, scanorama, bbknn
- R (optional): Seurat, SingleCellExperiment (for RDS input)
- Container engines: Docker, Singularity, or Conda

### Known Limitations
- Large datasets (>500k cells) may require significant memory
- RDS conversion requires R and rpy2 packages
- Interactive Shiny dashboard not yet implemented (planned for v1.1.0)

## [Unreleased]

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
