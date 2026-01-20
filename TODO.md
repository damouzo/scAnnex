# scAnnex SLC (Single-Cell LifeCycle) - TODO List

**Vision:** Build a Simple, Lovable, Complete single-cell RNA-seq analysis pipeline that takes users from raw counts to interactive dashboard with best-practice methods built-in.

**Last Updated:** January 20, 2026  
**Version:** SLC v1.0

---

## ✅ COMPLETED - Core SLC Pipeline (v1.0)

### Input & Quality Control
- [x] Unified input module (H5, MTX, RDS → H5AD conversion)
- [x] Quantile-based QC filtering (10th-90th percentile)
- [x] MAD-based filtering (fallback option)
- [x] Cell Attrition Log (CSV + TXT reports)
- [x] Before/after QC visualization
- [x] Configurable QC thresholds via nextflow.config

### Doublet Detection
- [x] Scrublet integration
- [x] Toggle for doublet removal (--doublet-removal flag)
- [x] Attrition tracking (JSON output)
- [x] Configurable expected doublet rate

### Standard Processing
- [x] Normalization (target sum 10,000)
- [x] Log transformation
- [x] HVG selection (top 2,000 genes)
- [x] Scaling (max value 10)
- [x] PCA (50 components)
- [x] Neighbor graph (k=15)
- [x] UMAP generation
- [x] Multi-resolution Leiden clustering (5 resolutions: 0.1, 0.3, 0.5, 0.7, 0.9)
- [x] Dashboard-ready outputs (umap_coordinates.csv, cell_metadata.csv)

### Auto-Annotation
- [x] CellTypist integration
- [x] Configurable model selection (default: Immune_All_Low.pkl)
- [x] Majority voting option
- [x] Confidence scores export
- [x] CSV output for downstream analysis

### Integration (Batch Correction)
- [x] Harmony integration method
- [x] BBKNN integration method
- [x] Optional execution (conditional on batch_key)
- [x] Runs after annotation (pipeline order optimized)

### Configuration & Documentation
- [x] Complete nextflow.config with SLC parameters
- [x] Low memory profile (conf/low_memory.config)
- [x] Base configuration (conf/base.config)
- [x] Module configuration (conf/modules.config)
- [x] Docker container documentation (docker/README.md)
- [x] Container tag resolution docs (docs/CONTAINER_TAG_RESOLUTION.md)
- [x] Low memory usage guide (docs/LOW_MEMORY_USAGE.md)
- [x] SLC Quick Start guide (SLC_QUICKSTART.md)
- [x] Complete implementation summary (docs/summary_2026-01-20.md)
- [x] Workflow sync documentation (docs/workflow_sync_2026-01-20.md)

### Docker & Containers
- [x] Scanpy extended container (Dockerfile.scanpy-extended)
- [x] Base scanpy container (Dockerfile.scanpy)
- [x] Build scripts (build-scanpy-extended.sh)
- [x] check_max function fix for low memory systems

### Test Infrastructure
- [x] Test data download script (bin/download_test_data.py)
- [x] PBMC 1k v3 dataset support
- [x] Automatic samplesheet generation
- [x] Batch1/batch2 test setup for integration testing

---

## 🚧 IN PROGRESS - Current Sprint

### Pipeline Testing & Validation
- [ ] **END-TO-END TEST:** Run complete pipeline on PBMC 1k dataset
  - [ ] Download test data successfully
  - [ ] Execute full workflow without errors
  - [ ] Verify all outputs are generated
  - [ ] Check Cell Attrition Log accuracy
  - [ ] Validate multi-resolution clustering results
  - [ ] Confirm CellTypist annotations are correct

### Integration Module Enhancements
- [ ] Add diagnostic UMAPs (pre/post batch correction)
- [ ] Implement batch quality metrics (kBET)
- [ ] Implement LISI (Local Inverse Simpson's Index)
- [ ] Add batch mixing visualization
- [ ] Generate integration report (JSON + plots)

---

## 📋 HIGH PRIORITY - Next Phase

### Dashboard Implementation (SLC v1.1)
- [ ] **Tab 1: Samples & QC**
  - [ ] Sample selector dropdown
  - [ ] Pre-filter vs Post-filter metrics comparison
  - [ ] Interactive Cell Attrition Table
  - [ ] QC threshold adjustment sliders
  - [ ] Dynamic QC plot updates

- [ ] **Tab 2: Clustering**
  - [ ] Interactive UMAP (plotly)
  - [ ] Resolution selector (0.1 - 0.9)
  - [ ] Color by options:
    - [ ] Leiden clusters (all resolutions)
    - [ ] CellTypist annotations
    - [ ] Batch
    - [ ] QC metrics (MT%, nCount, nFeature)
  - [ ] Cluster summary statistics table

- [ ] **Tab 3: Expression**
  - [ ] Gene search functionality
  - [ ] Feature plots (gene expression on UMAP)
  - [ ] Dot plots for marker genes
  - [ ] Top marker genes per cluster table

- [ ] **Dashboard Infrastructure**
  - [ ] Shiny app structure (ui.R, server.R, global.R)
  - [ ] Data loading from pipeline outputs
  - [ ] Performance optimization for large datasets
  - [ ] Export functionality (plots as PNG/PDF)

### Pipeline Robustness
- [ ] Add comprehensive error handling
- [ ] Implement checkpointing (resume failed runs)
- [ ] Add input validation (samplesheet format checking)
- [ ] Memory profiling for all modules
- [ ] Optimize for datasets > 100k cells

---

## 📊 MEDIUM PRIORITY - Feature Expansion

### Marker Gene Detection
- [ ] Implement Wilcoxon rank-sum test for marker genes
- [ ] Add t-test option
- [ ] Export top N markers per cluster
- [ ] Generate heatmaps of top markers
- [ ] CSV output for downstream enrichment analysis

### Differential Expression Analysis
- [ ] DESeq2 integration for pseudo-bulk DE
- [ ] MAST implementation for single-cell DE
- [ ] Volcano plots
- [ ] MA plots
- [ ] Export DE results to CSV

### Batch Correction Quality Metrics
- [ ] kBET (k-nearest neighbor Batch Effect Test)
- [ ] LISI (Local Inverse Simpson's Index)
- [ ] Silhouette scores
- [ ] Mixing metrics
- [ ] Generate comprehensive integration report

### Advanced Annotation Options
- [ ] SingleR integration (alternative to CellTypist)
- [ ] scAnnotate support
- [ ] Manual annotation refinement interface
- [ ] Hierarchical annotation (coarse → fine)
- [ ] Custom marker gene database support

### Quality Reports
- [ ] MultiQC-style HTML report
- [ ] Per-sample QC summary
- [ ] Clustering quality metrics (silhouette, Davies-Bouldin)
- [ ] Annotation confidence distribution
- [ ] Integration quality assessment

---

## 🔬 LOW PRIORITY - Advanced Features (Future)

### Trajectory Analysis
- [ ] PAGA (Partition-based Graph Abstraction)
- [ ] Monocle3 integration
- [ ] RNA velocity (scVelo)
- [ ] Pseudotime calculation
- [ ] Trajectory plotting

### Cell-Cell Communication
- [ ] CellChat integration
- [ ] NicheNet implementation
- [ ] Ligand-receptor database
- [ ] Communication network visualization
- [ ] Signaling pathway enrichment

### Multi-Modal Integration
- [ ] CITE-seq (RNA + Protein)
- [ ] scATAC-seq integration
- [ ] Spatial transcriptomics support
- [ ] MOFA+ for multi-omics factor analysis

### Performance Optimization
- [ ] Parallel processing optimization
- [ ] GPU acceleration for PCA/UMAP
- [ ] Sparse matrix optimizations
- [ ] Incremental data loading (for > 1M cells)

### Alternative Methods
- [ ] Seurat interoperability
- [ ] SCRAN normalization option
- [ ] Louvain clustering (alternative to Leiden)
- [ ] t-SNE as alternative to UMAP
- [ ] Additional integration methods (Scanorama, scGen)

---

## 🛠️ TECHNICAL DEBT & MAINTENANCE

### Code Quality
- [ ] Add unit tests for Python scripts
- [ ] Add integration tests for workflow
- [ ] Implement CI/CD pipeline (GitHub Actions)
- [ ] Code linting (flake8, black for Python)
- [ ] Nextflow DSL2 best practices audit

### Documentation
- [ ] Add inline code comments to all Python scripts
- [ ] Complete API documentation (Sphinx)
- [ ] Tutorial videos (YouTube)
- [ ] FAQ section
- [ ] Troubleshooting guide expansion

### Container Management
- [ ] Pin all software versions (reproducibility)
- [ ] Minimize container sizes
- [ ] Multi-stage builds
- [ ] Automated container builds (CI)
- [ ] Container registry setup (Docker Hub / Quay.io)

### User Experience
- [ ] Add progress indicators for long-running steps
- [ ] Better error messages with solutions
- [ ] Parameter validation and helpful warnings
- [ ] Example datasets library
- [ ] Interactive parameter selection wizard

---

## 🎯 SUCCESS CRITERIA (SLC v1.0)

**Simple:**
- [x] One-command execution
- [x] Sensible defaults (no parameter tweaking needed)
- [x] Clear documentation

**Lovable:**
- [x] Cell Attrition Log (users love transparency)
- [x] Multi-resolution clustering (exploratory freedom)
- [ ] Interactive dashboard (planned)

**Complete:**
- [x] Raw data → Processed data
- [x] QC → Clustering → Annotation
- [ ] Dashboard integration (planned)

---

## 📦 MIGRATION CHECKLIST (For New 16GB Machine)

### Before Push (Current Machine)
- [x] Review all staged changes
- [x] Update documentation (this TODO file)
- [x] Verify .gitignore excludes large files
- [x] Create comprehensive commit message
- [ ] Push to GitHub

### After Clone (New Machine)
- [ ] Clone repository: `git clone <repo-url>`
- [ ] Install Nextflow: `curl -s https://get.nextflow.io | bash`
- [ ] Install Docker (if needed)
- [ ] Download test data: `python bin/download_test_data.py --output-dir test_data`
- [ ] Run test pipeline: `nextflow run main.nf --input test_data/samplesheet.csv --outdir results/ -profile docker`
- [ ] Verify all outputs
- [ ] Continue dashboard development

---

## 📝 NOTES & DECISIONS

### Architecture Decisions
1. **SLC over MVP:** Chose end-to-end functionality over feature-by-feature development
2. **Quantile-based QC:** More robust than fixed thresholds for diverse datasets
3. **Multi-resolution clustering:** Avoids single-resolution bias
4. **CellTypist default:** Fast, accurate, pre-trained models
5. **Integration after annotation:** Preserves biological signal during batch correction

### Key Parameters
- **QC:** 10th-90th percentile (adjustable)
- **Doublet rate:** 5% expected (configurable)
- **Clustering resolutions:** 0.1, 0.3, 0.5, 0.7, 0.9
- **PCA components:** 50
- **HVGs:** 2,000 genes
- **UMAP neighbors:** 15

### Pipeline Philosophy
- **Transparency:** Cell Attrition Log shows every filtering decision
- **Flexibility:** All key parameters configurable via nextflow.config
- **Reproducibility:** Docker containers, pinned versions (in progress)
- **User-Friendly:** Dashboard-ready outputs, clear documentation

---

## 🙏 ACKNOWLEDGMENTS

**Strategy Inspiration:** SLC (Simple, Lovable, Complete) framework  
**Technical Foundation:** Nextflow, Scanpy, CellTypist  
**Community:** Seqera, nf-core best practices  

---

**Last Review:** January 20, 2026  
**Next Review:** After new machine setup and end-to-end test  
**Maintainer:** scAnnex Development Team

---

## 🚀 GETTING STARTED (New Contributors)

1. Read `SLC_QUICKSTART.md` for 3-step pipeline execution
2. Review `docs/summary_2026-01-20.md` for complete feature overview
3. Check this TODO for current priorities
4. Run test pipeline to familiarize with outputs
5. Pick a task from HIGH PRIORITY section

**Questions?** Check `docs/Troubleshooting.md` or open an issue on GitHub.
