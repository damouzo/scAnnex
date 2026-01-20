# scAnnex - Camino a Production Ready

**Objetivo:** Que cualquier persona pueda `git clone` → correr sus datos → ver en dashboard

**Estado actual:** ~80% completo  
**Tiempo estimado para completar:** 2-4 horas de trabajo enfocado

---

## 🚨 OBSTÁCULOS CRÍTICOS (Debes Resolver)

### 1. TEST END-TO-END COMPLETO 🔴 (MÁS IMPORTANTE)

**¿Por qué es crítico?**
- No sabes si el pipeline funciona de principio a fin
- Un usuario que clone el repo podría encontrarse con errores fatales
- El dashboard podría no cargar los resultados

**¿Qué hacer?**

#### Test 1: Con datos incluidos (30 min)
```bash
cd /home/damo/scAnnex

# Limpiar resultados previos
rm -rf test_results

# Correr pipeline completo con test profile
nextflow run main.nf \
  -profile test,docker \
  --outdir test_results \
  -resume

# Verificar que todo se generó
ls -lh test_results/
```

**Checklist de outputs esperados:**
```
test_results/
├── unified_input/
│   └── PBMC_MTX_quick_test.h5ad
├── qc/
│   ├── PBMC_MTX_quick_test_qc.h5ad
│   ├── qc_report.json
│   └── plots/*.png
├── doublet/  (si run_doublet_detection = true)
├── standard_processing/
│   └── PBMC_MTX_quick_test_processed.h5ad
├── integration/
│   └── PBMC_MTX_quick_test_integrated.h5ad
├── auto/  (si run_auto_annotation = true)
│   └── PBMC_MTX_quick_test_annotated.h5ad
└── pipeline_info/
    ├── execution_report.html
    └── execution_timeline.html
```

**Si falla:**
- Revisar `test_results/.nextflow.log`
- Revisar `test_results/pipeline_info/execution_trace.txt`
- Identificar qué proceso falló
- Debuggear ese proceso específico

#### Test 2: Dashboard con resultados del test (10 min)
```bash
cd dashboard
conda activate scannex-dashboard
./launch_dashboard.sh

# En el dashboard:
# 1. Cargar: ../test_results/auto/PBMC_MTX_quick_test_annotated.h5ad
# 2. Verificar que UMAP se muestra
# 3. Verificar que cell types aparecen
# 4. Probar buscar un gen (ej: CD3D)
```

**Si dashboard falla:**
- Revisar que el h5ad tiene `.obsm['X_umap']`
- Revisar que tiene `predicted_labels` (si corriste annotation)
- Verificar logs en la consola de R

#### Test 3: Con datos externos REALES (1 hora)

**Descargar PBMC 3k de 10X:**
```bash
mkdir -p external_test
cd external_test

# Descargar datos públicos (3k PBMCs)
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz

cd ..
```

**Crear samplesheet:**
```bash
cat > external_test/samplesheet.csv << 'EOF'
sample_id,data_path,input_type,batch
PBMC_3k,external_test/filtered_gene_bc_matrices/hg19,mtx,batch1
EOF
```

**Correr pipeline:**
```bash
nextflow run main.nf \
  --input external_test/samplesheet.csv \
  --outdir external_test_results \
  --run_auto_annotation true \
  --annotation_celltypist_model Immune_All_Low.pkl \
  -profile docker
```

**Este es EL TEST DEFINITIVO** - si esto funciona, tu pipeline está listo.

---

### 2. VERIFICAR CONTAINERS PÚBLICOS 🔴

**¿Por qué es crítico?**
- Si un container no existe, el pipeline FALLA para cualquier usuario
- Tú estás usando conda localmente, pero otros usuarios usarán containers

**¿Qué hacer?**

```bash
cd /home/damo/scAnnex

# Verificar manualmente cada container
docker pull quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0
docker pull docker.io/satijalab/seurat:5.0.0
docker pull quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0

# O usar el script
python3 bin/check_containers.py
```

**Si alguno NO existe:**

**Opción A (Rápido):** Cambiar a una versión disponible
```bash
# Buscar versiones disponibles
# Para scanpy:
# https://quay.io/repository/biocontainers/scanpy?tab=tags

# Actualizar en cada módulo .nf
```

**Opción B (Mejor):** Construir y publicar tus propios containers
```bash
# Crear Dockerfile optimizado
# Pushear a tu Docker Hub o GitHub Container Registry
# Actualizar referencias en módulos
```

**Opción C (Fallback):** Usar solo conda
```bash
# En nextflow.config, hacer que docker profile use conda
# Más lento pero garantizado
```

---

### 3. DOCUMENTACIÓN "GETTING STARTED" 🟡

**¿Por qué es importante?**
- Primera impresión del usuario
- Si no es clara, abandonan el proyecto

**¿Qué necesitas?**

Crear **UN SOLO ARCHIVO** que el usuario lea primero:

<function_calls>
<invoke name="write">
<parameter name="content"># scAnnex - Getting Started

**Single-cell RNA-seq analysis pipeline from raw data to interactive dashboard in 3 steps.**

---

## Quick Start (5 minutes)

### Prerequisites
- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- [Docker](https://docs.docker.com/get-docker/) or [Conda](https://docs.conda.io/en/latest/miniconda.html)

### Test the Pipeline

```bash
# 1. Clone repository
git clone https://github.com/YOUR-USERNAME/scAnnex.git
cd scAnnex

# 2. Run test (1k cells, ~5 minutes)
nextflow run main.nf -profile test,docker --outdir test_results

# 3. View results in dashboard
cd dashboard
./setup_dashboard.sh
./launch_dashboard.sh
# Open http://localhost:8888
# Load: ../test_results/auto/SAMPLE_annotated.h5ad
```

**If it works → Your setup is correct! Continue to "Analyze Your Data" below.**

---

## Analyze Your Data

### Step 1: Prepare Input

#### Option A: H5AD file (simplest)
```bash
# If you have a .h5ad file with raw counts:
nextflow run main.nf \
  --input your_data.h5ad \
  --input_type h5ad \
  --outdir results \
  -profile docker
```

#### Option B: 10X MTX format
```bash
# If you have filtered_feature_bc_matrix/ from cellranger:
nextflow run main.nf \
  --input path/to/filtered_feature_bc_matrix \
  --input_type mtx \
  --outdir results \
  -profile docker
```

#### Option C: Multiple samples (samplesheet)
```bash
# Create samplesheet.csv:
cat > samplesheet.csv << 'EOF'
sample_id,data_path,input_type,batch
sample1,/path/to/sample1.h5ad,h5ad,batch1
sample2,/path/to/sample2.h5ad,h5ad,batch2
EOF

# Run pipeline:
nextflow run main.nf \
  --input samplesheet.csv \
  --run_integration true \
  --batch_key batch \
  --outdir results \
  -profile docker
```

### Step 2: Wait for Results

**Expected time:**
- 1k cells: ~5 minutes
- 10k cells: ~15 minutes
- 100k cells: ~1-2 hours

**Monitor progress:**
```bash
tail -f .nextflow.log
```

### Step 3: Explore in Dashboard

```bash
cd dashboard
./launch_dashboard.sh

# In browser:
# 1. Load: ../results/auto/SAMPLE_annotated.h5ad
# 2. Explore UMAP colored by cell types
# 3. Search genes (CD3D, CD14, MS4A1, etc.)
# 4. Export plots
```

---

## Pipeline Features

**What the pipeline does automatically:**

1. **Quality Control** (MAD-based adaptive filtering)
2. **Doublet Detection** (Scrublet)
3. **Normalization** (log1p, 10k target sum)
4. **Feature Selection** (2k highly variable genes)
5. **Dimensionality Reduction** (PCA + UMAP)
6. **Clustering** (Leiden, multiple resolutions)
7. **Cell Type Annotation** (CellTypist with pre-trained models)
8. **Batch Correction** (optional, Harmony/BBKNN)

**Output files:**
```
results/
├── qc/
│   ├── SAMPLE_qc.h5ad                # After QC filtering
│   └── qc_report.json                # QC metrics
├── auto/
│   └── SAMPLE_annotated.h5ad         # Final annotated data
└── pipeline_info/
    ├── execution_report.html         # Run statistics
    └── execution_timeline.html       # Timeline visualization
```

---

## Common Issues

### Pipeline fails with "Container not found"

**Problem:** Docker can't pull containers

**Solution 1:** Use Singularity instead
```bash
nextflow run main.nf -profile test,singularity --outdir results
```

**Solution 2:** Use Conda (slower but always works)
```bash
nextflow run main.nf -profile test,conda --outdir results
```

### Dashboard won't load data

**Problem:** File path incorrect or backed mode issue

**Solution:**
```bash
# Use absolute path
/home/user/scAnnex/results/auto/SAMPLE_annotated.h5ad

# Disable backed mode for small files (<50k cells)
# In dashboard: uncheck "Use backed mode"
```

### Out of memory error

**Problem:** Dataset too large for available RAM

**Solution:** Use low-memory profile
```bash
nextflow run main.nf -profile docker,low_memory --outdir results
```

---

## Configuration

### Key Parameters

Edit `nextflow.config` or use command-line flags:

```bash
# QC thresholds
--min_genes 200                     # Minimum genes per cell
--min_counts 500                    # Minimum UMI per cell
--max_mito_percent 20               # Max % mitochondrial genes

# Clustering
--clustering_resolutions '0.1,0.3,0.5,0.7,0.9'

# Annotation
--run_auto_annotation true
--annotation_celltypist_model Immune_All_Low.pkl

# Batch correction
--run_integration true
--batch_key batch
--integration_method harmony
```

### Available CellTypist Models

- `Immune_All_Low.pkl` - Human immune cells (recommended)
- `Immune_All_High.pkl` - Human immune cells (higher resolution)
- `Pan_Tissue.pkl` - General human tissues
- Custom models: place in `assets/celltypist_models/`

---

## Need Help?

- **Documentation:** [docs/](docs/)
- **Dashboard Guide:** [docs/dashboard/README.md](docs/dashboard/README.md)
- **Issues:** https://github.com/YOUR-USERNAME/scAnnex/issues
- **Examples:** [docs/examples/](docs/examples/)

---

## Citation

If you use scAnnex in your research:

```bibtex
@software{scannex2026,
  author = {Your Name},
  title = {scAnnex: Automated Single-Cell RNA-seq Analysis Pipeline},
  year = {2026},
  url = {https://github.com/YOUR-USERNAME/scAnnex}
}
```

---

## Next Steps

- [Dashboard Documentation](docs/dashboard/README.md) - Complete dashboard guide
- [Pipeline Details](docs/pipeline.md) - Technical documentation
- [Examples](docs/examples/) - Analysis walkthroughs
- [Troubleshooting](docs/troubleshooting.md) - Common issues

---

**Quick Links:**
- [Test Profile Configuration](conf/test.config)
- [Parameter Reference](docs/parameters.md)
- [HPC Usage](docs/hpc.md)
- [Development Guide](docs/development.md)
