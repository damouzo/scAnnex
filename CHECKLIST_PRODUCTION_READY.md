# CHECKLIST: scAnnex Production Ready

**Fecha:** Enero 20, 2026  
**Objetivo:** Pipeline funcionando para cualquier usuario  
**Tiempo estimado:** 2-4 horas

---

## ✅ YA COMPLETADO (No tocar)

- ✅ Pipeline implementado con todos los módulos
- ✅ Dashboard funcionando con fix de UMAP
- ✅ Documentación organizada en docs/
- ✅ Test data incluido
- ✅ Containers definidos en módulos
- ✅ Scripts de deployment (Conda/Docker/Apptainer)
- ✅ GitHub Actions configurado
- ✅ TODO.md actualizado

---

## 🔴 CRÍTICO - HACER AHORA (Orden de prioridad)

### 1. TEST END-TO-END ⏱️ 30 min

**¿Qué?** Correr el pipeline completo y verificar que funciona

**¿Cómo?**
```bash
cd /home/damo/scAnnex
rm -rf test_e2e_results

# Test con perfil incluido
nextflow run main.nf \
  -profile test,docker \
  --outdir test_e2e_results \
  --run_auto_annotation true \
  -resume

# Checklist:
# [ ] Pipeline termina sin errores
# [ ] Se generan todos los outputs
# [ ] QC plots existen
# [ ] H5AD final tiene X_umap
# [ ] H5AD final tiene predicted_labels (si annotation = true)
```

**Si falla:**
- Revisar `.nextflow.log`
- Identificar proceso que falló
- Corregir el error
- **NO SEGUIR** hasta que esto funcione

**Resultado esperado:**
```
test_e2e_results/
├── qc/
│   ├── *.h5ad
│   └── qc_report.json
├── standard_processing/ (o auto/ si annotation está activo)
│   └── *_annotated.h5ad  ← ESTE ES EL IMPORTANTE
└── pipeline_info/
    └── execution_report.html
```

---

### 2. DASHBOARD CON RESULTADOS DEL TEST ⏱️ 10 min

**¿Qué?** Verificar que el dashboard carga los resultados del pipeline

**¿Cómo?**
```bash
cd dashboard
conda activate scannex-dashboard
./launch_dashboard.sh

# En el navegador (http://localhost:8888):
# 1. En "H5AD File Path", poner ruta ABSOLUTA:
#    /home/damo/scAnnex/test_e2e_results/auto/SAMPLE_annotated.h5ad
#
# 2. Click "Load Data"
#
# 3. Ir a "Clustering & UMAP"
#
# 4. Verificar:
#    [ ] UMAP se muestra
#    [ ] Dropdown "Color by" tiene opciones
#    [ ] Si seleccionas "predicted_labels", se ven cell types
#
# 5. Ir a "Gene Expression"
#
# 6. Buscar gen "CD3D" y click "Plot Expression"
#
# 7. Verificar:
#    [ ] UMAP se muestra con colores de expresión
```

**Si falla:**
- Revisar que el h5ad tiene las keys necesarias:
  ```python
  import anndata as ad
  adata = ad.read_h5ad('test_e2e_results/auto/SAMPLE_annotated.h5ad')
  print(adata.obs.columns)  # Debe tener predicted_labels
  print(adata.obsm.keys())  # Debe tener X_umap
  ```

---

### 3. VERIFICAR CONTAINERS PÚBLICOS ⏱️ 15 min

**¿Qué?** Asegurar que los containers que usa el pipeline están disponibles

**¿Cómo?**
```bash
# Intentar pull de cada container
docker pull quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0
docker pull docker.io/satijalab/seurat:5.0.0
docker pull quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0
```

**Si alguno falla:**

#### Opción A: Buscar versión alternativa
```bash
# Buscar en https://quay.io/repository/biocontainers/scanpy?tab=tags
# Encontrar versión disponible
# Actualizar en TODOS los módulos que lo usan
```

#### Opción B: Cambiar a conda (más seguro)
```bash
# En nextflow.config, línea ~124:
docker {
    docker.enabled = false  # ← Cambiar a false
    conda.enabled = true    # ← Agregar esto
}
```

---

## 🟡 IMPORTANTE - HACER DESPUÉS

### 4. TEST CON DATOS EXTERNOS REALES ⏱️ 1 hora

**¿Por qué?** El test incluido podría tener quirks que oculten bugs

**¿Cómo?**
```bash
# Descargar PBMC 3k de 10X (datos públicos)
mkdir -p external_test && cd external_test
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
cd ..

# Crear samplesheet
cat > external_test/samplesheet.csv << 'EOF'
sample_id,data_path,input_type,batch
PBMC3k,external_test/filtered_gene_bc_matrices/hg19,mtx,batch1
EOF

# Correr pipeline
nextflow run main.nf \
  --input external_test/samplesheet.csv \
  --outdir external_test_results \
  --run_auto_annotation true \
  -profile docker

# Verificar que funciona y dashboard carga los resultados
```

---

### 5. CREAR GETTING_STARTED.md ⏱️ Ya hecho ✅

**Status:** Ya creado arriba, solo necesitas:
1. Reemplazar `YOUR-USERNAME` con tu username de GitHub
2. Revisar que los comandos son correctos
3. Agregarlo al repo

---

### 6. ACTUALIZAR README.md PRINCIPAL ⏱️ 15 min

**¿Qué?** Hacer que el README sea más conciso y apunte a GETTING_STARTED.md

**Cambios sugeridos:**
```markdown
# scAnnex

Single-cell RNA-seq analysis from raw data to interactive dashboard.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Quick Start

```bash
# Test the pipeline (5 minutes)
git clone https://github.com/YOUR-USERNAME/scAnnex.git
cd scAnnex
nextflow run main.nf -profile test,docker --outdir results

# View results
cd dashboard && ./launch_dashboard.sh
# Open http://localhost:8888
```

**👉 [Complete Getting Started Guide](GETTING_STARTED.md)**

## Features

✅ Automated QC with MAD-based filtering  
✅ Doublet detection (Scrublet)  
✅ Batch correction (Harmony/BBKNN)  
✅ Cell type annotation (CellTypist)  
✅ Interactive dashboard (R Shiny)  
✅ Multiple input formats (H5AD, MTX, RDS)

## Documentation

- 📖 [Getting Started](GETTING_STARTED.md) - Complete setup guide
- 📊 [Dashboard Guide](docs/dashboard/README.md) - Interactive visualization
- ⚙️  [Pipeline Details](docs/pipeline.md) - Technical documentation
- 🔧 [Troubleshooting](docs/troubleshooting.md) - Common issues

## Citation

[Your citation info]

## License

MIT
```

---

## 🟢 OPCIONAL - NICE TO HAVE

### 7. CI/CD Testing
- GitHub Actions que corran el test automáticamente en cada push
- Asegura que no rompes el pipeline con cambios futuros

### 8. Ejemplos Adicionales
- `docs/examples/` con análisis completos
- Diferentes tipos de datasets (immune, brain, tumor, etc.)

### 9. Video Tutorial
- Screencast de 5-10 min mostrando el pipeline + dashboard
- Súbelo a YouTube y enlázalo en README

### 10. Publicar Containers
- Construir tus propios containers optimizados
- Publicar en GitHub Container Registry
- Más control sobre dependencias

---

## 📋 CHECKLIST FINAL ANTES DE "RELEASE"

Antes de considerar el proyecto "production-ready", verifica:

### Pipeline
- [ ] Test profile funciona sin errores
- [ ] Test con datos externos (PBMC 3k) funciona
- [ ] Todos los outputs esperados se generan
- [ ] Pipeline soporta múltiples samples
- [ ] Batch correction funciona (si aplica)
- [ ] Cell type annotation funciona

### Dashboard
- [ ] Carga h5ad correctamente
- [ ] UMAP se visualiza
- [ ] Cell types se colorean correctamente
- [ ] Gene expression funciona
- [ ] QC plots se muestran (si disponibles)
- [ ] Metadata table muestra datos

### Documentación
- [ ] README.md claro y conciso
- [ ] GETTING_STARTED.md completo
- [ ] docs/dashboard/ organizado
- [ ] Parámetros documentados
- [ ] Troubleshooting básico incluido

### Reproducibilidad
- [ ] Containers especificados y accesibles
- [ ] O conda environments disponibles
- [ ] Versiones de software documentadas
- [ ] Test data incluido en repo

### Usabilidad
- [ ] Usuario puede clonar y correr test en <10 min
- [ ] Mensajes de error son claros
- [ ] No requiere configuración manual compleja

---

## 🎯 CRITERIO DE ÉXITO

**El proyecto está "production-ready" cuando:**

1. ✅ Un colega tuyo puede hacer `git clone` → `nextflow run` → funciona
2. ✅ El test profile termina sin errores
3. ✅ Dashboard carga y muestra los resultados
4. ✅ Documentación explica cómo analizar sus propios datos
5. ✅ Si algo falla, los mensajes de error son útiles

**Prueba final (Gold Standard):**
- Dale el link de GitHub a alguien que NO haya visto el código
- Dile: "Analiza tus datos de single-cell con esto"
- Si lo logran en <1 hora → **SUCCESS** ✅
- Si se frustran → identifica qué los trabó y arréglalo

---

## 📞 ¿Necesitas Ayuda?

Si te trabas en algún paso:
1. Revisa `.nextflow.log` para errores del pipeline
2. Revisa `docs/troubleshooting.md` (si existe)
3. Pregunta en Issues del GitHub (si es público)

---

**Próximo paso inmediato:**  
👉 **Ejecutar Test End-to-End (#1)** - No continúes hasta que esto funcione
