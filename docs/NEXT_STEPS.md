# Próximos Pasos - scAnnex Development

**Fecha:** 2026-01-21  
**Estado Actual:** Test 3.1 COMPLETADO ✅  
**Pipeline Status:** FUNCIONAL con conda profile

---

## 🎯 Estado Actual del Proyecto

### ✅ Completado (Week 1 + Inicio Week 2)

**Week 1: Critical Fixes** ✅ DONE
- Fix #1: base.config syntax error
- Fix #2: nf-validation plugin
- Fix #3: Input flexibility
- Fix #4: Container updates (scanpy 1.10.0)
- Fix #5: Conda profile creation
- Fix #6: CI/CD setup
- Fix #7: Documentation organization

**Week 2: Testing & Validation** 🚧 25% DONE
- ✅ Phase 1: Environment prep
- ✅ Phase 2: Compilation tests
- ✅ Phase 3.1: Single file input test (H5AD)
- ⏳ Phase 3.2-3.4: Otros formatos de input
- ⏳ Phase 4-8: Tests completos

---

## 📋 Próximos Pasos Inmediatos

### 1. Test 3.2: Samplesheet Input (PRIORIDAD ALTA)

**Objetivo:** Validar que el pipeline puede procesar múltiples samples desde CSV

**Tiempo Estimado:** 30-45 min

**Pasos:**
```bash
# 1. Crear samplesheet con múltiples samples
cat > test_data/samplesheet_multi.csv << 'EOF'
sample,file,batch,condition
sample1,test_data/outputs/PBMC_MTX_quick_test.h5ad,batch1,control
sample2,test_data/outputs/PBMC_MTX_quick_test.h5ad,batch2,treated
EOF

# 2. Ejecutar pipeline
nextflow run main.nf \
  --input test_data/samplesheet_multi.csv \
  --outdir test_results/test_2_samplesheet \
  --run_integration true \
  --batch_key batch \
  --max_memory 8.GB \
  --max_cpus 4 \
  -profile conda \
  -resume
```

**Validaciones:**
- ✅ CSV parseado correctamente
- ✅ Múltiples samples procesados en paralelo
- ✅ Batch correction ejecutado (Harmony)
- ✅ Outputs separados por sample + integrado

---

### 2. Test 3.3: MTX Input Format (PRIORIDAD MEDIA)

**Objetivo:** Validar input desde 10X MTX format

**Tiempo Estimado:** 20-30 min

**Pasos:**
```bash
nextflow run main.nf \
  --input test_data/mtx/filtered_feature_bc_matrix \
  --input_type mtx \
  --outdir test_results/test_3_mtx \
  --max_memory 8.GB \
  --max_cpus 4 \
  -profile conda \
  -resume
```

**Validaciones:**
- ✅ MTX directory detectado
- ✅ matrix.mtx, barcodes.tsv, features.tsv leídos
- ✅ Conversión a H5AD correcta
- ✅ Pipeline continúa normalmente

---

### 3. Documentar Fixes en Conda Specs (PRIORIDAD ALTA)

**Objetivo:** Crear registro permanente de cambios en versioning

**Archivo a Crear:** `docs/CONDA_VERSIONING_FIXES.md`

**Contenido:**
- Problema original (versiones exactas no disponibles)
- Solución implementada (versioning flexible)
- Lista de todos los módulos actualizados
- Recomendaciones para futuro mantenimiento

---

### 4. Actualizar env/scanpy.yml (PRIORIDAD MEDIA)

**Objetivo:** Sincronizar environment file con module specs

**Archivo:** `env/scanpy.yml`

**Cambios Necesarios:**
```yaml
dependencies:
  - scanpy>=1.9        # Era: scanpy=1.10.0
  - celltypist>=1.6    # Era: celltypist=1.6.2
  - scrublet>=0.2      # Era: scrublet=0.2.3
  - harmonypy>=0.0.9   # Era: harmonypy=0.0.9 (OK)
  - python-igraph      # AÑADIR (faltaba)
  - leidenalg          # AÑADIR (faltaba)
```

---

### 5. Commit de Fixes (PRIORIDAD ALTA)

**Objetivo:** Guardar todos los fixes aplicados hoy

**Archivos Modificados:**
- `modules/local/*.nf` (6 files - conda specs)
- `docs/WSL2_SINGULARITY_NOTES.md` (nuevo)
- `test_results/TEST_3.1_REPORT.md` (nuevo)

**Mensaje de Commit Sugerido:**
```
fix: Update conda environment specifications for compatibility

- Changed from exact versions to flexible versioning (>=)
- Removed channel prefixes to avoid mamba conflicts
- Added missing dependencies: python-igraph, leidenalg
- Fixed all modules: unify_input, quality_control, doublet_detection,
  standard_processing, auto_annot_celltypist, normalize_integrate

Fixes compatibility issues with conda 25.11.0 channel availability.

Test 3.1 (single file H5AD input) now passes completely.
```

---

## 🔄 Week 2 Testing Plan (Resto de Fases)

### Phase 4: Edge Cases & Error Handling (3-4h)
- Test con archivo corrupto
- Test con sample sin metadata
- Test con parámetros inválidos
- Validar mensajes de error

### Phase 5: Resource Limits (2-3h)
- Test con max_memory reducido
- Test con max_cpus = 1
- Validar que resource checks funcionan

### Phase 6: Optional Features (3-4h)
- Test con doublet_removal = false
- Test con run_auto_annotation = false
- Test con diferentes métodos de integración (scanorama, bbknn)
- Test con diferentes clustering methods (louvain vs leiden)

### Phase 7: Output Validation (2-3h)
- Validar estructura de H5AD outputs
- Verificar que CSVs tienen formato correcto
- Comprobar que plots se generan (si aplicable)
- Validar pipeline_info reports

### Phase 8: Benchmarking (2-3h)
- Medir tiempos de ejecución por módulo
- Medir uso de memoria peak
- Medir tamaño de outputs
- Comparar conda vs containers (si disponible)

**Tiempo Total Restante Week 2:** ~15-20 horas

---

## 🐛 Known Issues (Documentar)

### 1. Singularity en WSL2
**Issue:** Securebits error en WSL2  
**Status:** Documentado en `docs/WSL2_SINGULARITY_NOTES.md`  
**Workaround:** Usar conda profile en WSL2  
**Solution:** Testear con Singularity en Linux nativo/HPC

### 2. Python 3.14 Warning
**Issue:** Conda creó env con Python 3.14 (muy reciente)  
**Impact:** Ninguno aparente, pero podría causar issues futuros  
**Recommendation:** Especificar `python=3.10` en conda specs

### 3. Report Rendering Warnings
**Issue:** "Failed to render execution report/timeline"  
**Causa:** Posible issue con Nextflow 25.04.6 o missing dependencies  
**Impact:** Mínimo - reports se generan de todas formas  
**Solution:** Investigar en próximo test

---

## 📚 Documentación Pendiente (Week 3)

### High Priority
1. **README.md completo**
   - Installation instructions
   - Quick start guide
   - Example commands
   - Expected outputs

2. **Parameter documentation**
   - Descripción de cada parámetro
   - Valores por defecto
   - Ejemplos de uso

3. **Troubleshooting guide**
   - Errores comunes
   - Soluciones
   - FAQ

### Medium Priority
4. **Output specification**
   - Descripción de cada output file
   - Formato de CSVs
   - Estructura de H5AD

5. **Best practices guide**
   - Recomendaciones de parámetros por tipo de dataset
   - Interpretación de resultados
   - Cómo elegir resolución de clustering

---

## 🚀 Week 4: Release Preparation

### Pre-Release Checklist
- [ ] Todos los tests de Week 2 pasando
- [ ] Documentación completa
- [ ] CI/CD funcionando
- [ ] Container images publicados (Docker Hub / Quay.io)
- [ ] Changelog creado
- [ ] README badges actualizados
- [ ] Versión tagged (v0.1.0)

### Release Artifacts
1. GitHub Release con changelog
2. Container images en Quay.io
3. Conda environment file validado
4. Test datasets disponibles
5. Ejemplo de samplesheet

---

## 💡 Recomendaciones Estratégicas

### 1. Priorizar Container Testing
**Razón:** Conda funciona, pero containers son la solución preferida  
**Acción:** En próxima sesión en Linux nativo, testear con Singularity  
**Beneficio:** Validar que containers funcionan en producción (HPC)

### 2. Añadir Test Data al Repo
**Razón:** Facilita CI/CD y tests de usuarios  
**Acción:** Incluir `test_data/` en repo (con .gitattributes LFS si >50MB)  
**Beneficio:** Tests reproducibles sin descargas externas

### 3. Crear Profile "test_quick"
**Razón:** Tests actuales toman 3-15 min  
**Acción:** Crear profile con parámetros ultra-reducidos (10 células, skip clustering)  
**Beneficio:** CI/CD más rápido, desarrollo más ágil

### 4. Documentar Decisiones de Diseño
**Razón:** Futuro mantenimiento y contribuciones  
**Archivo:** `docs/DESIGN_DECISIONS.md`  
**Contenido:**
  - Por qué multi-resolution clustering
  - Por qué CellTypist para annotation
  - Por qué Harmony para batch correction
  - Por qué estructura de outputs actual

---

## 🎯 Objetivo de Esta Semana

**Meta:** Completar Week 2 Testing Phase

**Entregables:**
1. ✅ Test 3.1 report (DONE)
2. ⏳ Test 3.2 report (samplesheet)
3. ⏳ Test 3.3 report (MTX)
4. ⏳ Complete testing summary
5. ⏳ Bug list (si hay)
6. ⏳ Performance benchmarks

**Criterio de Éxito:** Poder decir con confianza:
- "El pipeline funciona con todos los formatos de input"
- "Los outputs son científicamente válidos"
- "El performance es aceptable"
- "Los errores son manejados correctamente"

---

## 📞 Próxima Sesión - Checklist

Antes de empezar la próxima sesión:

1. **Revisar este documento** ✅
2. **Decidir qué test ejecutar** (recomiendo 3.2: samplesheet)
3. **Verificar que conda cache existe** (para re-usar environments)
4. **Tener git status limpio** (opcional, para commits incrementales)

Comando para empezar:
```bash
cd /home/damo/scAnnex
git status                    # Ver estado
cat docs/NEXT_STEPS.md        # Revisar este archivo
# Luego ejecutar test elegido
```

---

**Última Actualización:** 2026-01-21 10:45 UTC  
**Preparado por:** OpenCode Assistant  
**Session Status:** Test 3.1 COMPLETADO ✅
