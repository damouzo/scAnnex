# Manual Dashboard Launch Guide

Esta guía te muestra cómo lanzar el dashboard manualmente paso a paso, útil para debugging o ejecución personalizada.

---

## 🎯 Método 1: Usando el Script Automático (RECOMENDADO)

```bash
cd /home/damo/scAnnex/dashboard

# Activa el entorno conda
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard

# Lanza el dashboard
./launch_dashboard.sh
```

✅ **Ventajas**: Auto-detecta todo, maneja puertos, muestra instrucciones
❌ **Desventajas**: Menos control sobre configuración

---

## 🔧 Método 2: Lanzamiento Manual Completo

### Paso 1: Preparar el Entorno

```bash
# Navegar al directorio del dashboard
cd /home/damo/scAnnex/dashboard

# Activar conda (miniforge en tu caso)
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard

# Verificar que estás en el entorno correcto
which python3
# Debe mostrar: /home/damo/miniforge3/envs/scannex-dashboard/bin/python3

which R
# Debe mostrar: /home/damo/miniforge3/envs/scannex-dashboard/bin/R
```

### Paso 2: Configurar Variables de Entorno (IMPORTANTE)

```bash
# Ruta a tus resultados
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run"

# Forzar a R que use el Python correcto (CRÍTICO)
export RETICULATE_PYTHON="/home/damo/miniforge3/envs/scannex-dashboard/bin/python3"

# Verificar
echo $SCANNEX_DATA_PATH
echo $RETICULATE_PYTHON
```

### Paso 3: Lanzar Shiny Server

**Opción A: Desde R interactivo** (mejor para debugging)
```bash
R
```

Dentro de R:
```r
# Verificar que el Python correcto está configurado
cat("RETICULATE_PYTHON:", Sys.getenv("RETICULATE_PYTHON"), "\n")

# Cargar y lanzar
shiny::runApp('.', host='0.0.0.0', port=3838)
```

**Opción B: Una sola línea** (más rápido)
```bash
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

**Opción C: Con log detallado** (para debugging)
```bash
R --vanilla << 'EOF'
# Mostrar configuración
cat("Working directory:", getwd(), "\n")
cat("RETICULATE_PYTHON:", Sys.getenv("RETICULATE_PYTHON"), "\n")
cat("SCANNEX_DATA_PATH:", Sys.getenv("SCANNEX_DATA_PATH"), "\n\n")

# Lanzar con opciones de debugging
options(shiny.trace = TRUE)
shiny::runApp('.', host='0.0.0.0', port=3838)
EOF
```

### Paso 4: Acceder desde tu Browser

En tu browser de Windows (Chrome, Edge, Firefox):
```
http://localhost:3838
```

Si el puerto 3838 está ocupado, usa otro:
```bash
R -e "shiny::runApp('.', host='0.0.0.0', port=8888)"
# Luego accede a: http://localhost:8888
```

### Paso 5: Usar el Dashboard

1. **En la tab "Data Input"**:
   - H5AD file path: `/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad`
   - QC directory: `/home/damo/scAnnex/results_slc_first_run/qc`
   - Click "Load Data"

2. **En la tab "Clustering & UMAP"**:
   - Selecciona color by: `predicted_labels`
   - Ajusta point size y opacity si quieres
   - Explora el UMAP interactivo

3. **En la tab "Gene Expression"**:
   - Busca genes: `CD3D`, `CD14`, `CD79A`, `MS4A1`
   - Click "Plot Expression"

4. **En la tab "Metadata"**:
   - Explora la tabla completa de células
   - Usa los filtros para buscar células específicas

### Paso 6: Detener el Dashboard

- Presiona `Ctrl + C` en la terminal donde corre R
- El servidor se detendrá limpiamente

---

## 🚨 Troubleshooting Manual

### Problema: "Error: ModuleNotFoundError: No module named 'anndata'"

**Causa**: R está usando el Python equivocado

**Solución**:
```bash
# ANTES de lanzar R, forzar el Python correcto
export RETICULATE_PYTHON="$(which python3)"
echo $RETICULATE_PYTHON

# Verificar que anndata está disponible
python3 -c "import anndata; print('OK')"

# Ahora lanzar R
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

### Problema: "Error in file(filename, 'r', encoding = encoding)"

**Causa**: R no encuentra los archivos del dashboard

**Solución**:
```bash
# Asegúrate de estar EN el directorio dashboard
cd /home/damo/scAnnex/dashboard
ls -la  # Debes ver: app.R, global.R, server.R, ui.R

# Ahora lanza
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

### Problema: "Port 3838 already in use"

**Solución 1**: Usar otro puerto
```bash
R -e "shiny::runApp('.', host='0.0.0.0', port=8888)"
```

**Solución 2**: Matar proceso en puerto 3838
```bash
# Encontrar qué está usando el puerto
ss -tlnp | grep 3838
# o
lsof -i :3838

# Matar el proceso (reemplaza PID)
kill <PID>
```

### Problema: Dashboard carga pero no muestra datos

**Causa**: Ruta incorrecta al archivo h5ad

**Solución**:
```bash
# Verificar que el archivo existe
ls -lh /home/damo/scAnnex/results_slc_first_run/auto/*annotated*.h5ad

# Verificar que tiene predicted_labels
python3 << 'EOF'
import anndata as ad
adata = ad.read_h5ad('/home/damo/scAnnex/results_slc_first_run/auto/PBMC_TEST_annotated.h5ad')
print("Columns:", adata.obs.columns.tolist())
print("Has predicted_labels:", 'predicted_labels' in adata.obs.columns)
EOF
```

### Problema: "Can't access from Windows browser"

**Causa**: WSL2 networking

**Solución**:
```bash
# En WSL, verifica la IP
hostname -I

# Usa esa IP en tu browser Windows
# Ejemplo: http://172.x.x.x:3838
```

O simplemente usa `localhost` (debería funcionar en WSL2):
```
http://localhost:3838
```

---

## 📋 Checklist Pre-Launch

Antes de lanzar manualmente, verifica:

```bash
# ✓ Conda environment activado
conda info --envs | grep '*'

# ✓ Python correcto
which python3
python3 -c "import anndata; print('anndata OK')"

# ✓ R disponible
which R
R --version

# ✓ Archivos del dashboard presentes
ls -1 app.R global.R server.R ui.R

# ✓ Datos existen
ls -lh ../results_slc_first_run/auto/*annotated*.h5ad

# ✓ Variable RETICULATE_PYTHON configurada
echo $RETICULATE_PYTHON
```

Si TODOS muestran OK → Lanza con confianza ✅

---

## 🔬 Modo Debug Avanzado

Para investigar problemas detalladamente:

```bash
cd /home/damo/scAnnex/dashboard

# Activar entorno
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard

# Variables de entorno
export RETICULATE_PYTHON="$(which python3)"
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run"

# Lanzar R con debug output
R --vanilla << 'EOF'
# Debug info
cat("═══════════════════════════════════════\n")
cat("DEBUG INFO\n")
cat("═══════════════════════════════════════\n")
cat("Working dir:", getwd(), "\n")
cat("RETICULATE_PYTHON:", Sys.getenv("RETICULATE_PYTHON"), "\n")
cat("SCANNEX_DATA_PATH:", Sys.getenv("SCANNEX_DATA_PATH"), "\n")
cat("R version:", R.version.string, "\n")
cat("═══════════════════════════════════════\n\n")

# Cargar global.R con traceback completo
options(error = traceback)
source("global.R")

cat("\n✓ global.R loaded\n\n")

# Intentar importar Python modules manualmente
library(reticulate)
cat("Python config:\n")
print(py_config())
cat("\n")

tryCatch({
    ad <- import("anndata")
    cat("✓ anndata imported\n")
    
    sc <- import("scanpy")
    cat("✓ scanpy imported\n")
    
    cat("\n✓ All Python modules OK\n")
}, error = function(e) {
    cat("✗ Python import error:\n")
    print(e)
})

cat("\n Ready to launch. Starting Shiny...\n\n")

# Lanzar con trace
options(shiny.trace = TRUE)
shiny::runApp('.', host='0.0.0.0', port=3838)
EOF
```

---

## 🎓 Explicación: ¿Por Qué Estos Pasos?

### 1. ¿Por qué activar conda?
- Necesitas el entorno `scannex-dashboard` con R y Python packages
- Sin activar conda, usarías el R/Python del sistema (sin packages)

### 2. ¿Por qué export RETICULATE_PYTHON?
- R package `reticulate` conecta R con Python
- Por defecto, busca Python en `/usr/bin/python3`
- Necesitamos que use el Python de conda (que tiene scanpy/anndata)
- **CRÍTICO**: Debe configurarse ANTES de cargar reticulate

### 3. ¿Por qué host='0.0.0.0'?
- Permite acceso desde fuera de WSL (tu browser Windows)
- Con `127.0.0.1` solo funcionaría dentro de WSL

### 4. ¿Por qué port=3838?
- Puerto estándar para Shiny
- Puedes usar cualquier puerto libre (8000-9999)

---

## 📚 Comandos Útiles Durante Sesión

### Ver logs en tiempo real
Si lanzas en background:
```bash
R -e "shiny::runApp('.', port=3838)" > dashboard.log 2>&1 &
tail -f dashboard.log
```

### Verificar que está corriendo
```bash
curl http://localhost:3838
# Debería devolver HTML de Shiny
```

### Ver memoria/CPU usado
```bash
ps aux | grep shiny
htop  # presiona F4, busca "shiny"
```

### Matar elegantemente
```bash
# Si está en foreground: Ctrl+C
# Si está en background:
pkill -f "shiny::runApp"
```

---

## 💡 Pro Tips

### Tip 1: Alias para launch rápido
Agrega a tu `~/.bashrc`:
```bash
alias dashboard='cd /home/damo/scAnnex/dashboard && eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)" && conda activate scannex-dashboard && export RETICULATE_PYTHON=$(which python3) && R -e "shiny::runApp(\".\", host=\"0.0.0.0\", port=3838)"'
```

Luego solo:
```bash
dashboard
```

### Tip 2: Launch script personalizado
Crea `~/launch_scannex.sh`:
```bash
#!/bin/bash
cd /home/damo/scAnnex/dashboard
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard
export RETICULATE_PYTHON=$(which python3)
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run"
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

```bash
chmod +x ~/launch_scannex.sh
~/launch_scannex.sh
```

### Tip 3: Auto-abrir browser
```bash
R -e "shiny::runApp('.', host='0.0.0.0', port=3838, launch.browser=TRUE)"
```

---

## ✅ Quick Command Summary

**Setup una vez:**
```bash
cd /home/damo/scAnnex/dashboard
./setup_dashboard.sh
```

**Lanzar siempre (opción fácil):**
```bash
cd /home/damo/scAnnex/dashboard
./launch_dashboard.sh
```

**Lanzar manual (control total):**
```bash
cd /home/damo/scAnnex/dashboard
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard
export RETICULATE_PYTHON="$(which python3)"
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

**Acceder:**
```
http://localhost:3838
```

**Detener:**
```
Ctrl + C
```

---

¿Listo? ¡Ahora sabes todos los métodos posibles! 🚀
