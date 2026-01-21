# Dashboard Integration - User Guide

## 🎯 Overview

El pipeline scAnnex ahora incluye integración automática con el dashboard interactivo. Al finalizar la ejecución del pipeline, recibirás instrucciones en consola sobre cómo lanzar el dashboard para explorar tus resultados.

---

## ✨ Características

### Automático por Defecto
- El dashboard está **habilitado por defecto** (`--enable_dashboard true`)
- Al finalizar el pipeline, verás un mensaje en consola con:
  - 📂 Ubicación de resultados
  - 🚀 Comando para lanzar el dashboard
  - 🌐 URL donde estará disponible
  - 💡 Tips y documentación

### Personalizable
Puedes controlar el comportamiento del dashboard con parámetros:

```bash
nextflow run main.nf \
  --input datos.h5ad \
  --enable_dashboard true \       # Activar/desactivar (default: true)
  --dashboard_port 3838 \         # Puerto del servidor (default: 3838)
  --dashboard_host localhost      # Host (default: localhost)
```

---

## 📋 Flujo de Trabajo

### 1. Ejecutar Pipeline

```bash
nextflow run main.nf \
  --input test_data/outputs/PBMC_MTX_quick_test.h5ad \
  --outdir results/my_analysis
```

### 2. Pipeline Completa

Al finalizar exitosamente, verás en consola:

```
════════════════════════════════════════════════════════════════
  🎉 Pipeline Completed Successfully!
════════════════════════════════════════════════════════════════

📊 Your interactive dashboard is ready to launch!

📂 Results location:
   results/my_analysis

🚀 To launch the dashboard, run:

   cd /home/damo/scAnnex/dashboard
   bash launch_dashboard.sh results/my_analysis

════════════════════════════════════════════════════════════════
  🌐 Dashboard URL (after launch):
  http://localhost:3838
════════════════════════════════════════════════════════════════

💡 Tips:
   - Click or copy the URL to open in your browser
   - Use Ctrl+C to stop the dashboard
   - Dashboard info saved to: results/my_analysis/dashboard_info.txt
```

### 3. Lanzar Dashboard

**Opción A: Script Conveniente (Recomendado)**
```bash
cd dashboard
bash launch_dashboard.sh ../results/my_analysis
```

El script automáticamente:
- Detecta el mejor método disponible (Apptainer, Singularity, Docker, Conda+R)
- Configura el entorno
- Lanza el dashboard
- Muestra la URL

**Opción B: R Directo (Si tienes R + Shiny instalado)**
```bash
cd dashboard
Rscript -e "shiny::runApp(port=3838, host='localhost')"
```

**Opción C: Docker**
```bash
cd dashboard
docker build -t scannex-dashboard .
docker run -p 3838:3838 -v $(pwd)/../results:/data scannex-dashboard
```

### 4. Acceder al Dashboard

1. Abre tu navegador
2. Navega a: `http://localhost:3838`
3. Explora tus resultados interactivamente

---

## 🎨 Características del Dashboard

Una vez lanzado, el dashboard te permite:

### Tab 1: Data Overview
- Resumen del dataset
- Estadísticas de células y genes
- Información de QC

### Tab 2: UMAP Visualization
- Visualización interactiva de UMAP
- Colorear por:
  - Clustering (múltiples resoluciones)
  - QC metrics (genes, counts, mito%)
  - Anotaciones CellTypist
  - Expresión de genes específicos
- Zoom, pan, seleccionar células

### Tab 3: Quality Control
- Distribuciones de QC metrics
- Filtros interactivos
- Células antes/después de filtering

### Tab 4: Cluster Analysis
- Estadísticas por cluster
- Genes marcadores por cluster
- Comparación entre clusters
- Heatmaps de expresión

### Tab 5: Gene Expression
- Búsqueda de genes
- Violin plots por cluster
- Feature plots (UMAP + expresión)
- Exportar listas de genes

---

## ⚙️ Configuración Avanzada

### Cambiar Puerto

Si el puerto 3838 está ocupado:

```bash
nextflow run main.nf \
  --input data.h5ad \
  --dashboard_port 8080
```

Dashboard estará en: `http://localhost:8080`

### Acceso Remoto (HPC/Servidor)

Para acceder al dashboard desde tu computadora local cuando el pipeline corre en un servidor remoto:

**1. Ejecutar pipeline con host 0.0.0.0:**
```bash
nextflow run main.nf \
  --input data.h5ad \
  --dashboard_host 0.0.0.0 \
  --dashboard_port 3838
```

**2. Crear SSH tunnel desde tu computadora:**
```bash
ssh -L 3838:localhost:3838 usuario@servidor.com
```

**3. Acceder en tu navegador local:**
```
http://localhost:3838
```

### Desactivar Dashboard

Si no quieres el mensaje del dashboard:

```bash
nextflow run main.nf \
  --input data.h5ad \
  --enable_dashboard false
```

---

## 📁 Archivos Generados

El proceso de dashboard crea:

```
results/my_analysis/
└── dashboard_info.txt    # Información y comandos para lanzar dashboard
```

**Contenido de dashboard_info.txt:**
```
Dashboard Configuration
======================
Results Directory: results/my_analysis
Dashboard Port: 3838
Dashboard Host: localhost
Dashboard URL: http://localhost:3838

Launch Command:
  cd /home/damo/scAnnex/dashboard
  bash launch_dashboard.sh results/my_analysis

Documentation:
  /home/damo/scAnnex/dashboard/README.md
  /home/damo/scAnnex/dashboard/QUICKSTART.md
```

---

## 🐛 Troubleshooting

### Dashboard no Lanza

**Síntoma:** Script de launch falla

**Soluciones:**
1. Verifica que tienes al menos uno de: Apptainer, Singularity, Docker, o R+Shiny
2. Lee el mensaje de error del script
3. Consulta `dashboard/TROUBLESHOOTING_WSL2.md` si estás en WSL2

### Puerto Ocupado

**Síntoma:** Error "port already in use"

**Solución:**
```bash
# Usa otro puerto
bash launch_dashboard.sh results/ --port 8080
```

### No Puedo Acceder desde Browser

**Síntoma:** "Connection refused" en navegador

**Soluciones:**
1. Verifica que el dashboard está corriendo (debe ver mensaje "Listening on http://...")
2. Si estás en servidor remoto, necesitas SSH tunnel
3. Verifica firewall no bloquea el puerto

### Dashboard Muestra Error al Cargar Datos

**Síntoma:** Error en dashboard "Cannot find H5AD file"

**Solución:**
```bash
# Asegúrate de pasar el directorio correcto de resultados
bash launch_dashboard.sh /ruta/completa/a/results
```

---

## 📚 Documentación Adicional

- **Dashboard README:** `dashboard/README.md`
- **Quick Start:** `dashboard/QUICKSTART.md`
- **WSL2 Issues:** `dashboard/TROUBLESHOOTING_WSL2.md`
- **Manual Launch:** `dashboard/MANUAL_LAUNCH.md`

---

## 💡 Ejemplos de Uso

### Ejemplo 1: Pipeline Completo con Dashboard

```bash
# Ejecutar pipeline
nextflow run main.nf \
  --input pbmc_data.h5ad \
  --outdir results/pbmc_analysis \
  --run_auto_annotation true

# Después de que termine, lanzar dashboard
cd dashboard
bash launch_dashboard.sh ../results/pbmc_analysis

# Abrir en navegador: http://localhost:3838
```

### Ejemplo 2: Múltiples Análisis con Diferentes Puertos

```bash
# Análisis 1 en puerto 3838
cd dashboard
bash launch_dashboard.sh ../results/analysis1 --port 3838 &

# Análisis 2 en puerto 3839
bash launch_dashboard.sh ../results/analysis2 --port 3839 &

# Ahora puedes comparar en:
# http://localhost:3838 (análisis 1)
# http://localhost:3839 (análisis 2)
```

### Ejemplo 3: HPC con SSH Tunnel

```bash
# En servidor HPC:
nextflow run main.nf \
  --input data.h5ad \
  --outdir results \
  --dashboard_host 0.0.0.0

cd dashboard
bash launch_dashboard.sh ../results

# En tu computadora local:
ssh -L 3838:localhost:3838 usuario@hpc.universidad.edu

# Abrir navegador: http://localhost:3838
```

---

## 🎓 Best Practices

1. **Siempre usar el script de launch:** `launch_dashboard.sh` detecta automáticamente el mejor método
2. **No dejar dashboards corriendo:** Usa Ctrl+C para detener cuando termines
3. **SSH tunnels para HPC:** Más seguro que exponer puertos públicamente
4. **Revisar dashboard_info.txt:** Contiene toda la info de configuración
5. **Usar diferentes puertos:** Para comparar múltiples análisis simultáneamente

---

## ❓ FAQ

**P: ¿El dashboard se lanza automáticamente?**  
R: No, el pipeline solo muestra las instrucciones. Debes ejecutar manualmente el comando de launch.

**P: ¿Puedo desactivar el mensaje del dashboard?**  
R: Sí, usa `--enable_dashboard false`

**P: ¿El dashboard modifica mis resultados?**  
R: No, es solo visualización. Los archivos H5AD son leídos en modo read-only.

**P: ¿Puedo compartir el dashboard con colaboradores?**  
R: Sí, pero necesitas configuración de red apropiada. Ver sección "Acceso Remoto".

**P: ¿Funciona en Windows/Mac/Linux?**  
R: Sí, el dashboard funciona en todas las plataformas donde corra el pipeline.

---

**Última Actualización:** 2026-01-21  
**Versión:** v0.1.0  
**Maintainer:** damouzo
