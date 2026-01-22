# scAnnex Version Compatibility Fix

## Problema Resuelto

Se identificó y corrigió una incompatibilidad entre las versiones de `anndata` utilizadas por el pipeline de procesamiento y el dashboard interactivo, que impedía que el dashboard cargara correctamente los datos y mostrara los UMAPs.

## Solución Aplicada

**Sincronización de versiones entre environments:**

- **Pipeline (scannex)**: Actualizado a Python 3.11, anndata 0.12.7, scanpy 1.11.5
- **Dashboard (scannex-dashboard)**: Ya tenía Python 3.11, anndata 0.12.7, scanpy 1.11.5

Ambos environments ahora utilizan las mismas versiones, garantizando compatibilidad completa.

## Cómo Usar el Dashboard

### 1. Procesar Datos con el Pipeline

```bash
nextflow run main.nf -profile conda \
    --input samplesheet.csv \
    --outdir results
```

### 2. Iniciar el Dashboard

```bash
cd dashboard
SCANNEX_DATA_PATH=../results bash launch_dashboard.sh
```

El dashboard estará disponible en: **http://localhost:3838**

### 3. Cargar Datos en el Dashboard

1. Abre el navegador en http://localhost:3838
2. Ve al tab **"Data Input"**
3. Ingresa las rutas:
   - **H5AD File Path**: `/home/damo/scAnnex/results/standard/PBMC_1k_processed.h5ad`
   - **QC Results Directory**: `/home/damo/scAnnex/results/qc`
4. Haz clic en **"Load Data"**

### 4. Funcionalidades Disponibles

Una vez cargados los datos, todas las funcionalidades están disponibles:

- **QC Overview**: Métricas de calidad antes y después del filtrado
- **Clustering & UMAP**: Visualización interactiva de UMAP coloreada por metadata
- **Gene Expression**: Búsqueda de genes y visualización de expresión en UMAP
- **Annotation Station**: Anotación manual de clusters

## Archivos Generados

El pipeline genera los siguientes archivos:

```
results/
├── unified_input/
│   └── PBMC_1k_unified.h5ad
├── qc/
│   ├── PBMC_1k_qc.h5ad
│   └── qc_results/
│       ├── qc_report.json
│       └── *.png (plots)
├── doublet_detection/
│   └── PBMC_1k_doublets.h5ad
└── standard/
    ├── PBMC_1k_processed.h5ad  ← Usar este en el dashboard
    └── standard_processing_results/
        ├── cell_metadata.csv
        ├── umap_coordinates.csv
        └── *.png (plots)
```

## Verificación de Compatibilidad

Para verificar que las versiones están sincronizadas:

```bash
# Pipeline
conda run -n scannex python3 -c "import anndata; print(anndata.__version__)"

# Dashboard  
conda run -n scannex-dashboard python3 -c "import anndata; print(anndata.__version__)"
```

Ambos deben mostrar: **0.12.7**

## Notas Importantes

- **No es necesario** convertir archivos H5AD manualmente
- Los archivos generados por el pipeline son **directamente compatibles** con el dashboard
- Si encuentras errores, verifica que las versiones de ambos environments estén sincronizadas
- El dashboard auto-detecta archivos en el directorio especificado por `SCANNEX_DATA_PATH`

## Mantenimiento Futuro

Para evitar problemas de compatibilidad en el futuro:

1. Mantener las mismas versiones de `anndata` y `scanpy` en ambos environments
2. Al actualizar packages, actualizar ambos environments simultáneamente
3. Probar la compatibilidad después de cualquier actualización

## Solución de Problemas

### Dashboard no carga datos

1. Verifica que el archivo H5AD existe
2. Verifica que las rutas son correctas (absolutas o relativas al directorio del dashboard)
3. Revisa los logs del dashboard para errores específicos

### UMAPs no se muestran

1. Verifica que el archivo H5AD contiene `X_umap` en `.obsm`
2. Asegúrate de hacer clic en "Load Data" después de ingresar las rutas
3. Verifica que la columna de metadata seleccionada existe (ej: 'batch', 'leiden_0.5')

### Búsqueda de genes no funciona

1. Verifica que los nombres de genes están en el archivo H5AD
2. Los nombres de genes son case-sensitive (ej: 'CD3D' no 'cd3d')
3. El dashboard debe haber cargado completamente los datos antes de buscar genes
