# 🔧 Solución de Problemas - Dashboard en WSL2

## ❓ Problema: "No se puede conectar a localhost:3838"

Este es un problema común en WSL2. Aquí están TODAS las soluciones:

---

## ✅ SOLUCIÓN 1: Usar la IP de WSL (MÁS CONFIABLE)

### Paso 1: Obtén la IP de WSL
```bash
hostname -I
# Ejemplo output: 172.24.160.12
```

### Paso 2: Lanza el dashboard normalmente
```bash
cd /home/damo/scAnnex/dashboard
./launch_dashboard_fixed.sh
```

### Paso 3: Accede con la IP (NO con localhost)
En tu browser Windows:
```
http://172.24.160.12:3838
```
(Usa TU IP, no esta exacta)

---

## ✅ SOLUCIÓN 2: Matar procesos anteriores y reintentar

```bash
# Matar cualquier R/Shiny corriendo
pkill -f shiny
pkill -f "runApp"

# Verificar que el puerto está libre
ss -tln | grep 3838
# Si no muestra nada = puerto libre ✓

# Lanzar de nuevo
cd /home/damo/scAnnex/dashboard
./launch_dashboard_fixed.sh
```

Luego prueba AMBOS:
- `http://localhost:3838`
- `http://TU_IP_WSL:3838`

---

## ✅ SOLUCIÓN 3: Launch manual con output visible

```bash
cd /home/damo/scAnnex/dashboard

# Activar entorno
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard

# Configurar variables
export RETICULATE_PYTHON="$(which python3)"
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run"

# Lanzar (verás TODOS los mensajes)
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

Cuando veas "Listening on http://0.0.0.0:3838":
1. Obtén tu IP: `hostname -I` (en OTRA terminal)
2. Accede con esa IP en browser Windows

---

## ✅ SOLUCIÓN 4: Usar otro puerto

Si 3838 da problemas, usa 8888:

```bash
cd /home/damo/scAnnex/dashboard
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard
export RETICULATE_PYTHON="$(which python3)"
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run"

R -e "shiny::runApp('.', host='0.0.0.0', port=8888)"
```

Accede con:
- `http://localhost:8888`
- `http://TU_IP_WSL:8888`

---

## ✅ SOLUCIÓN 5: Configurar port forwarding de Windows

Si ninguna funciona, configura port forwarding:

### En PowerShell (como Administrador en Windows):
```powershell
netsh interface portproxy add v4tov4 listenport=3838 listenaddress=0.0.0.0 connectport=3838 connectaddress=TU_IP_WSL
```

Ejemplo:
```powershell
netsh interface portproxy add v4tov4 listenport=3838 listenaddress=0.0.0.0 connectport=3838 connectaddress=172.24.160.12
```

Luego accede a `http://localhost:3838` en Windows.

### Para remover después:
```powershell
netsh interface portproxy delete v4tov4 listenport=3838 listenaddress=0.0.0.0
```

---

## 🔍 Verificación: ¿Está corriendo el dashboard?

### Dentro de WSL:
```bash
# Ver si R está corriendo
ps aux | grep -i shiny
ps aux | grep -i "runApp"

# Ver si el puerto está escuchando
ss -tln | grep 3838
# Debería mostrar: LISTEN 0.0.0.0:3838

# Test con curl
curl http://localhost:3838
# Debería devolver HTML
```

### Si curl funciona pero browser no
→ Problema de networking de WSL2, usa la IP directamente

---

## 📋 Checklist de Debugging

```bash
# 1. ¿Conda environment activo?
conda info --envs | grep '*'
# Debe mostrar: * scannex-dashboard

# 2. ¿Python correcto?
echo $RETICULATE_PYTHON
which python3
# Deben ser iguales y apuntar a scannex-dashboard

# 3. ¿Puerto libre?
ss -tln | grep 3838
# No debería mostrar nada ANTES de lanzar

# 4. ¿Datos existen?
ls -lh /home/damo/scAnnex/results_slc_first_run/auto/*annotated*.h5ad
# Debe mostrar el archivo

# 5. ¿IP de WSL?
hostname -I
# Anota esta IP

# 6. ¿Firewall de Windows?
# Verifica en Windows: Settings → Privacy & Security → Windows Security → Firewall
```

---

## 🎯 Método Más Simple (RECOMENDADO)

**UN solo comando que siempre funciona:**

```bash
cd /home/damo/scAnnex/dashboard && \
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)" && \
conda activate scannex-dashboard && \
export RETICULATE_PYTHON="$(which python3)" && \
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run" && \
echo "" && \
echo "═══════════════════════════════════════════════════" && \
echo "  Dashboard Starting..." && \
echo "═══════════════════════════════════════════════════" && \
echo "" && \
echo "  Your WSL IP: $(hostname -I | awk '{print $1}')" && \
echo "" && \
echo "  Access at:" && \
echo "    http://localhost:3838" && \
echo "    http://$(hostname -I | awk '{print $1}'):3838" && \
echo "" && \
echo "  Press Ctrl+C to stop" && \
echo "═══════════════════════════════════════════════════" && \
echo "" && \
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

Copia y pega TODO ese bloque en tu terminal.

---

## 🚨 Errores Comunes y Soluciones

### Error: "Port already in use"
```bash
pkill -f shiny
sleep 2
# Reintentar
```

### Error: "ModuleNotFoundError: anndata"
```bash
# Verificar environment
conda activate scannex-dashboard
python3 -c "import anndata; print('OK')"

# Si falla, reinstalar
conda env remove -n scannex-dashboard
conda env create -f environment_dashboard.yml
```

### Error: "Cannot open connection"
```bash
# Usar IP en lugar de localhost
WSL_IP=$(hostname -I | awk '{print $1}')
echo "Try: http://${WSL_IP}:3838"
```

### Dashboard inicia pero se cierra inmediatamente
```bash
# Ver logs completos
cd /home/damo/scAnnex/dashboard
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)"
conda activate scannex-dashboard
export RETICULATE_PYTHON="$(which python3)"

# Lanzar con debug
R --vanilla << 'EOF'
options(shiny.trace = TRUE)
options(error = traceback)
shiny::runApp('.', host='0.0.0.0', port=3838)
EOF
```

---

## ✅ Prueba Rápida

Ejecuta este test para verificar todo:

```bash
cd /home/damo/scAnnex/dashboard
./test_dashboard.sh
```

Si dice "All tests passed", el problema es solo de networking WSL2 → usa la IP.

---

## 💡 Pro Tip: Alias Permanente

Agrega a tu `~/.bashrc`:

```bash
alias dashboard='cd /home/damo/scAnnex/dashboard && eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)" && conda activate scannex-dashboard && export RETICULATE_PYTHON=$(which python3) && export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run" && echo "WSL IP: $(hostname -I | awk '\''{print $1}'\'')" && echo "Access: http://localhost:3838 or http://$(hostname -I | awk '\''{print $1}'\'''):3838" && R -e "shiny::runApp(\".\", host=\"0.0.0.0\", port=3838)"'
```

Luego solo escribe:
```bash
dashboard
```

---

## 📞 Si Nada Funciona

1. Verifica versión de WSL: `wsl --version` (en PowerShell Windows)
2. Considera actualizar a WSL2 si estás en WSL1
3. Prueba acceder desde la misma WSL con navegador lynx/w3m
4. Como último recurso: export dashboard para acceder via Jupyter

---

¿Cuál método te funcionó? Prueba primero el "Método Más Simple" copiando todo el comando.
