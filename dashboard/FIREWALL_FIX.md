### SOLUCIÓN RÁPIDA - Firewall de Windows para WSL2

El dashboard ESTÁ corriendo, pero Windows está bloqueando la conexión.

## ✅ OPCIÓN 1: Permitir Puerto en Firewall (RECOMENDADO)

### En PowerShell de Windows (como Administrador):

```powershell
# Permitir tráfico en puerto 3838
New-NetFirewallRule -DisplayName "WSL2 Shiny Dashboard" -Direction Inbound -LocalPort 3838 -Protocol TCP -Action Allow
```

Luego vuelve a intentar acceder en el browser:
```
http://169.254.76.190:3838
```

---

## ✅ OPCIÓN 2: Usar localhost con Port Forwarding

### En PowerShell de Windows (como Administrador):

```powershell
# Crear port forwarding
netsh interface portproxy add v4tov4 listenport=3838 listenaddress=127.0.0.1 connectport=3838 connectaddress=169.254.76.190

# Verificar que se creó
netsh interface portproxy show all
```

Luego accede con localhost:
```
http://localhost:3838
```

### Para remover después (si quieres):
```powershell
netsh interface portproxy delete v4tov4 listenport=3838 listenaddress=127.0.0.1
```

---

## ✅ OPCIÓN 3: Usar Puerto Diferente (sin privilegios Admin)

Si no tienes permisos de administrador, usa un puerto alto (>8000):

### En tu terminal WSL (presiona Ctrl+C primero para detener):

```bash
cd /home/damo/scAnnex/dashboard && \
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)" && \
conda activate scannex-dashboard && \
export RETICULATE_PYTHON="$(which python3)" && \
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run" && \
echo "" && \
echo "══════════════════════════════════════════════════════════" && \
echo "  Dashboard Starting on PORT 8888..." && \
echo "══════════════════════════════════════════════════════════" && \
echo "" && \
echo "  Try these URLs:" && \
echo "     http://localhost:8888" && \
echo "     http://169.254.76.190:8888" && \
echo "" && \
echo "══════════════════════════════════════════════════════════" && \
echo "" && \
R -e "shiny::runApp('.', host='0.0.0.0', port=8888)"
```

Luego prueba:
- `http://localhost:8888`
- `http://169.254.76.190:8888`

---

## ✅ OPCIÓN 4: Desactivar Firewall Temporalmente (TESTING)

**SOLO PARA TESTING - No recomendado para uso permanente**

### En Windows:
1. Buscar "Firewall de Windows Defender"
2. Click "Activar o desactivar Firewall de Windows Defender"
3. Desactivar para "Redes privadas" temporalmente
4. Probar acceso
5. IMPORTANTE: Reactivar después

---

## ✅ OPCIÓN 5: Usar .wslconfig para Networking

Crear/editar archivo en Windows: `C:\Users\TU_USUARIO\.wslconfig`

```ini
[wsl2]
networkingMode=mirrored
firewall=false
```

Luego en PowerShell:
```powershell
wsl --shutdown
```

Reinicia WSL y vuelve a lanzar el dashboard.

---

## 🎯 MI RECOMENDACIÓN

Prueba en este orden:

### 1. PRIMERO: Opción 3 (Puerto 8888)
- No requiere admin
- Más probable que funcione
- Copia el comando de arriba

### 2. SI TIENES ADMIN: Opción 2 (Port Forwarding)
- Permite usar localhost
- Solución limpia
- Ejecuta el comando PowerShell

### 3. SI NADA FUNCIONA: Opción 5 (.wslconfig)
- Solución permanente
- Requiere reiniciar WSL

---

## 🔍 Verificar Estado Actual

En otra terminal WSL (sin cerrar el dashboard):

```bash
# Ver si el servidor está escuchando
ss -tlnp | grep 3838

# Test interno (desde WSL)
curl http://localhost:3838
# Si esto funciona, es problema de firewall Windows
```

---

## 📋 Comando Completo para Opción 3 (Puerto 8888)

```bash
# Detén el dashboard actual (Ctrl+C)
# Luego ejecuta:

cd /home/damo/scAnnex/dashboard && \
eval "$(/home/damo/miniforge3/bin/conda shell.bash hook)" && \
conda activate scannex-dashboard && \
export RETICULATE_PYTHON="$(which python3)" && \
export SCANNEX_DATA_PATH="/home/damo/scAnnex/results_slc_first_run" && \
pkill -f shiny 2>/dev/null && \
sleep 2 && \
echo "" && \
echo "══════════════════════════════════════════════════════════" && \
echo "  Dashboard Starting on PORT 8888..." && \
echo "══════════════════════════════════════════════════════════" && \
echo "" && \
echo "  🌐 Try in your Windows browser:" && \
echo "" && \
echo "     http://localhost:8888  (try first)" && \
echo "     http://127.0.0.1:8888  (if localhost fails)" && \
echo "" && \
echo "  Presiona Ctrl+C para detener" && \
echo "══════════════════════════════════════════════════════════" && \
echo "" && \
R -e "shiny::runApp('.', host='127.0.0.1', port=8888)"
```

**Nota:** Cambié `host='0.0.0.0'` a `host='127.0.0.1'` para evitar problemas de firewall.

---

## ✨ Resultado Esperado

Una vez que funcione, verás en el browser:

```
┌─────────────────────────────────────────────┐
│ scAnnex Dashboard                           │
├─────────────────────────────────────────────┤
│ [Data Input] [QC Overview] [Clustering]    │
├─────────────────────────────────────────────┤
│ Load your H5AD file here...                 │
└─────────────────────────────────────────────┘
```

---

¿Qué opción quieres probar primero? Te recomiendo la **Opción 3 (puerto 8888)** copiando el último comando.
