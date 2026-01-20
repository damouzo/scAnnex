# GitHub Actions - Costos y Configuración

## 💰 Resumen de Costos

### Repositorio Público
- **GRATIS ilimitado** - Sin cargos nunca
- Recomendado para proyectos open source

### Repositorio Privado
- **2,000 minutos gratis/mes** (plan Free)
- Después: $0.008 por minuto
- Cada build: ~20 minutos

## ⚙️ Configuración Actual

El workflow está configurado para ejecutarse **SOLO**:

1. ✅ Cuando publicas un release (`git tag v1.0.0` + push)
2. ✅ Cuando lo activas manualmente (GitHub UI → Actions → Run workflow)
3. ❌ **NO** en cada push/commit (esto ahorra minutos)

## 📊 Uso Estimado

| Escenario | Minutos/Mes | Costo |
|-----------|-------------|-------|
| Desarrollo diario | 0 | $0 |
| Release mensual | ~20 | $0 (dentro de free tier) |
| 2 releases/mes | ~40 | $0 (dentro de free tier) |
| 10 releases/mes | ~200 | $0 (dentro de free tier) |

**Conclusión:** Muy por debajo del límite de 2000 min/mes

## 🎯 Estrategia Recomendada

### Durante Desarrollo (AHORA)
Usuarios usan **Conda environment** (método por defecto):
```bash
./setup_dashboard.sh  # Crea conda env (5-10 min, una vez)
./launch_dashboard.sh # Lanza dashboard
```
**Costo:** $0

### En Release Oficial (CUANDO PUBLIQUES)
GitHub Actions construye containers automáticamente:
```bash
git tag -a v1.0.0 -m "First release"
git push origin v1.0.0
```
→ Esperar ~20 min
→ Container disponible en GitHub Container Registry

Usuarios pueden:
- Opción A: Pull container pre-construido (~2 min)
- Opción B: Usar Conda (5 min, sigue funcionando)

**Costo:** ~$0 (1 release/mes = 20 min)

## 🔧 Opciones de Configuración

### Opción 1: Mantener Actual (RECOMENDADO) ✅
- Solo builds en releases
- Costo: ~$0/mes
- Ya configurado

### Opción 2: Desactivar Completamente
```bash
git mv .github/workflows/build-containers.yml \
       .github/workflows/build-containers.yml.disabled
```
- Solo usar Conda
- Costo: $0/mes
- Sin containers pre-construidos

### Opción 3: Build en Cada Push (NO RECOMENDADO para tesis)
Editar `.github/workflows/build-containers.yml`:
```yaml
on:
  push:
    branches: [ main ]
```
- Costo: ~$5-10/mes (depende de cuántos pushes)
- Útil solo para desarrollo muy activo con múltiples colaboradores

## 📈 Monitorear Uso

Ver uso actual:
1. GitHub → Settings → Billing and plans
2. Usage this month → Actions

## ❓ FAQ

**Q: ¿Qué pasa si me paso del límite?**
A: GitHub te cobra $0.008/min automáticamente. Te avisa al 75% y 90% del límite.

**Q: ¿Puedo establecer un límite de gasto?**
A: Sí, en Settings → Billing → Spending limit → Set limit ($5, $10, etc.)

**Q: ¿Repositorio público o privado?**
A: Si es público: GRATIS ilimitado. Si es privado: 2000 min/mes gratis.

**Q: ¿Vale la pena el costo?**
A: Para tu caso (1-2 releases/mes): No hay costo. Todo dentro del free tier.

## ✅ Recomendación Final

**Mantén la configuración actual:**
- Builds solo en releases
- Usuarios usan Conda por defecto
- Container pre-construido como opción alternativa
- Costo: $0 (bien dentro del límite gratuito)

**No necesitas cambiar nada más.**
