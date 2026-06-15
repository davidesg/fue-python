# fue ↔ fuf — Flujo de previsión y formato de archivos

Este documento describe el programa **fuf** (Forecast Using Fue), el formato de
archivo `.fuf.inp` que lo distingue del `.inp` estándar de fue, y la API Python
para generar y consumir previsiones mediante ese flujo.

---

## 1. Contexto histórico

En el sistema C original, la metodología Box-Jenkins-Treadway usa dos binarios
independientes:

| Binario | Rol |
|---------|-----|
| `fue`   | Estimación ML exacta → escribe `.out` + `.pre` |
| `fuf`   | Previsión L-pasos usando parámetros pre-estimados → `.out` |

El analista estima el modelo con `fue`, obtiene un `.pre` con los parámetros
estimados, lo convierte a formato fuf (añadiendo el horizonte y σ²) y lanza `fuf`
para generar el informe de previsión sin re-estimar.

En fue Python esta separación se mantiene intacta: `m.fit()` estima, `m.write_fuf()`
produce el archivo fuf, `fue.load_fuf()` lo lee, y `m.forecast_fuf()` calcula las
previsiones con los parámetros fijos del archivo.

---

## 2. Diferencia de formato: `.inp` fue vs `.inp` fuf

Ambos archivos usan la misma gramática (`**` como separador de secciones, datos
en líneas sin `**`). La única diferencia es una sección adicional que el archivo
fuf inserta **después de la línea de observaciones** y **antes de las variables
deterministas**:

```
** Forecast horizon and estimated innovation variance:
12  0.0000605632
```

La primera columna es el horizonte L (entero); la segunda es σ² (float de 10 dígitos).

El parser `_InpParser` detecta esta sección comprobando si la siguiente clave de
separador contiene `"forecast"`. Si está presente, almacena los valores en
`model._fuf_horizon` y `model._fuf_sigma2`. `fue.load_fuf()` exige que exista
(error si no); `fue.load()` la ignora silenciosamente si aparece.

### Otras diferencias visuales (no semánticas)

| Aspecto | `.inp` (fue) | `.inp` (fuf) |
|---------|-------------|-------------|
| Cabecera `*…*` | `Input file for program FUE` | `Input file for program FUF` |
| Parámetros | Pueden ser iniciales o estimados | Siempre estimados (fijos) |
| Sección fuf | Ausente | `** Forecast horizon and estimated...` |

El resto del formato (frecuencia, observaciones, det-vars, ARMA, BoxCox, ifadf,
refactor, datos) es idéntico en ambos.

---

## 3. API Python

### 3.1 Generar el archivo fuf desde un modelo estimado

```python
ts, m = fue.load("modelo.pre")
m.fit()
m.write_fuf(horizon=12, path="forecast_modelo.inp")
```

`Model.write_fuf(horizon, sigma2=None, path=None)`:
- Usa `model._result.sigma2` si el modelo está estimado.
- Acepta `sigma2` explícito para sobrescribir.
- `path` debe terminar en `.inp` (requisito de `load_fuf`).
- Sin `path`, devuelve el contenido como `str`.

### 3.2 Cargar y calcular previsiones

```python
ts_fuf, m_fuf = fue.load_fuf("forecast_modelo.inp")
fr = m_fuf.forecast_fuf()
```

`fue.load_fuf(path)`:
- Añade `.inp` al path si no termina en `.inp` ni `.pre`.
- Almacena `model._fuf_horizon` y `model._fuf_sigma2`.
- Devuelve el modelo sin ajustar (`_result is None`).

`Model.forecast_fuf(horizon=None, sigma2=None)`:
- Usa `model._fuf_horizon` / `model._fuf_sigma2` si no se pasan explícitamente.
- Llama a `eval_at_params(model)` (evaluación sin reoptimizar) para obtener
  residuos y sigma2 de los datos completos con los parámetros fijos.
- Sustituye `sigma2` del resultado con `model._fuf_sigma2` (σ² del archivo fuf),
  que es el calculado en el período de estimación — no el re-evaluado en la
  muestra ampliada (comportamiento deliberado para el seguimiento SPS).
- Construye `FitResult` sintético y lo asigna a `model._result` para que
  `write_forecast_report` pueda acceder a él.
- Devuelve `ForecastResult`.

### 3.3 `forecast()` vs `forecast_fuf()` — diferencia clave

| | `m.forecast(horizon)` | `m.forecast_fuf(horizon)` |
|---|---|---|
| Requiere | `m.fit()` previo (`_result` con sigma2 ML) | `_fuf_horizon` + `_fuf_sigma2` en el modelo |
| sigma2 usado | ML de la muestra de estimación | Fijado del archivo fuf (período de estimación) |
| Residuos | Del último ajuste MLE | `eval_at_params` con parámetros fijos |
| Uso correcto | Diagnóstico post-estimación | Informe de previsión y seguimiento SPS |

**Usar `forecast()` para los informes SPS es incorrecto**: el σ² MLE varía con
cada ampliación de muestra, mientras que el σ² fuf permanece constante (fijado
en el momento de la estimación) para que las bandas de previsión sean comparables
entre actualizaciones.

### 3.4 Generar el informe HTML

```python
from fue.report_forecast import write_forecast_report

write_forecast_report(m_fuf, fr, path="forecast_modelo.html",
                      title="IPC España", source="INE",
                      sps_name="Spain CPI")
```

Requiere que `m_fuf._result is not None` — satisfecho automáticamente por
`forecast_fuf()`.

---

## 4. Flujo completo (pipeline fuf)

```
fue.load(".pre")
    │
    ▼
m.fit()                         # estimación ML exacta
    │
    ▼
m.write_fuf(horizon, path)      # escribe .inp con sección fuf
    │
    ▼
fue.load_fuf(path)              # → (ts_fuf, m_fuf)
    │   almacena: _fuf_horizon, _fuf_sigma2
    ▼
m_fuf.forecast_fuf()            # eval_at_params + sigma2 fijo → ForecastResult
    │   almacena: m_fuf._result (sintético)
    ▼
write_forecast_report(...)      # HTML con figura 2-panel + tabla SPS
```

Para el seguimiento SPS (añadir nuevas observaciones):

```
fue.load_fuf(path_anterior)     # modelo + sigma2 fijo del período de estimación
    │
    ▼
construir ts_nuevo con datos ampliados
m_nuevo = fue.Model(ts_nuevo, **spec_del_modelo_anterior)
m_nuevo._fuf_sigma2 = sigma2_original   # mantener sigma2 de estimación
m_nuevo._fuf_horizon = horizon_original
    │
    ▼
m_nuevo.forecast_fuf()          # residuos en la muestra ampliada, sigma2 fijo
    │
    ▼
m_nuevo.write_fuf(...)          # actualizar el .inp fuf
write_forecast_report(...)      # nuevo informe con datos actualizados
```

---

## 5. `ForecastResult` — campos

```python
@dataclass
class ForecastResult:
    horizon: int
    level: np.ndarray           # previsión en nivel (unidades originales)  [L]
    level_std: np.ndarray       # desv. estándar del nivel (unidades orig.) [L]
    diff1: np.ndarray           # variación mensual/trimestral (%)          [L]
    diff1_std: np.ndarray
    seasonal_diff: np.ndarray   # variación anual (%)                       [L]
    seasonal_diff_std: np.ndarray
    sigma2: float               # varianza de innovaciones usada
```

`level_std` es **refactor-invariante**: `level_std[h] = sqrt(sigma2_1 · Σψ²_j)`
donde `sigma2_1` es la varianza en la escala original. El factor `refactor` se
cancela algebraicamente. Multiplicar por 100 da el error estándar en % directamente,
sin depender del valor de `refactor`:

```python
# Correcto para cualquier refactor:
std_pct = fr.level_std[h] * 100   # siempre en %
```

---

## 6. Escalado en el informe HTML (`refactor`)

`report_forecast.py` usa `model.refactor` para convertir las magnitudes internas
(escala transformada Box-Cox) a porcentajes. Las fórmulas correctas son:

| Cantidad | Fórmula en report_forecast.py | Nota |
|----------|-------------------------------|------|
| ERR tabla | `100 * residuals[i] / refactor` | residuos en escala ln·refactor → % |
| Std nivel | `fr.level_std[h] * 100` | level_std refactor-invariante |
| σ gráfico (banda ±2σ ERR) | `sqrt(fr.sigma2) * 100 / refactor` | coherente con ERR en % |

Estas fórmulas producen valores idénticos con `refactor=1` o `refactor=100`.
El analista puede usar cualquier valor de `refactor` sin afectar al informe.

---

## 7. CLI `fuf`

```bash
# Generar informe de previsión desde un archivo fuf
fuf forecast_modelo

# Con título, fuente y nombre SPS (genera .html además de .out)
fuf --title "IPC España" --source "INE" --sps "Spain" forecast_modelo
```

`fuf_cli.py` espera el archivo con extensión `.inp` (la añade si falta).
Escribe `forecast_modelo.out` (ASCII) y, si se pasa `--sps`, también
`forecast_modelo.html`.

---

## 8. Extensión `.inp` en rutas fuf

`fue.load_fuf(path)` añade `.inp` automáticamente si el path no termina en `.inp`
ni `.pre`. **El archivo escrito por `m.write_fuf(path=...)` debe terminar en `.inp`**
para que `load_fuf` lo encuentre sin transformación. Si se pasa un path sin
extensión (e.g. `"forecast_modelo"`), `write_fuf` escribe ese archivo exactamente,
pero `load_fuf("forecast_modelo")` buscará `"forecast_modelo.inp"` y fallará.

Patrón recomendado:

```python
fuf_path = "cases/serie/forecast_serie.inp"  # siempre terminar en .inp
m.write_fuf(horizon=12, path=fuf_path)
ts_fuf, m_fuf = fue.load_fuf(fuf_path)       # encuentra el archivo directamente
```
