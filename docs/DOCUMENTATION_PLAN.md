# Sistema de documentación de fue — propuesta

*12 de agosto de 2026. Redactada tras medir lo que hay, lo que existía antes del
porte y lo que el código promete sin poder demostrarlo.*

---

## 1. La tesis, que es lo que hay que poder decir en una frase

> **fue estima ARMAX por máxima verosimilitud exacta no condicionada, con
> estacionalidad resuelta frecuencia a frecuencia.**

Las dos mitades importan y ninguna es común:

- **ARMAX por ML exacta.** No condicionada, no Kalman con inicialización difusa.
  Con regresores deterministas que entran por funciones de transferencia
  racionales —impulso, escalón, rampa, armónicos, alternador, y definidos por el
  usuario— y transformación Box-Cox con retransformación de las previsiones.
- **Estacionalidad frecuencia a frecuencia.** La diferencia anual factorizada en
  sus factores de frecuencia fija, cada uno determinista o estocástico de forma
  independiente: la clase MEG/HSM. `statsmodels` no tiene nada equivalente —
  `SARIMAX` da estacionalidad multiplicativa entera o nada, y `VARMAX` ni eso.

Todo lo demás de la documentación cuelga de esa frase. Si un lector se va con
ella entendida y con un ejemplo corriendo, la documentación ha cumplido.

---

## 2. El objetivo declarado, y por qué cambia el diseño

> Ayudar a programadores, científicos y estadísticos a **evaluar y verificar** los
> algoritmos de fue.

Eso no es un tutorial con más páginas. Es una **cadena de custodia numérica**, y
el diseño entero se ordena alrededor de ella. Un manual de uso explica qué
teclear; esto tiene que permitir que alguien que no se fía compruebe.

### La cadena, medida — y es mucho mejor de lo que parecía

La procedencia de fue tiene cuatro eslabones, no tres, y el que yo daba por
ausente resulta estar documentado en el propio código fuente:

```
    algoritmo publicado  →  DRVUS  →  FUE  →  porte a Python
    Melard (1984) AS 197      Mauricio    Treadway lo       2026
    Mauricio (1995) JASA      2000        rebautiza
    Mauricio (1997) AS 311                («free»)
```

**DRVUS 1.01, de José Alberto Mauricio (2000)**, está en `SRC/drvus`, con
versiones previas en `SRC/drvus-source` (1.0 … 1.2.03). Su `readme` mapea cada
módulo a su publicación, y eso es la tabla de procedencia ya escrita por los
autores:

| módulo | qué calcula | fuente |
|---|---|---|
| `ELFVARMA.C` | log-verosimilitud de un VARMA(p,q) | **Mauricio (1995) JASA 90, 282-291** |
| `USMELARD.C` | log-verosimilitud de un ARMA(p,q) escalar | **Melard (1984) AS 197, 104-114** |
| `QNEWTOPT.C` | optimización BFGS quasi-Newton factorizada | |
| `NLATOOLS.C` | álgebra lineal numérica y memoria dinámica | (NR) |
| `DRVMLEST.C` | driver de estimación por MV exacta | |
| `DIAGNOSE.C` | estadísticos muestrales y gráficos | |

### Y el núcleo numérico es el de Mauricio, verbatim

Comparados los ficheros de `drvus/src` con la copia empotrada en
`fue/csrc/internal`, tras normalizar los finales de línea (los de Mauricio son
CRLF y latin-1; los de fue, LF y UTF-8):

| fichero | líneas distintas | qué son |
|---|---|---|
| `elfvarma.c` | **23** | cabecera GPL (14), «José» en UTF-8, el `#include`, **un cambio funcional**, y un comentario de cierre |
| `usmelard.c` | **21** | cabecera GPL, codificación, `#include`. **Ningún cambio funcional** |
| `drvmlest.c` | **22** | ídem. **Ningún cambio funcional** |
| `nlatools.c` | 1129 | limpieza aparte y GSL — es el módulo que sí se reescribió |

El único cambio funcional en todo el núcleo de la verosimilitud es, en
`elfvarma.c` línea 513:

```c
-    eigenqr( a, n, wr, wi );
+    if ( n>1 ) {gsl_eigenqr( a, n, wr, wi );}
```

La rutina de autovalores de Numerical Recipes sustituida por la de GSL, con
guarda para `n>1`. **Eso es todo.** El código que implementa AS 311 y AS 197 es
el que Mauricio publicó, letra por letra.

**Consecuencia de diseño, y es la buena noticia de este documento:** la
verificación no hay que construirla de cero. Hay que **exhibirla** — y una
afirmación de procedencia que se puede comprobar con un `diff` de tres líneas
vale más que cualquier párrafo de prosa.

Lo que sí falta, y sigue faltando:

- **Contrastar contra los valores publicados en los propios artículos.** AS 311 y
  AS 197 traen casos de prueba con datos y resultados. Ningún test los ejecuta:
  `test_estimation.py` usa *«reference values obtained from the reference
  binary»*, así que verifica porte-contra-binario, no binario-contra-artículo.
  Con el código siendo el de Mauricio esto es casi una formalidad — pero es la
  formalidad que convierte «confía» en «comprueba».
- **El `diff` como test.** Que el núcleo siga siendo el de Mauricio es una
  propiedad que puede romperse en silencio, exactamente como se rompió la
  sincronía de las copias del C en drvarma (`drvarma/bugs/BUG-0002`). Un test que
  corra ese `diff` contra `drvus/src` con la lista de excepciones declarada
  convierte la afirmación en un invariante.

### La capa 1.0: el manual primitivo, y lo que resuelve

`drvus-source/1.0` es el original, y trae dos documentos que hoy no lee nadie:

**`README.txt` (2002)** añade un módulo que el de 1.01 ya no menciona:

| módulo | fuente | estado en 1.0 |
|---|---|---|
| `MULTSHEA.C` | **Shea, B.L. (1989) ALGORITHM AS 242, 161-184** | *no implementado* |

Ese mismo `multshea.c` sigue compilado y sin conectar en `drvarma` hoy, bajo el
epígrafe «Out of scope for this port — Shea (AS 242)». **Veintiséis años de la
misma decisión, tomada dos veces sin que la segunda supiera de la primera.**

**`ABTreadway-Drvus.doc` (SMP 3/2001, con notas de ABT 5/2001)** es el manual
primitivo. Describe el `.inp` parte por parte —intervalo muestral, deterministas,
estructura ARMA, operadores de frecuencia fija, media, transformaciones, datos— y
cierra con dos secciones que valen más que el resto:

*Comentarios finales*, punto 2 — el contrato del parser, dicho en 2001:

> *El programa no interpreta los comentarios escritos entre asteriscos, solamente
> lee los números que espera encontrar en la posición correcta. Esto quiere decir
> que si se cometiera el error de introducir una línea en blanco, por ejemplo, al
> comienzo del fichero de entrada, el programa bloquearía al ordenador.*

Un parser **posicional**, que ignora los comentarios y no tolera desviaciones.
Eso es la especificación del formato y hay que citarla literalmente en §4.

*NOTAS ADICIONALES, por ABT 5/2001* — tres consejos, y el segundo es el origen de
una convención que atraviesa toda la suite:

1. Imprimir siempre el fichero de entrada antes del de salida en el mismo
   listado, para localizar los errores de entrada.
2. > *Es deseable que la norma del gradiente sea cero hasta toda la precisión que
   > ofrece el programa. Cuando esto no ocurre, aunque se obtenga el mensaje de
   > convergencia […] es muchas veces útil **escalar los datos** para que tengan
   > más precisión de entrada. Esto se hace multiplicando todos los valores de la
   > variable lnY […] por, p.e., **100** antes de introducirlos como entrada al
   > programa. Por supuesto, también se multiplicará cada parámetro de
   > intervención por el mismo factor, y la salida presentará una sigma
   > multiplicado por el mismo factor.*
3. Usar otro programa para factorizar operadores.

**El consejo (2) es el origen documentado de `refactor=100`.** La convención que
recorre `_RESCALE_FACTOR`, `docs/RESCALING_ARCHITECTURE.md` de art, y los
defectos BUG-0001 (μ colapsa al reescalar) y BUG-0007 (μ leída 100× fuera de
escala) **no es una decisión de modelización: es un remedio de precisión numérica
del optimizador**, escrito por Treadway en 2001 y aplicado desde entonces sin que
su porqué estuviera en ningún sitio legible.

Documentar eso cambia cómo se lee el reescalado: no es una rareza heredada, es
una respuesta a la norma del gradiente. Y explica por qué los dos defectos que
produjo eran de *escala en el informe*, no de estimación.

El consejo (3) también tiene descendencia: es el programa `Root`, que acabó en
`art/roots.py`.

Y `1.0/Drvus/M2.inp` más el `S.inp` del manual son los **ficheros de ejemplo
originales**. El del manual es un ARMAX con dos rampas consecutivas, exactamente
el caso que el ejemplo 3 del §5 tiene que enseñar. Se recuperan.

---

## 3. Punto de partida, medido

| | líneas | dónde |
|---|---|---|
| `README.md` | 113 | la página de PyPI **es** esto |
| `FUF.md` | 244 | previsión |
| `PERFORMANCE.md` | 280 | arquitectura y rendimiento, en español |
| `docs/` | **0** | el directorio existe y está vacío |
| **total fue** | **637** | |
| **total art** | **2.825** | 11 documentos |

fue tiene el **22 %** de la documentación de art, siendo la pieza sobre la que
art, drtran y drvarma se apoyan.

Y existe material previo que no se ha aprovechado:

- **Manual de Treadway y Guerrero**, `atws/manuales/fue/` — 61 KB de LaTeX en
  español, latin-1, con índice ya bueno: *Lo que hace FUE · Cómo aprender a usar
  FUE · Ejecución · Nociones de series temporales (componente determinista /
  estocástico) · Especificación del `.inp` · Contenido del `.out` · Contrastes
  formales de hipótesis · Otros procedimientos*. La sección de contrastes cubre
  no estacionariedad, no invertibilidad, frecuencia fija de un AR(2),
  estacionalidad estocástica y simplificación.
- **Manual de FUF**, `atws/manuales/fuf/` — versión 1.06.2.
- **El manual primitivo de DRVUS**, `drvus-source/1.0/Drvus/ABTreadway-Drvus.doc`
  — Word de 2001, ilegible sin convertir, con el contrato del parser y el origen
  del reescalado. Ver arriba.
- **`biblio.bib`** con las referencias ya reunidas: Box-Cox 1964, Box-Jenkins(-Reinsel),
  Mauricio 1992/1995/1997, Melard 1984 (AS 197), Gallego 1995/1996, Hylleberg,
  Hannan, Cleveland-Tiao, Bell-Hillmer, Treadway 1994.

**No se reescribe lo de Treadway: se traduce, se actualiza al porte y se le
añade lo que el porte trajo.** Su índice es mejor punto de partida que uno nuevo.

---

## 4. La estructura propuesta

Ocho secciones. El orden es el de un lector que llega sin conocer nada y acaba
auditando el código.

### 0 · Portada — «qué es y por qué es distinto»
Una página. La frase del §1, la ecuación que fue ajusta de verdad, un ejemplo de
diez líneas que corre, y la tabla de comparación con `statsmodels` (existe ya en
`atsw-suite/docs/COMPARISON_STATSMODELS.md`; se enlaza, no se duplica).

### 1 · Empezar
Instalación —incluida la distinción rueda compilada / respaldo en Python puro y
cómo saber cuál se está usando—, el primer modelo, y la lectura del `.out`.
Rescata «Cómo aprender a usar FUE» de Treadway.

### 2 · El modelo, formalmente
La clase de modelos, escrita sin ambigüedad: `y = D + N`, aditiva **en el
nivel**; el componente determinista y sus tipos; el estocástico
ARIMA(p,d,q)(P,D,Q)ₛ y la forma MEG; Box-Cox y la retransformación; y la
factorización de la diferencia anual en factores de frecuencia con la tabla de
Abraham–Box. Es la traducción de §«Nociones de series temporales» de Treadway,
con la notación unificada con la del paper SF_MEG.

### 3 · **Verificación** ← el capítulo que justifica el resto
Cuatro partes, y las tres primeras hay que **construirlas**:

1. **Procedencia por algoritmo.** Una tabla: qué función implementa qué
   algoritmo, de qué publicación, y qué test lo comprueba.

   | función | algoritmo | fuente | verificado por |
   |---|---|---|---|
   | `elfvarma.elf_scalar` | verosimilitud exacta VARMA | Mauricio (1995, JASA 90, 282-291); AS 311 (1997) | C: idéntico a `drvus/src/elfvarma.c` salvo GSL. Python: *(pendiente: caso publicado)* |
   | `elfvarma.flikam_scalar` | verosimilitud ARMA rápida | Melard (1984), AS 197, 104-114 | C: **idéntico** a `drvus/src/usmelard.c`. Python: *(pendiente)* |
   | `cast_us.cast_us_py` | forma de innovaciones | Ansley (1979, 1982) | `test_cast_us.py` |
   | `qnewtopt` | quasi-Newton BFGS | Dennis & Schnabel | `test_qnewtopt.py` |

   Las casillas «pendiente» son honestas y son el trabajo.

2. **Contra valores publicados.** Los casos de prueba de AS 311 y AS 197 vienen
   con datos y resultados en el propio artículo: se ejecutan y se comparan con
   tolerancia declarada. Esto es lo que cierra el eslabón que falta.

3. **C contra Python.** Ya existe y está bien; lo que falta es **publicarlo como
   resultado**, con la tolerancia (~1e-11) y el conjunto de casos, en vez de
   dejarlo dentro de la carpeta de tests.

4. **Los límites conocidos, dichos.** El salto del perfil en la frontera de
   invertibilidad (SF_MEG, apéndice), el óptimo espurio según el build
   (`bugs/BUG-0005`), y `converged=True` sin diagnóstico sobre ajustes absurdos.
   Documentar los límites es parte de permitir verificar; ocultarlos lo impide.

### 4 · El formato `.inp` / `.out` / `.pre`
El contrato, con gramática completa y campo a campo. Es la sección que más se
consulta y la que hace posible auditar: los artefactos son texto inspeccionable.
Incluye el invariante del `.pre` —correr fue sobre un `.pre` no mueve los
números— como propiedad comprobable, no como afirmación.

### 5 · Contrastes formales
Traducción y actualización de §«Contrastes formales de hipótesis» de Treadway:
no estacionariedad, no invertibilidad (DCD), frecuencia fija de un AR(2),
estacionalidad estocástica (MEG) y simplificación. Con los valores críticos del
paper SF_MEG donde ahora los hay, y diciendo cuáles siguen siendo interpolados.

### 6 · El porte
Lo que David pide explícitamente y no existe en ningún sitio: qué se hizo, qué
cambió y qué no. Motor C empotrado (`csrc/`) frente al respaldo en Python puro
(`elfvarma.py`, `cast_us.py`, `qnewtopt.py`: 1.895 líneas que reimplementan el
motor); qué es idéntico y qué se decidió distinto; los defectos que el porte
encontró en el original; el reescalado; y las ruedas y su CI. La bitácora
`PORTE.md` de drtran es el modelo — es la mejor pieza documental de la suite.

### 7 · Migración desde el FUE en C
Para los usuarios de Treadway: equivalencias de línea de órdenes, qué
`.inp` se leen tal cual —y el que no, por la codificación latin-1,
`bugs/BUG-0010`—, y qué salidas difieren.

### 8 · Referencia de la API
Generada de los docstrings. 24 símbolos públicos; hoy no hay referencia.

---

## 5. Ejemplos, graduados

Cada uno corre solo, con datos que viajan en el paquete, y se prueba en CI.

| # | ejemplo | enseña |
|---|---|---|
| 1 | ARIMA(1,1,0) sobre un índice de precios | el flujo mínimo: cargar, ajustar, leer |
| 2 | **Airline (0,1,1)(0,1,1)₁₂** | el canónico de Box-Jenkins con estacionalidad |
| 3 | **ARMAX**: intervenciones escalón e impulso con fecha | la mitad «X», que es la tesis |
| 4 | **Armónicos deterministas frente a ∇₁₂** | la otra mitad: estacionalidad por frecuencias |
| 5 | **Un MEG mixto**: unas frecuencias deterministas y otras estocásticas | la clase que sólo fue tiene |
| 6 | Box-Cox y retransformación de previsiones | la banda asimétrica, que se lee mal a menudo |
| 7 | Previsión con `fuf` y origen anclado | `-estwin`, parámetros fijos |
| 8 | **Reproducir un caso publicado de AS 311** | el ejemplo que ES verificación |
| 9 | **El `S.inp` original del manual de 2001** | ARMAX con dos rampas consecutivas — el ejemplo que Mauricio y Treadway eligieron para explicar el programa. Vale como continuidad y como prueba de que los `.inp` de 2001 siguen leyéndose |

Los ejemplos 3, 4 y 5 son la razón de ser del paquete y hoy no están escritos en
ningún sitio.

---

## 6. Publicación

- **PyPI**: el README pasa a ser una portada de verdad —la tesis, el ejemplo de
  diez líneas, el enlace al sitio— en vez de la documentación entera.
- **GitHub Pages con mkdocs**, como `atsw-suite`. La infraestructura ya está
  resuelta ahí, incluido el escollo de que mkdocs resuelve `-d` respecto al
  fichero de configuración y no al directorio de trabajo.
- **Enlazado desde `atsw`**, para que el paraguas lleve al motor.
- Y los PDF de Treadway **accesibles tal cual**, como material histórico: son la
  especificación original y tienen valor de archivo.

---

## 7. Orden y esfuerzo

Por dependencia, no por tamaño.

| | qué | por qué primero | esfuerzo |
|---|---|---|---|
| 1 | §3.1 tabla de procedencia | Es un inventario: se hace leyendo el código, y **decide todo lo demás** porque enseña dónde están los huecos | 1 día |
| 2 | §4 el contrato de ficheros | Es lo más consultado y ya está casi escrito en el manual de Treadway | 1 día |
| 3 | §3.2 casos publicados de AS 311 / AS 197 | Cierra el eslabón que falta. Es código, no prosa | 2-3 días |
| 4 | §2 el modelo + §5 contrastes | Traducción y actualización de Treadway | 2 días |
| 5 | Ejemplos 1-5 | Los tres del medio son la tesis | 2 días |
| 6 | §6 el porte | Sólo lo puede escribir quien lo hizo | 1-2 días |
| 7 | §0/§1 portada y arranque | Se escriben al final: se sabe qué prometer cuando el resto existe | medio día |
| 8 | §8 API generada, §7 migración, publicación | mecánico | 1 día |

Del orden de dos semanas de trabajo. La pieza 3 es la que más valor añade por
día y la única que no se puede escribir «con lo que ya se sabe».

---

## 8. Idioma: inglés — decidido

**Inglés**, y hay un argumento mejor que el alcance: **el código original ya está
en inglés.** Los comentarios de `elfvarma.c`, `usmelard.c` y `drvmlest.c` —los de
Mauricio, de 1995 y 1996— están en inglés, igual que el `readme` de DRVUS
describe los módulos en inglés. Las versiones iniciales lo estaban también. Lo
que está en español es el manual de usuario de Treadway, que es una capa
posterior.

Documentar en inglés no es traducir el proyecto: es volver al idioma en el que se
escribió el núcleo. Los manuales de Treadway se enlazan como material histórico,
en su idioma, que es donde tienen valor de archivo.
