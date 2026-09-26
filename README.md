https://github.com/mramirezres/st-2026-2-tarea1-herramientas-ramirez_restrepo-manuela

# Tarea 1 

## Qué contiene el repositorio

| Archivo | Contenido | Dependencias |
|---|---|---|
| `R/00-lectura.R` | `leer_serie()`: lee un `ts` o un csv y devuelve un tibble con `t`, `fecha`, `y` | `tibble` |
| `R/01-graficos.R` | `graficar_serie()`, `correlograma()`, `grafico_pronostico()` | `ggplot2`, `patchwork` |
| `R/02-metodos.R` | Los ocho métodos (`ajustar_media`, `ajustar_mm`, `ajustar_ses`, `ajustar_dmm`, `ajustar_tendencia`, `ajustar_holt`), `optimizar()` y sus bloques de verificación | R base |
| `R/03-evaluacion.R` | `medidas()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()`, `validar_errores()` | `R/01-graficos.R` (usa `acf_manual()` y `correlograma()`) |
| `ejemplos/ejemplos.R` | Los nueve ejemplos (uno por método, más el contraejemplo), paso a paso | Todo lo anterior |
| `informe/informe.qmd` | El análisis escrito de cada ejemplo | `ejemplos/ejemplos.R` |
| `figs/` | Las figuras que genera `ejemplos.R` | — |

## Cómo se corre

Desde la raíz del repositorio, en R 4.3 o posterior:

```r
source("ejemplos/ejemplos.R")
```

No requiere ningún argumento ni archivo externo: las series vienen incluidas en R. Tiempo de ejecución observado: unos 20 segundos, en su mayoría el tiempo que toma guardar las 50 figuras en `figs/`.

## Cómo se usan las funciones

Ejemplo mínimo con la media móvil:

```r
source("R/00-lectura.R")
source("R/01-graficos.R")
source("R/02-metodos.R")

datos <- leer_serie(Nile, "Durbin y Koopman (2001)", "10^8 m^3")
ajuste <- ajustar_mm(datos$y, k = 9)

ajuste$yhat[1:12]        # pronosticos de un paso (NA en el calentamiento)
ajuste$pronosticar(5)    # 5 pronosticos extramuestrales
ajuste$parametros        # k usado y la ultima media movil
```

## Convenciones que fijan los números

- **Media simple:** recursiva, `Ŷ_{t+1}` es el promedio de `y_1, ..., y_t` (no el promedio de toda la muestra). `Ŷ_2 = Y_1`.
- **Media móvil (orden `k`):** `Ŷ_{t+1} = MM_t(k)`; los primeros `k` periodos no tienen pronóstico.
- **Suavizamiento exponencial simple:** `Ŷ_2 = Y_1`; luego `Ŷ_{t+1} = Ŷ_t + α(Y_t − Ŷ_t)`.
- **Doble media móvil (orden `k ≥ 2`):** hacen falta `2k − 1` observaciones antes del primer pronóstico.
- **Tendencias (lineal, cuadrática, exponencial):** ajustadas por mínimos cuadrados con toda la serie de estimación; `ŷ_t` es el valor ajustado de la curva, no un pronóstico de un paso (no hay calentamiento).
- **Holt:** `L_1 = Y_1`, `T̂_1 = 0`.
- **ACF (`acf_manual()`):** divisor único `T` en la autocovarianza y en la varianza, rezago 0 excluido.
- **Bandas del correlograma:** `±qnorm(0.975)/√n`, la misma fórmula que usa `plot.acf()` de R.

## Resumen de resultados

| Ejemplo | Serie | Parámetros | MASE método (validación) | MASE ingenuo (validación) |
|---|---|---|---|---|
| Media simple | discoveries | — | 0.919 | 1.137 |
| Media móvil | discoveries | k = 9 | 0.588 | 1.137 |
| SES | Nile | α = 0.24 | 0.806 | 0.836 |
| Doble media móvil | PetrolPrice | k = 2 | 0.643 | 1.219 |
| Tendencia lineal | austres | β0 = 12963.3, β1 = 50.71 | 3.974 | 6.091 |
| Tendencia cuadrática | airmiles | β0 = 1420.1, β1 = −473.9, β2 = 74.3 | 1.206 | 4.496 |
| Tendencia exponencial | JohnsonJohnson | a = −0.697, θ = 0.0428 | 4.140 | 9.458 |
| Holt | austres | α = 0.95, β = 0.80 | 1.484 | 6.091 |

*(Contraejemplo: Holt sobre lynx, MASE del método 2.203 contra 1.703 del ingenuo — el método no lo supera, como se espera en una serie estacionaria sin tendencia.)*

## Declaración de uso de IA

Usé Claude (Anthropic) como apoyo para escribir y depurar el código en R de este repositorio. 
La ayuda más importante fue para identificar qué casos debía validar cada función con stopifnot(): 
por ejemplo, que la tendencia exponencial solo acepta series con valores positivos (porque toma logaritmo), 
que ninguna función acepta un vector con NA, y que una ventana como k no puede ser mayor que el largo de la serie.

También me ayudó a entender y programar partes puntuales del código: cómo construir las fechas de una serie a partir 
de su frecuencia y su fecha de inicio, y cómo guardar cada gráfico en la carpeta figs/ con ggsave() para que quedaran 
ahí de forma automática al correr ejemplos.R. En la parte de optimización, me explicó qué significa que el valor óptimo 
caiga en el borde de la rejilla evaluada (por ejemplo, un alpha en el límite superior) y por qué eso es una señal de que 
convendría ampliar el rango de búsqueda. También me ayudó a interpretar la prueba de Durbin-Watson, entendiendo por qué 
no tiene un valor p como las demás pruebas y que, en cambio, su resultado se compara contra las cotas dL y dU de una tabla.

Además me ayudó a encontrar y corregir errores de programación (variables mal escritas, paréntesis y comas en el lugar 
equivocado). El diseño de los métodos, la elección de las series para cada ejemplo y el análisis final son míos.