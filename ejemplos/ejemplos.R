for (archivo in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  source(archivo)
}
dir.create("figs", showWarnings = FALSE)

fuente_help <- function(tema) {
  ayuda <- do.call(help, list(tema, package = "datasets"))
  texto <- capture.output(tools::Rd2txt(utils:::.getHelpFile(ayuda),
                                        options = list(underline_titles = FALSE)))
  inicio <- which(texto == "Source:")
  encabezados <- which(grepl("^[A-Z][A-Za-z ]*:$", texto))
  fin <- min(c(encabezados[encabezados > inicio], length(texto) + 1)) - 1
  fuente <- paste(texto[(inicio + 1):fin], collapse = " ")
  fuente <- gsub("_", "", fuente)
  trimws(gsub("[[:space:]]+", " ", fuente))
}


ingenuo <- function(y_est, h, s = 1) {
  n <- length(y_est)
  if (s == 1) {
    return(rep(y_est[n], h))
  }
  y_est[n - s + ((1:h - 1) %% s) + 1]
}


mostrar_prueba <- function(nombre, h0, h1, distribucion, res, nota = NULL) {
  if (is.na(res$valor_p)) {
    decision <- "se compara con las cotas dL y dU de la tabla"
  } else if (res$valor_p < 0.05) {
    decision <- "se rechaza H0 al 5%"
  } else {
    decision <- "no se rechaza H0 al 5%"
  }
  critico <- if (is.na(res$valor_critico)) "no aplica" else format(res$valor_critico, digits = 5)
  valor_p <- if (is.na(res$valor_p)) "no aplica" else format(res$valor_p, digits = 3)
  cat("\n  [", nombre, "]\n", sep = "")
  cat("    H0:", h0, "\n    H1:", h1, "\n")
  cat("    Estadistico:", format(res$estadistico, digits = 5), "| Distribucion:", distribucion, "\n")
  cat("    Valor critico 5%:", critico, "| Valor p:", valor_p, "\n")
  cat("    Decision:", decision, "\n")
  if (!is.null(nota)) cat("    Nota:", nota, "\n")
}

# --------------
# Ejemplo 1

cat("\n\n", strrep("=", 30), "\nEJEMPLO 1: Media simple\n", strrep("=", 30), "\n", sep = "")

fuente <- fuente_help("discoveries")
datos <- leer_serie(discoveries, fuente, "descubrimiento por año")
n <- nrow(datos)
y <- datos$y

stopifnot(n == length(discoveries), all(y == as.numeric(discoveries)))

cat("\n1. DESCRIPCION\n")
cat("  Nombre en R: discoveries\n  Fuente:", fuente, "\n  Unidad: descubrimientos por año\n")
cat("  Frecuencia:", attr(datos, "frecuencia"), "\n")
cat("  Fechas:", format(min(datos$fecha)), "a", format(max(datos$fecha)), "\n")
cat("  Observaciones:", n, "\n")


## Graficos, contraste y Ljung-Box

cat("\n2. PATRON DE LA SERIE\n")

ggsave("figs/e1-serie.png", graficar_serie(datos, "Descubrimientos importantes por anio"), width = 8, height = 4.5)

ggsave("figs/e1-correlograma.png", correlograma(datos), width = 8, height = 6)

m <- min(floor(n / 4), 24)
r <- acf_manual(y, m)
limite <- qnorm(0.975) / sqrt(n)

contraste <- data.frame(rezago = 1:m, rh = r, z = r * sqrt(n), fuera_de_banda = abs(r) > limite)

print(data.frame(banda = limite, r1 = r[1], rezagos_fuera = sum(contraste$fuera_de_banda),
                 esperados_por_azar = 0.05 * m), row.names = FALSE)

print(contraste[contraste$fuera_de_banda, ], row.names = FALSE)

### Ljung-Box

lb_serie <- ljung_box(r, n, m, 0)

mostrar_prueba("Ljung-Box sobre la serie (p = 0)",
               "todas las autocorrelaciones hasta m son cero", "alguna es distinta de cero",
               paste0("chi2 con ", lb_serie$gl, " g.l."), lb_serie)

### Partición

h <- min(12, floor(0.2 * n))
n_est <- n - h
y_est <- y[1:n_est]
y_val <- y[(n_est + 1):n]

print(data.frame(h = h, estimacion = n_est, validacion = h), row.names = FALSE)

## Metodo: media simple

cat("\n4. AJUSTE (media)\n")

ajuste <- ajustar_media(y_est)

print(data.frame(media = ajuste$parametros$media), row.names = FALSE)

p <- 1
pron <- ajuste$pronosticar(h)
pron_ing <- ingenuo(y_est, h)

tabla <- rbind(
  data.frame(caso = "Metodo, un paso (estimacion)", medidas(y_est, ajuste$yhat, y_est)),
  data.frame(caso = "Metodo, h pasos (validacion)", medidas(y_val, pron, y_est)),
  data.frame(caso = "Ingenuo, h pasos (validacion)", medidas(y_val, pron_ing, y_est)))

print(tabla, row.names = FALSE, digits = 4)

### Validación errores

cat("\n6. ERRORES DE UN PASO\n")

errores <- y_est - ajuste$yhat
ve <- validar_errores(errores, p = p, fechas = datos$fecha[1:n_est])

ggsave("figs/e1-errores.png", ve$grafico_errores, width = 8, height = 3.5)
ggsave("figs/e1-correlograma-errores.png", ve$correlograma, width = 8, height = 6)

cat("  Numero de errores:", ve$n, "| p =", p, "\n")

mostrar_prueba("t de media cero", "E[e] = 0", "E[e] distinto de 0",
               paste0("t con ", ve$prueba_t$gl, " g.l."), ve$prueba_t)

mostrar_prueba("Ljung-Box sobre los errores", "no hay autocorrelacion hasta m", "hay autocorrelacion",
               paste0("chi2 con ", ve$ljung_box$gl, " g.l. (m - p)"), ve$ljung_box)

mostrar_prueba("Jarque-Bera", "los errores son normales", "los errores no son normales",
               "chi2 con 2 g.l. (asintotica)", ve$jarque_bera,
               nota = if (ve$n < 20) paste("solo hay", ve$n, "errores (< 20): la decision es poco fiable"))

mostrar_prueba("Durbin-Watson", "no hay autocorrelacion de primer orden", "hay autocorrelacion de primer orden",
               "d entre 0 y 4", ve$durbin_watson)

ggsave("figs/e1-final.png",
       grafico_pronostico(datos, n_est, ajuste$yhat, pron, "Descubrimientos importantes por año (media)", ingenuo = pron_ing),
       width = 8, height = 4.5)

# ----------------------------------
# Ejemplo 2

cat("\n\n", strrep("=", 30), "\nEJEMPLO 2: Media movil\n", strrep("=", 30), "\n", sep = "")

fuente <- fuente_help("discoveries")
datos <- leer_serie(discoveries, fuente, "descubrimientos por año")
n <- nrow(datos)
y <- datos$y

stopifnot(n == length(discoveries), all(y == as.numeric(discoveries)))

cat("\n1. DESCRIPCION\n")
cat("  Nombre en R: discoveries\n  Fuente:", fuente, "\n  Unidad: descubrimientos por año\n")
cat("  Frecuencia:", attr(datos, "frecuencia"), "\n")
cat("  Fechas:", format(min(datos$fecha)), "a", format(max(datos$fecha)), "\n")
cat("  Observaciones:", n, "\n")

## Graficos, contraste y Ljung-Box

cat("\n2. PATRON DE LA SERIE\n")

ggsave("figs/e2-serie.png", graficar_serie(datos, "Descubrimientos importantes por año"), width = 8, height = 4.5)

ggsave("figs/e2-correlograma.png", correlograma(datos), width = 8, height = 6)

m <- min(floor(n / 4), 24)
r <- acf_manual(y, m)
limite <- qnorm(0.975) / sqrt(n)

contraste <- data.frame(rezago = 1:m, rh = r, z = r * sqrt(n), fuera_de_banda = abs(r) > limite)

print(data.frame(banda = limite, r1 = r[1], rezagos_fuera = sum(contraste$fuera_de_banda),
                 esperados_por_azar = 0.05 * m), row.names = FALSE)

print(contraste[contraste$fuera_de_banda, ], row.names = FALSE)

### Ljung-Box

lb_serie <- ljung_box(r, n, m, 0)

mostrar_prueba("Ljung-Box sobre la serie (p = 0)",
               "todas las autocorrelaciones hasta m son cero", "alguna es distinta de cero",
               paste0("chi2 con ", lb_serie$gl, " g.l."), lb_serie)

### Partición

h <- min(12, floor(0.2 * n))
n_est <- n - h
y_est <- y[1:n_est]
y_val <- y[(n_est + 1):n]

print(data.frame(h = h, estimacion = n_est, validacion = h), row.names = FALSE)

## Metodo: Media movil

cat("\n4. AJUSTE (mm)\n")

### Ajuste sobre la estimación

opt <- optimizar(y_est, "mm")

print(cbind(opt$optimo, en_borde = opt$en_borde), row.names = FALSE)

ggsave("figs/e2-optimizacion.png", grafico_mse(opt, "MSE de un paso - Descubrimientos importantes por año"), width = 7, height = 5)

interpretar_borde(opt)

ajuste <- ajustar_mm(y_est, opt$optimo$k)
p <- 1 

pron <- ajuste$pronosticar(h)
pron_ing <- ingenuo(y_est, h)

### Medidas y comparación

tabla <- rbind(
  data.frame(caso = "Metodo, un paso (estimacion)", medidas(y_est, ajuste$yhat, y_est)),
  data.frame(caso = "Metodo, h pasos (validacion)", medidas(y_val, pron, y_est)),
  data.frame(caso = "Ingenuo, h pasos (validacion)", medidas(y_val, pron_ing, y_est)))
print(tabla, row.names = FALSE, digits = 4)

### Validación errores

cat("\n6. ERRORES DE UN PASO\n")

errores <- y_est - ajuste$yhat
ve <- validar_errores(errores, p = p, fechas = datos$fecha[1:n_est])

ggsave("figs/e2-errores.png", ve$grafico_errores, width = 8, height = 3.5)
ggsave("figs/e2-correlograma-errores.png", ve$correlograma, width = 8, height = 6)

cat("  Numero de errores:", ve$n, "| p =", p, "\n")

mostrar_prueba("t de media cero", "E[e] = 0", "E[e] distinto de 0",
               paste0("t con ", ve$prueba_t$gl, " g.l."), ve$prueba_t)

mostrar_prueba("Ljung-Box sobre los errores", "no hay autocorrelacion hasta m", "hay autocorrelacion",
               paste0("chi2 con ", ve$ljung_box$gl, " g.l. (m - p)"), ve$ljung_box)

mostrar_prueba("Jarque-Bera", "los errores son normales", "los errores no son normales",
               "chi2 con 2 g.l. (asintotica)", ve$jarque_bera,
               nota = if (ve$n < 20) paste("solo hay", ve$n, "errores (< 20): la decision es poco fiable"))

mostrar_prueba("Durbin-Watson", "no hay autocorrelacion de primer orden", "hay autocorrelacion de primer orden",
               "d entre 0 y 4", ve$durbin_watson)

ggsave("figs/e2-final.png",
       grafico_pronostico(datos, n_est, ajuste$yhat, pron, "Descubrimientos importantes por año (mm)", ingenuo = pron_ing),
       width = 8, height = 4.5)

# -----------
# Ejemplo 3

cat("\n\n", strrep("=", 40), "\nEJEMPLO 3: Suavizamiento exponencial simple\n", strrep("=", 40), "\n", sep = "")

fuente <- fuente_help("Nile")
datos <- leer_serie(Nile, fuente, "10^8 m^3")
n <- nrow(datos)
y <- datos$y

stopifnot(n == length(Nile), all(y == as.numeric(Nile)))

cat("\n1. DESCRIPCION\n")
cat("  Nombre en R: Nile\n  Fuente:", fuente, "\n  Unidad: 10^8 m^3\n")
cat("  Frecuencia:", attr(datos, "frecuencia"), "\n")
cat("  Fechas:", format(min(datos$fecha)), "a", format(max(datos$fecha)), "\n")
cat("  Observaciones:", n, "\n")

# Graficos, contrastes y Ljung-Box

cat("\n2. PATRON DE LA SERIE\n")

ggsave("figs/e3-serie.png", graficar_serie(datos, "Caudal anual del Nilo"), width = 8, height = 4.5)
ggsave("figs/e3-correlograma.png", correlograma(datos), width = 8, height = 6)

m <- min(floor(n / 4), 24)
r <- acf_manual(y, m)
limite <- qnorm(0.975) / sqrt(n)

contraste <- data.frame(rezago = 1:m, rh = r, z = r * sqrt(n), fuera_de_banda = abs(r) > limite)

print(data.frame(banda = limite, r1 = r[1], rezagos_fuera = sum(contraste$fuera_de_banda),
                 esperados_por_azar = 0.05 * m), row.names = FALSE)

print(contraste[contraste$fuera_de_banda, ], row.names = FALSE)

### Ljung-Box

lb_serie <- ljung_box(r, n, m, 0)

mostrar_prueba("Ljung-Box sobre la serie (p = 0)",
               "todas las autocorrelaciones hasta m son cero", "alguna es distinta de cero",
               paste0("chi2 con ", lb_serie$gl, " g.l."), lb_serie)

### Partición 

h <- min(12, floor(0.2 * n))
n_est <- n - h
y_est <- y[1:n_est]
y_val <- y[(n_est + 1):n]

print(data.frame(h = h, estimacion = n_est, validacion = h), row.names = FALSE)

# Metodo: SES

cat("\n4. AJUSTE (ses)\n")

## Ajuste sobre la estimación

opt <- optimizar(y_est, "ses")

print(cbind(opt$optimo, en_borde = opt$en_borde), row.names = FALSE)

ggsave("figs/e3-optimizacion.png", grafico_mse(opt, "MSE de un paso - Caudal anual del Nilo"), width = 7, height = 5)

interpretar_borde(opt)

ajuste <- ajustar_ses(y_est, opt$optimo$alpha)

p <- 1

pron <- ajuste$pronosticar(h)
pron_ing <- ingenuo(y_est, h)

## Medidas y comparacion

tabla <- rbind(
  data.frame(caso = "Metodo, un paso (estimacion)", medidas(y_est, ajuste$yhat, y_est)),
  data.frame(caso = "Metodo, h pasos (validacion)", medidas(y_val, pron, y_est)),
  data.frame(caso = "Ingenuo, h pasos (validacion)", medidas(y_val, pron_ing, y_est)))

print(tabla, row.names = FALSE, digits = 4)

## Validacion errores

cat("\n6. ERRORES DE UN PASO\n")

errores <- y_est - ajuste$yhat
ve <- validar_errores(errores, p = p, fechas = datos$fecha[1:n_est])

ggsave("figs/e3-errores.png", ve$grafico_errores, width = 8, height = 3.5)
ggsave("figs/e3-correlograma-errores.png", ve$correlograma, width = 8, height = 6)

cat("  Numero de errores:", ve$n, "| p =", p, "\n")

mostrar_prueba("t de media cero", "E[e] = 0", "E[e] distinto de 0",
               paste0("t con ", ve$prueba_t$gl, " g.l."), ve$prueba_t)

mostrar_prueba("Ljung-Box sobre los errores", "no hay autocorrelacion hasta m", "hay autocorrelacion",
               paste0("chi2 con ", ve$ljung_box$gl, " g.l. (m - p)"), ve$ljung_box)

mostrar_prueba("Jarque-Bera", "los errores son normales", "los errores no son normales",
               "chi2 con 2 g.l. (asintotica)", ve$jarque_bera,
               nota = if (ve$n < 20) paste("solo hay", ve$n, "errores (< 20): la decision es poco fiable"))

mostrar_prueba("Durbin-Watson", "no hay autocorrelacion de primer orden", "hay autocorrelacion de primer orden",
               "d entre 0 y 4", ve$durbin_watson)

ggsave("figs/e3-final.png",
       grafico_pronostico(datos, n_est, ajuste$yhat, pron, "Caudal anual del Nilo (ses)", ingenuo = pron_ing),
       width = 8, height = 4.5)

