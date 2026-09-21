library(ggplot2)
library(patchwork)

graficar_serie <- function(datos, titulo){
  fuente <- attr(datos, "fuente")
  unidad <- attr(datos, "unidad")
  
  # En el caso de tener una fuente larga
  fuente_texto <- paste(strwrap(paste0("Fuente: ", fuente),
                                width = 90), collapse = "\n")
  caption <- paste0(fuente_texto, "\nObservaciones: ", nrow(datos))
  
  anios <- as.numeric(diff(range(datos$fecha))) / 365
  
  if (anios >= 2) {
    salto <- paste(max(1, round(anios / 10)), "years")
    formato <- "%Y"
    
  } else {
    salto <- "1 month"
    formato <-  "%b %Y"
  }
  
  ggplot(datos, aes(x = fecha, y = y)) +
    geom_line(linewidth = 0.7, color = "deeppink") +
    scale_x_date(date_breaks = salto, date_labels = formato) +
    labs(title = titulo, 
         x = "Fecha",
         y = paste0("Valor (", unidad, ")"),
         caption = caption) +
    theme_minimal() +
    theme(plot.caption = element_text(hjust = 0),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
}

# Correlograma

extraer_sucesion <- function(datos) {
  y <- if (is.data.frame(datos)) datos$y else datos
  y[!is.na(y)]
}

acf_manual <- function(serie, max_lag) {
  t <- length(serie)
  media <- mean(serie)
  
  # Denominador
  varianza <- sum((serie - media)^2)
  
  # Nuemerador
  acf_valores <- numeric(max_lag)
  
  for (h in 1:max_lag) {
    covarianza <- sum((serie[(h + 1):t] - media) * (serie[1:(t - h)] - media))
    acf_valores[h] <- covarianza / varianza
  }
  
  return(acf_valores)
  
}

barras_rezagos <- function(valores, limite, titulo, etiqueta_y) {
  m <- length(valores)
  paso <- ceiling(m / 24)
  tabla <- data.frame(rezago = seq_len(m), valor = valores)
  
  ggplot(tabla, aes(x = rezago, y = valor)) +
    geom_col(width = 0.25, fill = "deeppink") +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-limite, limite), linetype = "dashed", color = "#00CDCD") +
    scale_x_continuous(breaks = seq(paso, m, by = paso)) +
    labs(title = titulo, x = "Rezago", y = etiqueta_y) +
    theme_minimal()
  
}

correlograma <- function(datos, m = NULL) {
  y <- extraer_sucesion(datos)
  n <- length(y)
  
  if (is.null(m)) m <- min(floor(n / 4), 24)
  
  r <- acf_manual(y,m)
  p <- as.numeric(pacf(y, lag.max = m, plot = F)$acf)
  
  limite <- qnorm((1 + 0.95) / 2) / sqrt(n)
  
  grafico_acf <- barras_rezagos(r, limite, "ACF", "Autocorrelación")
  grafico_pacf <- barras_rezagos(p, limite, "PACF", "Autocorrelación parcial")
  
  grafico_acf / grafico_pacf
  
  
}

# Verificación

verificar_autocorrelacion <- function(datos, m = NULL) {
  y <- extraer_sucesion(datos)
  n <- length(y)
  if (is.null(m)) m <- min(floor(m / 4), 24)
  
  r_propia <- acf_manual(y, m)
  r_de_r <- as.numeric(acf(y, lag.max = m, plot = F)$acf)[-1]
  
  dif <- max(abs(r_propia - r_de_r))
  
  format(dif, scientific = TRUE)
  
}

# Grafico pronostico

grafico_pronostico <- function(datos, n_est, yhat, pronostico, titulo, ingenuo = NULL) {
  n <- nrow(datos)
  h <- n - n_est
  
  fuente <- attr(datos, "fuente")
  unidad <- attr(datos, "unidad")
  
  fechas_val <- datos$fecha[(n_est + 1):n]
  tabla <- rbind(
    data.frame(fecha = datos$fecha, valor = datos$y, serie = "Serie observada"),
    data.frame(fecha = datos$fecha[1:n_est], valor = yhat, serie = "Ajuste de un paso"),
    data.frame(fecha = fechas_val, valor = pronostico, serie = "Pronostico del metodo"))
  if (!is.null(ingenuo)) {
    tabla <- rbind(tabla, data.frame(fecha = fechas_val, valor = ingenuo,
                                     serie = "Pronostico ingenuo"))
  }
  tabla <- tabla[!is.na(tabla$valor), ]
  niveles <- c("Serie observada", "Ajuste de un paso", "Pronostico del metodo", "Pronostico ingenuo")
  tabla$serie <- factor(tabla$serie, levels = niveles)
  
  fuente_texto <- paste(strwrap(paste0("Fuente: ", fuente), width = 90), collapse = "\n")
  caption <- paste0(fuente_texto, "\nObservaciones: ", n, " (estimacion: ", n_est,
                    "; validacion: ultimas ", h, ", zona sombreada)")
  
  ggplot(tabla, aes(x = fecha, y = valor, color = serie, linetype = serie)) +
    annotate("rect", xmin = fechas_val[1], xmax = fechas_val[h],
             ymin = -Inf, ymax = Inf, fill = "gray70", alpha = 0.25) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = c("Serie observada" = "gray20", "Ajuste de un paso" = "steelblue4",
                                  "Pronostico del metodo" = "firebrick", "Pronostico ingenuo" = "darkorange3")) +
    scale_linetype_manual(values = c("Serie observada" = "solid", "Ajuste de un paso" = "solid",
                                     "Pronostico del metodo" = "solid", "Pronostico ingenuo" = "dashed")) +
    labs(title = titulo, x = "Fecha", y = paste0("Valor (", unidad, ")"),
         caption = caption, color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
}
