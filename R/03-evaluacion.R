ljung_box <- function(r, T, m, p) {
  h <- 1:m
  Q <- T * (T + 2) * sum(r[h]^2 / (T - h))
  gl <- m - p
  
  list(estadistico = Q,
       gl = gl,
       valor_critico = qchisq(0.95, df = gl),
       valor_p = pchisq(Q, df = gl, lower.tail = F)
  )
  
}

jarque_bera <- function(e) {
  e <- e[!is.na(e)]
  N <- length(e)
  
  media <- mean(e)
  m2 <- mean((e - media)^2)
  m3 <- mean((e - media)^3)
  m4 <- mean((e - media)^4)
  
  A <- m3 / sqrt(m2)^3
  K <- m4 / m2^2
  JB <- N / 6 * (A^2 + (K - 3)^2 / 4)
  
  list(estadistico = JB,
       gl = 2,
       valor_critico = qchisq(0.95, df = 2),
       valor_p = pchisq(JB, df = 2, lower.tail = F)
  )
  
}

durbin_watson <- function(e) {
  e <- e[!is.na(e)]
  N <- length(e)
  
  d <- sum((e[2:N] - e[1:(N - 1)])^2) / sum(e^2)
  
  list(estadistico = d,
       gl = NA,
       valor_critico = NA,
       valor_p = NA
  )
  
  
}

# Verificación

verificar_ljung_box <- function(serie, m, p = 0) {
  serie <- serie[!is.na(serie)]
  n <- length(serie)
  
  propia <- ljung_box(acf_manual(serie, m), n, m, p)
  de_r <- Box.test(serie, lag = m, type = "Ljung-Box", fitdf = p)
  
  dif_Q <- abs(propia$estadistico - unname(de_r$statistic))
  dif_p <- abs(propia$valor_p - de_r$p.value)
  
  cat("Ljung-Box (m =", m, ", p =", p, ")\n",
      "  Q propio:", format(propia$estadistico, digits = 10),
      "| Q de R:", format(unname(de_r$statistic), digits = 10), "\n")
  
  invisible(max(dif_Q, dif_p))
  
}


medidas <- function(y, yhat, y_est, periodo = 1) {
  stopifnot("periodo debe ser un entero positivo" = is.numeric(periodo) && length(periodo) == 1 && periodo == round(periodo) && periodo >= 1)
  
  usar <- !is.na(yhat)
  y <- y[usar]
  yhat <- yhat[usar]
  error <- y - yhat
  
  escala <- mean(abs(diff(y_est, lag = periodo)))
  mape <- if (any(y == 0)) NA_real_ else 100 * mean(abs(error / y))
  
  list(MSE  = mean(error^2),
       MAD  = mean(abs(error)),
       MAPE = mape,
       MASE = mean(abs(error)) / escala,
       n    = length(error))
  
}

# Validar errores

validar_errores <- function(e, p = 0, m = NULL, fechas = NULL) {
  
  usar <- !is.na(e)
  errores <- e[usar]
  n <- length(e)
  
  if (is.null(m)) m <- min(floor(n / 4), 24)
  
  # Prueba t de media cero
  desviacion <- sd(errores)
  t_est <- mean(errores) / (desviacion / sqrt(n))
  prueba_t <- list(estadistico   = t_est,
                   gl            = n - 1,
                   valor_critico = qt(0.975, df = n - 1),
                   valor_p       = 2 * pt(abs(t_est), df = n - 1, lower.tail = FALSE))
  
  # Ljung-Box, Jarque-Bera y Durbin-Watson
  lb <- ljung_box(acf_manual(errores, m), n, m, p)
  jb <- jarque_bera(errores)
  dw <- durbin_watson(errores)
  
  # Grafico de los errores
  tiempo <- if (is.null(fechas)) seq_along(e) else fechas
  tabla <- data.frame(tiempo = tiempo[usar], error = errores)
  grafico_errores <- ggplot(tabla, aes(x = tiempo, y = error)) +
    geom_hline(yintercept = 0, color = "gray40") +
    geom_line(color = "deeppink") +
    geom_point(size = 1, color = "deeppink4") +
    labs(title = "Errores de un paso",
         x = if (is.null(fechas)) "Periodo" else "Fecha",
         y = "Error (y - yhat)") +
    theme_minimal()
  
  resumen <- data.frame(
    prueba        = c("t de media cero", "Ljung-Box", "Jarque-Bera", "Durbin-Watson"),
    estadistico   = c(prueba_t$estadistico, lb$estadistico, jb$estadistico, dw$estadistico),
    gl            = c(prueba_t$gl, lb$gl, jb$gl, dw$gl),
    valor_critico = c(prueba_t$valor_critico, lb$valor_critico, jb$valor_critico, dw$valor_critico),
    valor_p       = c(prueba_t$valor_p, lb$valor_p, jb$valor_p, dw$valor_p))
  
  list(n = n,
       grafico_errores = grafico_errores,
       correlograma = correlograma(errores, m),
       prueba_t = prueba_t,
       ljung_box = lb,
       jarque_bera = jb,
       durbin_watson = dw,
       resumen = resumen)
}