# Validación

validar_serie <- function(y) {
  stopifnot("y debe ser un vector numerico" = is.numeric(y),
            "y no puede tener NA" = !anyNA(y),
            "y debe tener valores finitos" = all(is.finite(y)),
            "y debe tener al menos 2 observaciones" = length(y) >= 2)
  
}

es_constante <- function(valor) {
  is.numeric(valor) && length(valor) == 1 && !is.na(valor) && valor > 0 && valor < 1
}

validar_ventana <- function(k, n, minimo) {
  stopifnot("k debe ser un numero entero" = is.numeric(k) && length(k) == 1 && k %% 1 == 0,
            "k es menor que el minimo permitido" = k >= minimo)
}

# Media simple

ajustar_media <- function(y) {
  validar_serie(y)
  n <- length(y)
  
  yhat <- rep(NA_real_, n)
  media <- y[1]
  yhat[2] <- media
  for (t in 2:n) {
    media <- media + (y[t] - media) / t
    if(t < n) {
      yhat[t + 1] <- media
    }
    
  }
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" =
                is.numeric(h) && length(h) == 1 && h == round(h) && h >= 1)
    rep(media, h)
  }
  
  list(yhat = yhat,
       pronosticar = pronosticar,
       parametros = list(media = media, n = n))
  
}

# Media Movil

ajustar_mm <- function(y, k) {
  validar_serie(y)
  n <- length(y)
  validar_ventana(k , n, minimo = 2)
  stopifnot("La ventana k no puede ser mayor que la longitud de la serie" = k <= n)
  
  yhat <- rep(NA_real_, n)
  suma <- 0
  mm <- NA_real_
  for (t in 1:n) {
    suma <- suma + y[t]
    if(t > k) {
      suma <- suma - y[t - k]
    }
    
    if (t >= k) {
      mm <- suma / k
      if(t < n) {
        yhat[t + 1] <- mm
      }
    }
  }
  
  ultima_mm <- mm
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" =
                is.numeric(h) && length(h) == 1 && h == round(h) && h >= 1)
    rep(ultima_mm, h)
    
  }
  
  list(yhat = yhat,
       pronosticar = pronosticar,
       parametros = list(k = k, ultima_mm = ultima_mm))
  
  
}

# SES

ajustar_ses <- function(y, alpha) {
  validar_serie(y)
  stopifnot("alpha debe ser un numero en el intervalo (0, 1)" = es_constante(alpha))
  n <- length(y)
  
  yhat <- rep(NA_real_, n)
  yhat[2] <- y[1]
  
  siguiente <- yhat[2]
  for (t in 2:n) {
    error <- y[t] - yhat[t]
    siguiente <- yhat[t] + alpha * error
    if (t < n) {
      yhat[t + 1] <- siguiente
      
    }
  }
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" =
                is.numeric(h) && length(h) == 1 && h == round(h) && h >= 1)
    rep(siguiente, h)
    
  }
  
  list(yhat = yhat,
       pronosticar = pronosticar,
       parametro = list(alpha = alpha, nivel_final = siguiente))
  
}

# Doble media movil

ajustar_dmm <- function(y, k) {
  validar_serie(y)
  n <- length(y)
  validar_ventana(k, n, minimo = 2)
  stopifnot("Se necesitan al menos 2k - 1 observaciones" = n >= 2 * k - 1)
  
  mm    <- rep(NA_real_, n)
  E     <- rep(NA_real_, n)
  beta1 <- rep(NA_real_, n)
  yhat  <- rep(NA_real_, n)
  
  suma_y  <- 0                          
  suma_mm <- 0                         
  for (t in 1:n) {
    suma_y <- suma_y + y[t]
    if (t > k) {
      suma_y <- suma_y - y[t - k]
    }
    if (t >= k) {
      mm[t] <- suma_y / k
      suma_mm <- suma_mm + mm[t]
      if (t >= 2 * k) {
        suma_mm <- suma_mm - mm[t - k]
      }
    }
    if (t >= 2 * k - 1) {
      dmm <- suma_mm / k
      E[t] <- 2 * mm[t] - dmm
      beta1[t] <- 2 / (k - 1) * (mm[t] - dmm)
      if (t < n) {
        yhat[t + 1] <- E[t] + beta1[t]
      }
    }
  }
  
  nivel_final <- E[n]
  pendiente_final <- beta1[n]
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" =
                is.numeric(h) && length(h) == 1 && h == round(h) && h >= 1)
    nivel_final + pendiente_final * (1:h)
  }
  
  list(yhat = yhat,
       pronosticar = pronosticar,
       parametros = list(k = k, E = E, beta1 = beta1,
                         nivel_final = nivel_final,
                         pendiente_final = pendiente_final))
}

# Holt lineal

ajustar_holt <- function(y, alpha, beta) {
  validar_serie(y)
  stopifnot("alpha debe ser un numero en el intervalo (0, 1)" = es_constante(alpha),
            "beta debe ser un numero en el intervalo (0, 1)" = es_constante(beta))
  
  n <- length(y)
  
  nivel <- numeric(n)
  tendencia <- numeric(n)
  nivel[1] <- y[1]
  tendencia[1] <- 0
  
  yhat <- rep(NA_real_, n)
  yhat[2] <- nivel[1] + tendencia[1]
  
  for (t in 2:n) {
    nivel[t] <- alpha * y[t] + (1 - alpha) * yhat[t]
    tendencia[t] <- beta * (nivel[t] - nivel[t - 1]) + (1 - beta) * tendencia[t - 1]
    
    if (t < n) {
      yhat[t + 1] <- nivel[t] + tendencia[t]
    }
    
  }
  
  nivel_final <- nivel[n]
  tendencia_final <- tendencia[n]
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" =
                is.numeric(h) && length(h) == 1 && h == round(h) && h >= 1)
    nivel_final + tendencia_final * (1:h)
    
  }
  
  list(yhat = yhat,
       pronosticar = pronosticar,
       parametros = list(alpha = alpha, beta = beta,
                         nivel = nivel, tendencia = tendencia,
                         nivel_final = nivel_final,
                         tendencia_final = tendencia_final))
  
}

# Tendencia

ajustar_tendencia <- function(y, tipo, corregir_sesgo = F) {
  validar_serie(y)
  n <- length(y)
  p <- if (tipo == "cuadratica") 3 else 2
  
  if (tipo == "exponencial") {
    stopifnot("La tendencia exponencial exige y > 0" = all(y > 0))
  }
  
  matriz_diseno <- function(tiempo) {
    if (tipo == "cuadratica") return(cbind(1, tiempo, tiempo^2))
    cbind(1, tiempo)
  }
  
  tiempo <- 1:n
  X <- matriz_diseno(tiempo)
  z <- if (tipo == "exponencial") log(y) else y
  
  coef <- as.numeric(solve(crossprod(X), crossprod(X, z)))
  residuos <- as.numeric(z - X %*% coef)
  gl <- n - p
  sigma2 <- sum(residuos^2) / gl
  r2 <- 1 - sum(residuos^2) / sum((z - mean(z))^2)
  durbin_watson_residuos <- sum(diff(residuos)^2) / sum(residuos^2)
  
  inversa <- solve(crossprod(X))
  ee_ordinario <- sqrt(diag(sigma2 * inversa))
  
  rezagos <- floor(4 * (n / 100)^(2 / 9))
  u <- X * residuos
  S <- crossprod(u)
  
  for (l in seq_len(rezagos)) {
    G <- crossprod(u[(l + 1):n, , drop = FALSE], u[1:(n - l), , drop = FALSE])
    S <- S + (1 - l / (rezagos + 1)) * (G + t(G))
  }
  
  ee_robusto <- sqrt(diag(inversa %*% S %*% inversa))
  
  nombres <- switch (tipo,
                     lineal = c("beta0", "beta1"),
                     cuadratica = c("beta0", "beta1", "beta2"),
                     exponencial = c("a", "teta")
  )
  
  t_ordinario <- coef / ee_ordinario
  t_robusto <- coef / ee_robusto
  
  
  
  tabla <- data.frame(estimacion = coef,
                      ee_ordinario = ee_ordinario,
                      ee_robusto = ee_robusto,
                      t_ordinario = t_ordinario,
                      p_ordinario = 2 * pt(abs(t_ordinario), df = gl, lower.tail = FALSE),
                      t_robusto = t_robusto,
                      p_robusto = 2 * pt(abs(t_robusto), df = gl, lower.tail = FALSE),
                      row.names = nombres)
  
  factor <- if (corregir_sesgo) exp(sigma2 / 2) else 1
  evaluar <- function(tiempo){
    valores <- as.numeric(matriz_diseno(tiempo) %*% coef)
    if (tipo == "exponencial") valores <- exp(valores) * factor
    valores
    
  }
  
  yhat <- evaluar(tiempo)
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" =
                is.numeric(h) && length(h) == 1 && h == round(h) && h >= 1)
    evaluar((n + 1):(n + h))
  }
  
  if (tipo == "lineal") {
    coeficientes <- c(beta0 = coef[1], beta1 = coef[2])
    
  }
  
  else if (tipo == "cuadratica") {
    coeficientes <- c(beta0 = coef[1], beta1 = coef[2], beta2 = coef[3])
  }
  
  else {
    coeficientes <- c(a = coef[1], teta = coef[2],
                      beta0 = exp(coef[1]), beta1 = exp(coef[2]))
  }
  
  list(yhat = yhat,
       pronosticar = pronosticar,
       parametros = list(tipo = tipo, coeficientes = coeficientes,
                         tabla = tabla, r2 = r2, sigma2 = sigma2,
                         durbin_watson = durbin_watson_residuos,
                         rezagos_hac = rezagos,
                         corregir_sesgo = corregir_sesgo, n = n))
}

# Optimización

mse_un_paso <- function(y, yhat) {
  error <- y - yhat
  usar <- !is.na(error)
  if (!any(usar)) return(NA_real_)
  mean(error[usar]^2)
  
}

optimizar <- function(y, metodo, rejilla = NULL){
  validar_serie(y)
  n <- length(y)
  
  if (is.null(rejilla)) {
    if (metodo %in% c("mm", "dmm")) rejilla <- 2:12
    if (metodo == "ses") rejilla <- seq(0.02, 0.98, by = 0.02)
    if (metodo == "holt") {
      valores <- seq(0.05, 0.95, by = 0.05)
      rejilla <- expand.grid(alpha = valores, beta = valores)
    }
  }
  
  if (metodo == "holt") {
    tabla <- data.frame(alpha = rejilla$alpha, beta = rejilla$beta)
  } else {
    tabla <- data.frame(rejilla)
    names(tabla) <- if (metodo == "ses") "alpha" else "k"
    
  }
  
  tabla$mse <- NA_real_
  
  for (i in seq_len(nrow(tabla))) {
    ajuste <- NULL
    
    if (metodo == "mm" && tabla$k[i] <= n - 1){
      ajuste <- ajustar_mm(y, tabla$k[i])
    } else if (metodo == "dmm" && 2 * tabla$k[i] <= n) {
      ajuste <- ajustar_dmm(y, tabla$k[i])
    } else if (metodo == "ses") {
      ajuste <- ajustar_ses(y, tabla$alpha[i])
    } else if (metodo == "holt") {
      ajuste <- ajustar_holt(y, tabla$alpha[i], tabla$beta[i])
    }
    
    if (!is.null(ajuste)) {
      tabla$mse[i] <- mse_un_paso(y, ajuste$yhat)
    }
    
  }
  
  mejor <- which.min(tabla$mse)
  optimo <- tabla[mejor, , drop = F]
  rownames(optimo) <- NULL
  
  parametros <- setdiff(names(tabla), "mse")
  en_borde <- FALSE
  for (col in parametros) {
    valores <- tabla[[col]]
    if (length(unique(valores)) > 1) {
      if (valores[mejor] == min(valores) || valores[mejor] == max(valores)) {
        en_borde <- TRUE
      }
    }
  }
  
  list(metodo = metodo, rejilla = tabla, optimo = optimo, en_borde = en_borde)
  
}