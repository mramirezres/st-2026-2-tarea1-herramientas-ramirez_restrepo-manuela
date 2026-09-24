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