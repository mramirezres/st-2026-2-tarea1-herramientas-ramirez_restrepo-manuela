library(tibble)

pasos <- c("1" = "1 year", 
           "4" = "3 months", 
           "12" = "1 month", 
           "52" = "1 week", 
           "7" = "1 day")

paso_de_frecuencia <- function(frecuencia) {
  clave <- as.character(frecuencia)
  if (!clave %in% names(pasos)) {
    stop("Frecuencia no soportada: ", frecuencia)
  }
  pasos[[clave]]
}

verificar_fechas <- function(fechas, frecuencia) {
  esperadas <- seq(fechas[1], 
                   by = paso_de_frecuencia(frecuencia),
                   length.out = length(fechas))
  stopifnot("Las fechas no son crecientes" = all(diff(fechas) > 0),
            "Las fechas no son equiespaciadas segun la frecuencia" = all(fechas == esperadas))
  
}

detectar_frecuencia <- function(fechas) {
  stopifnot("Se necesitan al menos 2 fechas" = length(fechas) >= 2)
  
  dias <- as.numeric(fechas[2] - fecha[1])
  if (dias == 1) return(7)
  if (dias == 7) return(52)
  if (dias %in% 28:31) return(12)
  if (dias %in% 365:366) return(1)
  
  stop("No se pudo deducir la frecuencia")
}

leer_serie <- function(x, fuente, unidad){
  if (is.ts(x)) {
    frecuencia <- frequency(x)
    
    st <- start(x)
    mes <- switch(as.character(frecuencia), 
                  "12" = st[2],
                  "4" = (st[2] - 1)*3+1,
                  "1" = 1
    )
    y <- as.numeric(x)
    fechas <- seq(as.Date(sprintf("%d-%02d-01", st[1], mes)),
                  by = paso_de_frecuencia(frecuencia),
                  length.out = length(y)
    )
  } else if (is.charecter(x) && length(x) == 1) {
    datos <- read.csv(x, stringsAsFactors = F)
    
    stopifnot("El csv debe de tener las columnas fecha y valor")
    
    fechas <- as.Date(datos$fecha)
    y <- datos$valor
    
    stopifnot("Hay fechas que no se pudieron leer",
              "La columna valor debe ser numerica")
    
    frecuencia <- detectar_frecuencia(fecha)
    
  }
  
  verificar_fechas(fechas, frecuencia)
  
  structure(
    tibble::tibble(t = seq_along(y), fecha = fechas, y = y),
    frecuencia = frecuencia,
    fuente = fuente,
    unidad = unidad
  )
  
}