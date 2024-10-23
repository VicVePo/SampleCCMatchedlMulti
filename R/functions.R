#' Cálculo de Tamaño Muestral para Regresión Logística con Datos Pareados
#' 
#' @param alfa Nivel de significación
#' @param potencia Potencia deseada (opcional si se proporciona n)
#' @param n Tamaño muestral (opcional si se proporciona potencia)
#' @param OR Odds Ratio a detectar
#' @param p_discordante Proporción de pares discordantes
#' @param m Número de controles por caso (por defecto 1 para pareamiento 1:1)
#' @param metodo Método de cálculo: "rho" o "ee" (error estándar)
#' @param valores_rho Vector de valores rho para calcular tamaños muestrales
#' @param EE Error estándar del coeficiente (no del OR)
#' @param IC_superior Límite superior del intervalo de confianza del OR
#' @param IC_inferior Límite inferior del intervalo de confianza del OR
#' @param n_previo Tamaño de muestra del estudio previo (requerido para método "ee")
#' @return Un data frame con resultados según el método elegido
#' @export
SampleCCMatchedlMulti <- function(alfa, potencia = NULL, n = NULL,
                                     OR = NULL,
                                     p_discordante,
                                     m = 1,
                                     metodo = "rho",
                                     valores_rho = seq(0, 0.9, by = 0.1),
                                     EE = NULL,
                                     IC_superior = NULL, 
                                     IC_inferior = NULL,
                                     n_previo = NULL) {
  
  if (metodo == "ee") {
    if (is.null(EE) && (is.null(IC_superior) || is.null(IC_inferior))) {
      stop("Para método 'ee' se requiere EE o ambos límites del IC")
    }
    if (is.null(n_previo)) {
      stop("Para método 'ee' se requiere n_previo")
    }
  }
  
  # Si se proporcionaron los IC, calcular el EE
  if (!is.null(IC_superior) && !is.null(IC_inferior)) {
    EE <- (log(IC_superior) - log(OR)) / qnorm(0.975)
    cat("Error Estándar del coeficiente calculado:", EE, "\n")
  }
  
  z_alfa <- qnorm(1 - alfa/2)
  
  if (metodo == "rho") {
    if (!is.null(potencia)) {
      z_beta <- qnorm(potencia)
      
      calcular_n <- function(rho) {
        numerador <- (z_alfa * sqrt((m + 1) / m) + z_beta)^2
        denominador <- p_discordante * (log(OR))^2 * (1 - rho^2)
        
        n_pares <- ceiling(numerador / denominador)
        return(n_pares)
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        n_pares = sapply(valores_rho, calcular_n)
      )
      
      resultados$tamaño_total <- resultados$n_pares * (m + 1)
      
    } else {
      # Cálculo de potencia
      calcular_potencia <- function(rho) {
        z_beta <- sqrt(n * p_discordante * (log(OR))^2 * (1 - rho^2)) / 
          sqrt((m + 1) / m) - z_alfa
        potencia <- pnorm(z_beta)
        return(potencia)
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        potencia = sapply(valores_rho, calcular_potencia)
      )
    }
  } else {
    # Método basado en error estándar
    if (!is.null(potencia)) {
      z_gamma <- qnorm(potencia)
      n <- ceiling((z_alfa + z_gamma)^2 * n_previo * EE^2 / (log(OR))^2)
      resultados <- data.frame(
        n_pares = n
      )
      resultados$tamaño_total <- n * (m + 1)
    } else {
      z_gamma <- sqrt(n * (log(OR))^2 / (n_previo * EE^2)) - z_alfa
      potencia <- pnorm(z_gamma)
      resultados <- data.frame(
        potencia = potencia,
        n_pares = n,
        tamaño_total = n * (m + 1)
      )
    }
  }
  
  return(resultados)
}

#' Logística para Estudio Pareado
#'
#' @param n_pares Número de pares requeridos
#' @param m Número de controles por caso
#' @param tasa_identificacion_pares Tasa de identificación de pares por día
#' @param tasa_rechazo_pares Tasa de rechazo esperada para pares
#' @param dias_laborables_mes Número de días laborables por mes
#' @return Una lista con los cálculos logísticos del estudio
#' @export
logistica_estudio_pareado <- function(n_pares,
                                     m = 1,
                                     tasa_identificacion_pares,
                                     tasa_rechazo_pares,
                                     dias_laborables_mes) {
  
  pares_a_contactar <- n_pares / (1 - tasa_rechazo_pares)
  
  dias_total <- ceiling(pares_a_contactar / tasa_identificacion_pares)
  meses_total <- dias_total / dias_laborables_mes
  
  resultados <- list(
    pares_requeridos = n_pares,
    participantes_total = n_pares * (m + 1),
    pares_a_contactar = ceiling(pares_a_contactar),
    dias_reclutamiento = dias_total,
    meses_reclutamiento = round(meses_total, 2)
  )
  
  # Imprimir resumen
  cat("
Resumen de logística del estudio pareado:
")
  cat("----------------------------------------
")
  cat("Pares requeridos:", n_pares, "
")
  cat("Total de participantes:", n_pares * (m + 1),
      "(", m, ":", 1, "pareamiento)
")
  cat("Pares a contactar:", ceiling(pares_a_contactar),
      "(", tasa_rechazo_pares*100, "% de rechazo)
")
  cat("Días necesarios:", dias_total,
      "(", tasa_identificacion_pares, "pares por día)
")
  cat("Meses necesarios:", round(meses_total, 2),
      "(", dias_laborables_mes, "días laborables por mes)
")
  cat("----------------------------------------
")
  
  return(invisible(resultados))
}
