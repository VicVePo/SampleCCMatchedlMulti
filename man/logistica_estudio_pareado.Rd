% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{logistica_estudio_pareado}
\alias{logistica_estudio_pareado}
\title{Logística para Estudio Pareado}
\usage{
logistica_estudio_pareado(
  n_pares,
  m = 1,
  tasa_identificacion_pares,
  tasa_rechazo_pares,
  dias_laborables_mes
)
}
\arguments{
\item{n_pares}{Número de pares requeridos}

\item{m}{Número de controles por caso}

\item{tasa_identificacion_pares}{Tasa de identificación de pares por día}

\item{tasa_rechazo_pares}{Tasa de rechazo esperada para pares}

\item{dias_laborables_mes}{Número de días laborables por mes}
}
\value{
Una lista con los cálculos logísticos del estudio
}
\description{
Logística para Estudio Pareado
}
