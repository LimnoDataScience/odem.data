% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{odem_static}
\alias{odem_static}
\title{run static odem}
\usage{
odem_static(
  input.values,
  sed,
  nep,
  min,
  wind = NULL,
  khalf = 500,
  elev = 450,
  startdate = NULL,
  enddate = NULL,
  field.values
)
}
\arguments{
\item{input.values}{input matrix of for instance thermocline depth}

\item{sed}{sediment oxygen demand}

\item{nep}{epilimnion net ecosystem production}

\item{min}{hypolimnion net ecosystem production}

\item{wind}{time series of wind velocity}

\item{khalf}{half-saturation coefficient}

\item{elev}{lake elevation}

\item{startdate}{start data}

\item{enddate}{end date}

\item{field.values}{observed oxygen data}
}
\value{
list of output data, fit metric, and plot
}
\description{
run static odem
}
