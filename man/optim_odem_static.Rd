% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{optim_odem_static}
\alias{optim_odem_static}
\title{run static odem}
\usage{
optim_odem_static(
  p,
  input.values,
  nep = 1000,
  min = 100,
  sed = 3000,
  wind,
  khalf = 500,
  elev = NULL,
  verbose,
  startdate = NULL,
  enddate = NULL,
  field.values
)
}
\arguments{
\item{p}{estimated model parameters values}

\item{input.values}{input matrix of for instance thermocline depth}

\item{nep}{epilimnion net ecosystem production}

\item{min}{hypolimnion net ecosystem production}

\item{sed}{sediment oxygen demand}

\item{wind}{time series of wind velocity}

\item{khalf}{half-saturation coefficient}

\item{elev}{lake elevation}

\item{verbose}{verbose statement}

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
