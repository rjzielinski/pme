% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpme.R
\name{new_lpme}
\alias{new_lpme}
\title{Create New LPME Object}
\usage{
new_lpme(
  embedding_map,
  msd,
  coefficients,
  times,
  initial_parameterization,
  optimal_parameterization,
  d,
  D,
  n_knots,
  lambda,
  gamma,
  coefficient_list,
  parameterization_list,
  smoothing_method
)
}
\arguments{
\item{embedding_map}{The optimal embedding function.}

\item{msd}{Vector of Mean Squared Distance Values.}

\item{coefficients}{sol_coef.}

\item{times}{A value.}

\item{initial_parameterization}{A value.}

\item{optimal_parameterization}{A value.}

\item{d}{A value.}

\item{D}{A value.}

\item{n_knots}{A value.}

\item{lambda}{A value.}

\item{gamma}{A value.}

\item{coefficient_list}{A value.}

\item{parameterization_list}{A value.}

\item{smoothing_method}{A value.}
}
\value{
An object with class "lpme".
}
\description{
Create New LPME Object
}
