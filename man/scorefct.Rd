% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scorefct.R
\name{scorefct}
\alias{scorefct}
\title{Internal functions for the MESS package}
\usage{
scorefct(object, beta = NULL, testidx = NULL, sas = FALSE)
}
\arguments{
\item{object}{input geepack object from a geeglm fit.}

\item{beta}{The estimated parameters. If set to \code{NULL} then
the parameter estimates are extracted from the model fit object
object.}

\item{testidx}{Indices of the beta parameters that should be tested
equal to zero}

\item{sas}{Logical. Should the SAS version of the score test be
computed. Defaults to \code{FALSE}.}
}
\description{
Internal functions for the MESS package
}
\author{
Claus Ekstrom \email{claus@ekstroem.dk}
}
\keyword{manip}
