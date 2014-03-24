\name{summary.olmm}

\alias{summary.olmm}
\alias{print.summary.olmm}
\alias{print.olmm}

\title{Summarizing ordinal linear mixed model fits}

\description{
  Generates summary results of a fitted \code{\link{olmm}} object.
}

\usage{
\method{summary}{olmm}(object, silent = FALSE, ...)

\method{print}{olmm}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object, x}{a fitted \code{\link{olmm}} object.}
  \item{silent}{logical: should a warning be reported if the computation
    of the covariance matrix for the estimated coefficients failed.}
  \item{digits}{minimal number of _significant_ digits, see
    \code{\link{print.default}}.} 
  \item{...}{additional arguments.}
}

\value{
  Returns a list of class \code{"summary.olmm"} with components
  including: 
  \item{\code{methTitle}:}{character string. The title to be returned
    at \code{\link{show}} or \code{\link{print}} calls.}
  \item{\code{formula}:}{character string. The called formula.}
  \item{\code{data}:}{character string. The name of the data set
    used.} 
  \item{\code{subset}:}{character string. A description of the used
    subset of \code{data}.}
  \item{\code{AICtab}:}{a data frame of values of \code{\link{AIC}},
    \code{\link{BIC}}, log-likelihood and \code{\link{deviance}}.}
  \item{\code{FEmatEtaInv}:}{the matrix of estimates, standard errors,
    etc. for the predictor-invariant fixed-effects coefficients.}
  \item{\code{FEmatEtaVar}:}{the matrix of estimates, standard errors,
    etc. for the predictor-variable fixed-effects coefficients.}
  \item{\code{REmat}:}{the formatted Random-Effects matrix.}
  \item{\code{dims}:}{a named integer of dimensions.}  
}

\seealso{\code{\link{olmm}}, \code{\link{olmm-class}}}

\author{Reto Buergin}

\keyword{methods}