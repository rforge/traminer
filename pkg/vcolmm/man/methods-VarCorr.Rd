\docType{methods}
\name{VarCorr}

\alias{VarCorr}

\title{
  Extract variance and correlation components
}

\description{
  This function calculates the estimated variances, standard
  deviations, and correlations between the random-effects terms in a
  mixed-effects model. The within-group error variance and standard
  deviation are also calculated.  
}

\usage{
  VarCorr(x, sigma, rdig)
}

\arguments{
  \item{x}{any fitted model object from which variance components
    of random-effects estimates can be extracted.}
  \item{sigma}{an optional numeric value used as a multiplier for the
    standard deviations. Default is \code{1}.}
  \item{rdig}{an optional integer value specifying the number of digits
    used to represent correlation estimates. Default is \code{3}.}
}

\references{
  Pinheiro, J.C., and Bates, D.M. (2000) Mixed-Effects Models in S and
  S-PLUS, Springer, esp. pp. 100, 461.   
}

\seealso{
  \code{\link{olmm}}, \code{\link{olmm-class}}
}

\keyword{methods}