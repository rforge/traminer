\name{olmm-control}

\alias{olmm_control}

\title{Control parameters for \code{\link{olmm}}.} 

\description{
  Various parameters that control aspects for \code{\link{olmm}}.  
}

\usage{
olmm_control(fit = c("ucminf", "nlminb", "optim"), 
             doFit = TRUE, numGrad = FALSE, 
             numHess = numGrad, nGHQ = 7L,
             start = NULL, restricted = NULL, verbose = FALSE, ...)
}

\arguments{
 \item{fit}{character string. The name of the function to be used for the 
    optimization.}
  \item{doFit}{logical scalar. When \code{FALSE} an unfitted \code{\link{olmm}}  
    object is returned. See details.}
  \item{numGrad}{logical scalar indicating whether the score function
    should be retrieved numerically.}
  \item{numHess}{logical scalar. Indicates whether the Hess matrix for
    the variance-covariance matrix should be estimated numerically,
    which is an approximation of the observed fisher information). Must
    be \code{TRUE} if \code{numGrad} is \code{TRUE}. See details.} 
  \item{nGHQ}{a positive integer specifying the number of quadrature
    points for the approximation of the marginal Likelihood by numerical
    integration.}
  \item{start}{a named numeric vector of initial values for the
    parameters. The parameter must be named in exactly in the way as
    they appear when the model is fitted.}
  \item{restricted}{a character vector of names of coefficients to be
    restricted to the initial values. The argument is ignored in case of
    adjacent category models.}
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the optimization of the parameter estimates.}
  \item{...}{further arguments to be passed to \code{fit}.} 
}

\details{
  Initial values may decrease the computation time and avoid
  divergence. The \code{start} argument accepts a vector with named
  elements according to the column names of the
  \code{\link{model.matrix}}. At the time being, initial values for
  adjacent-categories models must be transformed into the
  baseline-category model form.

  Notice that an additional argument \code{control}, e.g.,
  \code{control = list(trace = 1)}, can be passed access control
  parameters of the optimizers. For arguments, see \code{\link{ucminf}},
  \code{\link{nlminb}} or \code{\link{nlminb}}. 
}

\value{
  A list of class \code{olmm_control} containing
  the control parameters.
}

\author{Reto Buergin}

\seealso{\code{\link{olmm}}}

\examples{
olmm_control(doFit = FALSE)
}