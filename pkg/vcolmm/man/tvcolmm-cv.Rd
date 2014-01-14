\name{tvcolmm-cv}			

\alias{cv.tvcolmm}
\alias{plot.cv.tvcolmm}
\alias{tvcolmm_folds}
\alias{stabpath.tvcolmm}
\alias{plot.stabpath.tvcolmm}

\title{Cross-validation for \code{\link{tvcolmm}} objects.}

\description{
  Estimates of the loss or stability paths of \code{\link{tvcolmm}}
  objects for varying alpha's via cross-validation. 
}

\usage{

tvcolmm_folds(x, type = c("subsampling", "kfold", "bootstrap"),
              K = ifelse(type == "kfold", 10L, 50L),
              prob = 0.5)

cv.tvcolmm(object, folds = tvcolmm_folds(object),
           alpha.max = 0.5, verbose = FALSE, ...)

\method{plot}{cv.tvcolmm}(x, ...)

stabpath.tvcolmm(object, q = 2L, ...)

\method{plot}{stabpath.tvcolmm}(x, ...)

}

\arguments{
  \item{object}{an object of class \code{\link{tvcolmm}}.}
  \item{x}{for \code{\link{tvcolmm_folds}}, this must be either a
    \code{factor} that represents the subject vector or an object of
    class \code{\link{tvcolmm}}. For the \code{plot} methods, this
    must be either an object of class \code{cv.tvcolmm} or
    \code{stabpath.tvcolmm}.}
  \item{alpha.max}{scalar numeric. Maximal alpha to be evaluated.}
  \item{folds}{a matrix of integers. User specified folds where columns
    extract subsets of the data. See \code{\link{tvcolmm_folds}}.}
  \item{K}{integer. The number of folds to evaluate.} 
  \item{prob}{numeric. The probability for the \code{"subsampling"} 
    cross-validation scheme.}
  \item{q}{average number of selected variables.}
  \item{type}{character string. The type of sampling scheme to be used.}  
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the validation.}
  \item{...}{other arguments to be passed to \code{cv.vcolmm} or
    \code{plot.default}.} 
}

\details{\code{cv.tvcolmm} estimates the negative log-likelihood
  prediction error via repeated fitting on subsamples.

  The sampling of subsets may be user-defined (by assigning a matrix to 
  \code{folds}). By default, the \code{\link{tvcolmm_folds}} function 
  is called selecting randomly units from the grouping factor
  (i.e. subjects) defined in the random effects part of the model
  formula of the \code{\link{tvcolmm}} object.

  The function \code{stabpath.tvcolmm} implements the "stability paths" of
  Meinshausen and Buehlmann (2010). 
}

\value{An object of class \code{cv.tvcolmm} or \code{stabpath.tvcolmm}.}

\references{
  Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) 
  \emph{Classification and Regression Trees}. Wadsworth.

  Meinshausen, N. and Buehlmann, P. (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}(4).
}

\seealso{
  \code{\link{tvcolmm}}
}

\examples{

\dontrun{
## Validate Example 2 in '?tvcolmm'

## load the data
data(tvcolmm_2)

## fit the model
mod <- tvcolmm(y ~ time + treat + (1|id) | z1 + z2 + z3, data = tvcolmm_2,
               control = tvcolmm_control(
                alpha = 0.05, bonferroni = FALSE),
                linear = "time")

## cross-validation
folds <- tvcolmm_folds(mod, type = "kfold", K = 10)
cv <- cv.tvcolmm(mod, folds, control = list(alpha = 0.5))

## Overfitting is hardly a problem here.
plot(cv, log = "x")

## Stability paths
sl <- stabpath.tvcolmm(mod, q = 1, folds = folds,
                       control = list(alpha = 0.1))
plot(sl, log = "x")
}

}

\author{Reto Buergin}

\keyword{validation}
