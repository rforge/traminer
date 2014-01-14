\name{vcolmm2-forest}

\alias{fvcolmm}
\alias{logLik.fvcolmm}
\alias{neglogLik.fvcolmm}
\alias{plot.fvcolmm}
\alias{predict.fvcolmm}
\alias{print.fvcolmm}

\title{Bagging and random forests using the \code{\link{tvcolmm}} algorithm.}

\description{
  Implements bagging and random forests for the \code{\link{tvcolmm}}
  algorithm. 
}

\usage{

fvcolmm(..., folds = NULL,
               control = tvcolmm_control(alpha = 1.0,
                 bonferroni = FALSE, mtry = 5), add = TRUE)

\method{logLik}{fvcolmm}(object, ...)
\method{neglogLik}{fvcolmm}(object, ...)
\method{predict}{fvcolmm}(object, newdata = NULL,
                       type = c("prob", "class", "link",
                         "terms"), ranef = FALSE, 
                       verbose = FALSE, ...)
\method{print}{fvcolmm}(x, ...)
}

\arguments{
  \item{...}{for \code{\link{fvcolmm}}, this is either a
    specification of a model with the same arguments used for
    \code{\link{tvcolmm}}, i.e. \code{fvcolmm(formula = ..., data
      = ..., ...)} or a fitted \code{\link{fvcolmm}} object as
    first argument that should be added additional trees. For other
    functions, these are just supplementary arguments.}
  \item{object}{an object of class \code{\link{fvcolmm}}.}
  \item{x}{an object of class \code{fvcolmm}.}
  \item{folds}{a matrix of integers defining the folds for
    cross-validation. See \code{\link{cv.tvcolmm}}.}
  \item{control}{a \code{list} that specifies the maximum control
    parameters to evaluate. See \code{\link{tvcolmm_control}}.}
  \item{add}{logical scalar. Indicates if a existing forest should be
    extended with additional trees.}
  \item{newdata}{an optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{type}{character string indicating the type of prediction. See
    \code{\link{predict.tvcolmm}}.}  
  \item{ranef}{logical scalar. \code{FALSE} returns the marginal and
    \code{TRUE} the conditional prediction at 0.} 
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the validation.}
}

\details{
  Implements "bagging" (Breiman, 1996) and "random forests" (Breiman,
  2001) for \code{\link{tvcolmm}} objects.

  For bagging, use \code{control = list(alpha = 1L, mtry = Inf)}.

  Prediction for \code{forest.vcolmm} objects works in the same manner
  as \code{\link{predict.tvcolmm}}. The only current restriction is that
  the \code{ranef} argument must be logical.
}

\references{
  Leo Breiman (1996). Bagging Predictors. \emph{Machine Learning},
  123--140 
  Leo Breiman (2001). Random Forests. \emph{Machine Learning},
  \bold{45}(1), 5--32. 
}

\examples{

\dontrun{
## Extend Example 2 in '?tvcolmm'

## load the data
data(tvcolmm_2)

## fit the model

folds <- tvcolmm_folds(tvcolmm_2$id, type = "subsampling", K = 2)

forest <- fvcolmm(y ~ time + treat + (1|id) | z1 + z2 + z3 + z4,
                  data = tvcolmm_2,
                  control = tvcolmm_control(mtry = 2,
                    minsplit = 100, maxevalsplit = 5,
                    verbose = TRUE),
                  folds = folds)

}
}

\value{An object of class \code{fvcolmm}.}

\author{Reto Buergin}

\keyword{models}