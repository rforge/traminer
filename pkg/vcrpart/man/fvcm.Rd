\name{fvcm}

\alias{fvcm}
\alias{fvcm_control}
\alias{fvcolmm}
\alias{fvcglm}

\title{Bagging and random forests based on \code{\link{tvcm}}.}

\description{
  Implements the bagging (Breiman, 1996) and random forest (Breiman, 2001) ensemble 
  algorithms for \code{\link{tvcm}}.  
}

\usage{

fvcm(..., folds, control = fvcm_control())

fvcolmm(..., family = cumulative(), folds, control = fvcm_control())

fvcglm(..., folds, control = fvcm_control())

fvcm_control(alpha = 1.0, maxwidth = 10L, minbucket = 50L, mtry = 5L, ...)
}

\arguments{
  \item{...}{for \code{\link{fvcm}}, \code{\link{fvcolmm}} and \code{\link{fvcglm}}, 
    \code{...} are arguments to be passed to \code{\link{tvcm}}, see \code{\link{tvcm}}. 
    This includes necessarily the arguments \code{formula}, \code{data} and
    \code{family}. For \code{\link{fvcm_control}}, \code{...} are
    further arguments to be passed to \code{\link{tvcm_control}}.}  
  \item{folds}{a matrix of integers defining the rows of observations
    in \code{data} on which single trees are fitted, e.g., a matrix produced by
    \code{\link{cvfolds}}. See \code{\link{cvfolds}}.}  
  \item{control}{control parameters as produced by
    \code{\link{fvcm_control}}.} 
  \item{family}{the model family, e.g., an object of class 
    \code{\link{family.olmm}} or \code{\link{family}}.} 
  \item{alpha}{numeric significance threshold between 0 and 1. A node is 
    splitted when the smallest \eqn{p} value for any coefficient
    coefficient test in the current step falls below \code{alpha}.}
  \item{maxwidth}{integer. The maximum width of the tree.}
  \item{minbucket}{integer. The minimum sum of weights in each node.}
  \item{mtry}{integer. The number of input variables randomly sampled as
    candidates at each node. Is to be set to \code{Inf} for bagging. See details.} 
}

\details{
  Implements the \emph{bagging} (Breiman, 1996) and \emph{random
  forests} (Breiman, 2001) ensemble algorithms for the
  \code{\link{tvcm}} algorithm. The method consist in growing multiple
  trees using \code{\link{tvcm}} and aggregating the fitted coefficient
  functions. For bagging, simply set \code{mtry = Inf} in \code{\link{fvcm_control}}. 

  \code{\link{fvcolmm}} and \code{\link{fvcglm}} are convenience
  functions for whether a \code{\link{olmm}} or a \code{\link{glm}}
  model is fitted. 

  \code{\link{fvcm_control}} is essentially \code{\link{tvcm_control}}
  with other defaults for the main tuning parameters, which are
  \code{maxwidth}, \code{minbucket} and \code{mtry}. \code{alpha} is set
  to one so that the stopping parameter for \code{type = "mob"} is
  skipped. The parameters \code{maxwidth} and \code{minbucket} are
  discussed in the details section of \code{\link{tvcm}}. The
  \code{mtry} argument, i.e., the number of randomly chosen variables
  for splitting in each step, may be specified as the number of
  partitioning  variables divided by 3, see Hastie et al. (2001),
  sec. 15. 
}

\value{An object of class \code{fvcm}.}

\author{Reto Buergin}

\references{
  Leo Breiman (1996). Bagging Predictors. \emph{Machine Learning},
  123--140
  
  Leo Breiman (2001). Random Forests. \emph{Machine Learning},
  \bold{45}(1), 5--32.

  T. Hastie, R. Tibshirani, J. Friedman (2001), The elements of
  statistical learning, Springer.
}

\seealso{\code{\link{fvcm-methods}}, \code{\link{tvcm}},
  \code{\link{glm}}, \code{\link{olmm}}}

\examples{
## ------------------------------------------------------------------- #
## Dummy example 1:
##
## Fitting a random forest 'tvcm' on artificially generated data.
## ------------------------------------------------------------------- #

## simulate data
data(vcrpart_2)

## set parameters
control <- fvcm_control(mtry = 1, minbucket = 25, verbose = TRUE)
folds <- cvfolds(factor(1:100), "subsampling", 3, 0.75)

## fit the forest
model <- fvcm(y ~ vc(z1, z2, by = x1) + x2, data = vcrpart_2,
              family = gaussian(), 
              folds = folds, control = control)

model
}

\keyword{models}