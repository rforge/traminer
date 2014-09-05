\name{fvcm}

\alias{fvcm}
\alias{fvcm_control}
\alias{fvcolmm}
\alias{fvcglm}

\title{Bagging and random forests based on \code{\link{tvcm}}}

\description{
  Bagging (Breiman, 1996) and random forest (Breiman, 2001) ensemble
  algorithms for \code{\link{tvcm}}.
}

\usage{

fvcm(..., folds = folds_control("subsampling", 5), control = fvcm_control())

fvcolmm(..., family = cumulative(), folds = folds_control("subsampling", 5),
        control = fvcm_control())

fvcglm(..., family, folds = folds_control("subsampling", 5), control = fvcm_control())

fvcm_control(alpha = 1.0, maxstep = 10L, minsize = NULL, ptry = 1,
             ntry = 1, vtry = 5L, ...)
}

\arguments{
  \item{\ldots}{for \code{\link{fvcm}}, \code{\link{fvcolmm}} and
    \code{\link{fvcglm}}, \code{\ldots} are arguments to be passed to
    \code{\link{tvcm}}, see \code{\link{tvcm}}. This includes at least
    the arguments \code{formula}, \code{data} and \code{family} (see
    examples). For \code{\link{fvcm_control}}, \code{\ldots} are further
    arguments to be passed to \code{\link{tvcm_control}}.}   
  \item{folds}{a list of parameters for folds as created by
    \code{\link{folds_control}}.}   
  \item{control}{control parameters as produced by
    \code{\link{fvcm_control}}.} 
  \item{family}{the model family, e.g., \code{\link{cumulative}} for
    fitting a cumulative link model.} 
  \item{alpha}{numeric significance threshold between 0 and 1. A node is 
    splitted when the smallest \eqn{p} value for any coefficient
    coefficient test in the current step falls below \code{alpha}.}
  \item{maxstep}{integer. The maximum number of steps for when growing
    individual trees.}
  \item{minsize}{integer. The minimum sum of weights in each node.}
  \item{ptry}{positive scalar numeric. The number of components randomly
    sampled as candidates at each step. If \code{0 < ptry < 1},
    \code{ptry} is interpreted as the relative number of components to 
    select relative to the total number of components.}
  \item{ntry}{positive numeric, either a scalar or a vector of length
    equal the number of \code{vc} components. The number(s) of input
    nodes of each component randomly sampled as candidates at each
    step. If \code{0 < ntry < 1}, \code{ntry} is interpreted as the
    relative number of nodes to select relative to the total number of
    nodes in the corresponding component.} 
  \item{vtry}{positive numeric, either a scalar or a vector of length
    equal the number of \code{vc} components. The number(s) of input
    variables of each component randomly sampled as candidates at each
    step. If \code{0 < vtry < 1}, \code{ntry} is interpreted as the
    relative number of variables to select relative to the total number
    of variables for the corresponding component.} 
}

\details{
  Implements the \emph{bagging} (Breiman, 1996) and \emph{random
    forests} (Breiman, 2001) ensemble algorithms for the
  \code{\link{tvcm}} algorithm. The method consist in growing multiple 
  trees by using \code{\link{tvcm}} and aggregating the fitted
  coefficient functions. For bagging, simply set \code{vtry = Inf} in
  \code{\link{fvcm_control}}.  

  \code{\link{fvcolmm}} and \code{\link{fvcglm}} are convenience
  functions for whether a \code{\link{olmm}} or a \code{\link{glm}}
  model is fitted. 

  \code{\link{fvcm_control}} is a wrapper for \code{\link{tvcm_control}}
  with other defaults for the main tuning parameters. The principal
  tuning parameters for \code{\link{fvcm}} are \code{maxstep},
  \code{ptry}, \code{ntry} and \code{vtry}. \code{alpha} is set to 1 so
  that the algorithm is not stopped because of high p-values in the
  coefficient constancy tests (the parameter is ignored if \code{sctest
    = FALSE}). The parameter \code{maxstep} is discussed in the details
  section of \code{\link{tvcm}}. The arguments \code{ptry}, \code{ntry}
  and \code{vtry} determine random selections in each
  step. Specifically, \code{ptry} specifies the number of
  \code{\link{vc}} components to be randomly selected, \code{ntry} the
  number of nodes and \code{vtry} the number of variables. The latter
  two parameters can be specified componentwise. The \code{vtry} is
  commonly specified as the number of partitioning  variables divided
  by 3, see Hastie et al. (2001), sec. 15.
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
control <- fvcm_control(sctest = TRUE, vtry = 1, minsize = 25, verbose = TRUE)

## fit the forest
model <- fvcm(y ~ -1 + vc(z1, z2, by = x1, intercept = TRUE) + x2,
              data = vcrpart_2, family = gaussian(), 
              folds = folds_control("subsampling", 3, 0.75), control = control)
model
}

\keyword{models}