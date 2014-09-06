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

fvcm(..., control = fvcm_control())

fvcolmm(..., family = cumulative(), control = fvcm_control())

fvcglm(..., family, control = fvcm_control())

fvcm_control(maxstep = 10, folds = folds_control("subsampling", 5),
             ptry = 1, ntry = 1, vtry = 5, 
             alpha = 1.0, maxoverstep = Inf, ...)

}

\arguments{
  \item{\ldots}{for \code{\link{fvcm}}, \code{\link{fvcolmm}} and
    \code{\link{fvcglm}}, \code{\ldots} are arguments to be passed to
    \code{\link{tvcm}}, see \code{\link{tvcm}}. This includes at least
    the arguments \code{formula}, \code{data} and \code{family} (see
    examples). For \code{\link{fvcm_control}}, \code{\ldots} are further
    arguments to be passed to \code{\link{tvcm_control}}.}      
  \item{control}{control parameters as produced by
    \code{\link{fvcm_control}}.} 
  \item{family}{the model family, e.g., \code{\link{cumulative}} for
    fitting a cumulative link model.} 
  \item{maxstep}{integer. The maximum number of steps for when growing
    individual trees.}
  \item{folds}{a list of parameters for folds as created by
    \code{\link{folds_control}}.}
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
  \item{maxoverstep, alpha}{These two parameters are specified to disable
    the default stopping rules. See \code{\link{tvcm_control}} for details.}
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

  \code{\link{fvcm_control}} is a wrapper of \code{\link{tvcm_control}} and 
  the arguments indiciated specify modified defaults and parameters for 
  randomizing the split selection in random forests. Notice that, relative 
  to \code{\link{tvcm_control}}, also the \code{cv} \code{prune} arguments are 
  internally disabled. The default arguments for \code{alpha} and \code{maxoverstep} 
  essentially disable the stopping rules of \code{\link{tvcm}}, where the argument 
  \code{maxstep} fully controls the stopping. The three parameters \code{ptry}, 
  \code{ntry} and \code{vtry} control the randomization for selecting the  
  \code{\link{vc}} term, the node and the variable. The default of \code{vtry = 5} 
  should be adjusted in applications. For this, the number of partitioning variables 
  divided by 3 is a common rule-of-thumb, see Hastie et al. (2001).

  The performance of the implementation has not been investigated yet. For 
  computational reasons we also suggest to use a single \code{\link{vc}} term
  combined with \code{sctest = TRUE} in \code{\link{fvcm_control}}, see the 
  example below.
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
control <- fvcm_control(sctest = TRUE, maxstep = 2, 
                        vtry = 1, verbose = TRUE, minsize = 10,
                        folds = folds_control("subsampling", 3, 0.75))

## fit the forest
model <- fvcm(y ~ -1 + vc(z1, z2, by = x1, intercept = TRUE) + x2,
              data = vcrpart_2, family = gaussian(), 
              folds = folds_control("subsampling", 3, 0.75), control = control)
model


}

\keyword{models}