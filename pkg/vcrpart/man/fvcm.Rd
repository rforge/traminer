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
    \code{\link{fvcglm}} arguments to be passed to \code{\link{tvcm}}. 
    This includes at least the arguments \code{formula}, \code{data} and 
    \code{family}, see examples below. For \code{\link{fvcm_control}} 
    further control arguments to be passed to \code{\link{tvcm_control}}.}      
  \item{control}{a list of control parameters as produced by
    \code{\link{fvcm_control}}.} 
  \item{family}{the model family, e.g., \code{\link{binomial}} or 
    \code{\link{cumulative}}.} 
  \item{maxstep}{integer. The maximum number of steps for when growing
    individual trees.}
  \item{folds}{a list of parameters to control the extraction of subsets, 
    as created by \code{\link{folds_control}}.}
  \item{ptry}{positive numeric scalar. The number of \code{\link{vc}} terms
    to be randomly sampled as candidates in each iteration. If 
    \code{0 < ptry < 1}, \code{ptry} is interpreted as the relative number 
    of terms to select regarding the total number of terms.}
  \item{ntry}{positive numeric, either a scalar or a vector of length
    equal the number of \code{\link{vc}} terms. The number(s) of
    nodes of each term to be randomly sampled as candidates in each
    iteration. If \code{0 < ntry < 1}, \code{ntry} is interpreted as the
    relative number of nodes to select regarding the current total number 
    of nodes in the term.} 
  \item{vtry}{positive numeric, either a scalar or a vector of length
    equal the number of \code{\link{vc}} terms. The number(s) of input
    variables of each term to be randomly sampled as candidates in each 
    iteration. If \code{0 < vtry < 1}, \code{ntry} is interpreted as the
    relative number of variables to select regarding the total number
    of variables of the corresponding term.}
  \item{maxoverstep, alpha}{These two parameters are merely specified to 
    disable the default stopping rules for \code{\link{tvcm}}. See also 
    \code{\link{tvcm_control}} for details.}
}

\details{
  Implements the \emph{bagging} (Breiman, 1996) and \emph{random
   forests} (Breiman, 2001) ensemble algorithms for \code{\link{tvcm}}. The 
  method consist in growing multiple trees by using \code{\link{tvcm}} and
  aggregating the fitted coefficient functions. To enable bagging, use 
  \code{ptry = Inf}, \code{ntry = Inf} and \code{vtry = Inf} in 
  \code{\link{fvcm_control}}.  

  \code{\link{fvcolmm}} and \code{\link{fvcglm}} are convenience
  functions for whether a \code{\link{olmm}} or a \code{\link{glm}}
  model is fitted. 

  \code{\link{fvcm_control}} is a wrapper of \code{\link{tvcm_control}} and 
  the arguments indicated specify modified defaults and parameters for 
  randomizing split selections. Notice that, relative to \code{\link{tvcm_control}}, 
  also the \code{cv} \code{prune} arguments are internally disabled. The 
  default arguments for \code{alpha} and \code{maxoverstep} essentially disable 
  the stopping rules of \code{\link{tvcm}}, where the argument \code{maxstep} 
  (the number of iterations i.e. the maximum number of splits) fully controls the 
  stopping. The three parameters \code{ptry}, \code{ntry} and \code{vtry} control 
  the randomization for selecting the \code{\link{vc}} term, the node and the 
  variable for splitting. The default of \code{vtry = 5} is arbitrary. It should 
  be adjusted in applications, e.g., to the number of partitioning variables 
  (for each term) divided by 3, see Hastie et al. (2001).
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
## Bagging 'tvcm' on the artificially generated data 'vcrpart_3'. The 
## true coefficient function is a sinus curve between -pi/2 and pi/2. 
## The parameters 'maxstep = 3' and 'K = 5' are chosen to restrict the 
## computations.
## ------------------------------------------------------------------- #

## simulated data
data(vcrpart_3)

## setting parameters
control <- fvcm_control(maxstep = 3, minsize = 10,
                        folds = folds_control("subsampling", K = 5, 0.5, seed = 3))

## fitting the forest
model <- fvcm(y ~ vc(z1, by = x1), data = vcrpart_3, 
              family = gaussian(), control = control)

## plot the first two trees
plot(model, "coef", 1:2)

## plotting the partial dependency of the coefficient for 'x1'
plot(model, "partdep")
}

\keyword{models}