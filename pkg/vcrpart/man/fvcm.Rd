\name{fvcm}

\alias{fvcm}
\alias{fvcm_control}
\alias{fvcolmm}
\alias{fvcolmm_control}
\alias{fvcglm}
\alias{fvcglm_control}

\title{Bagging and Random Forests based on \command{\link{tvcm}}}

\description{
  Bagging (Breiman, 1996) and Random Forest (Breiman, 2001) ensemble
  algorithms for \command{\link{tvcm}}.
}

\usage{

fvcm(..., control = fvcm_control())

fvcm_control(maxstep = 10, folds = folds_control("subsampling", 5),
             ptry = 1, ntry = 1, vtry = 5, 
             alpha = 1.0, mindev = 0.0, ...)

fvcolmm(..., family = cumulative(), control = fvcolmm_control())

fvcolmm_control(maxstep = 10, folds = folds_control("subsampling", 5),
                ptry = 1, ntry = 1, vtry = 5, alpha = 1.0, ...)

fvcglm(..., family, control = fvcglm_control())

fvcglm_control(maxstep = 10, folds = folds_control("subsampling", 5),
               ptry = 1, ntry = 1, vtry = 5, mindev = 0, ...)
}

\arguments{
  \item{\ldots}{for \command{\link{fvcm}}, \command{\link{fvcolmm}} and
    \command{\link{fvcglm}} arguments to be passed to
    \command{\link{tvcm}}. This includes at least the arguments
    \code{formula}, \code{data} and \code{family}, see examples below. For
    \command{\link{fvcm_control}} further control arguments to be passed
    to \command{\link{tvcm_control}}.}       
  \item{control}{a list of control parameters as produced by
    \command{\link{fvcm_control}}.} 
  \item{family}{the model family, e.g., \command{\link{binomial}} or 
    \command{\link{cumulative}}.} 
  \item{maxstep}{integer. The maximum number of steps for when growing
    individual trees.}
  \item{folds}{a list of parameters to control the extraction of subsets, 
    as created by \command{\link{folds_control}}.}
  \item{ptry}{positive numeric scalar. The number of \command{\link{vc}}
    terms to be randomly sampled as candidates in each iteration. If
    \code{0 < ptry < 1}, \code{ptry} is interpreted as the relative number
    of terms to select regarding the total number of terms.} 
  \item{ntry}{positive numeric, either a scalar or a vector of length
    equal the number of \command{\link{vc}} terms. The number(s) of
    nodes of each term to be randomly sampled as candidates in each
    iteration. If \code{0 < ntry < 1}, \code{ntry} is interpreted as the
    relative number of nodes to select regarding the current total number 
    of nodes in the term.} 
  \item{vtry}{positive numeric, either a scalar or a vector of length
    equal the number of \command{\link{vc}} terms. The number(s) of input
    variables of each term to be randomly sampled as candidates in each
    iteration. If \code{0 < vtry < 1}, \code{ntry} is interpreted as the
    relative number of variables to select regarding the total number of
    variables of the corresponding term.} 
  \item{mindev, alpha}{These two parameters are merely specified to
    disable the default stopping rules for \command{\link{tvcm}}. See
    also \command{\link{tvcm_control}} for details.} 
}

\details{
  Implements the \emph{Bagging} (Breiman, 1996) and \emph{Random
  Forests} (Breiman, 2001) ensemble algorithms for
  \command{\link{tvcm}}. The method consist in growing multiple trees by
  using \command{\link{tvcm}} and aggregating the fitted coefficient
  functions in the scale of the predictor function. To enable bagging,
  use \code{ptry = Inf}, \code{ntry = Inf} and \code{vtry = Inf} in
  \command{\link{fvcm_control}}.    

  \command{\link{fvcolmm}} and \command{\link{fvcglm}} are the
  extensions for \command{\link{tvcolmm}} and
  \command{\link{tvcglm}}.

  \command{\link{fvcm_control}} is a wrapper of
  \command{\link{tvcm_control}} and the arguments indicated specify
  modified defaults and parameters for randomizing split
  selections. Notice that, relative to \command{\link{tvcm_control}},
  also the \code{cv} \code{prune} arguments are internally disabled. The
  default arguments for \code{alpha} and \code{maxoverstep} essentially
  disable the stopping rules of \command{\link{tvcm}}, where the
  argument \code{maxstep} (the number of iterations i.e. the maximum
  number of splits) fully controls the stopping. The three parameters
  \code{ptry}, \code{ntry} and \code{vtry} control the randomization for
  selecting the \command{\link{vc}} term, the node and the variable for
  splitting. The default of \code{vtry = 5} is arbitrary. It should be
  adjusted in applications, e.g., to the number of partitioning
  variables (for each term) divided by 3, see Hastie et al. (2001). 
}

\value{An object of class \code{fvcm}.}
\author{Reto Buergin}

\references{
  Breiman, L. (1996). Bagging Predictors. \emph{Machine Learning},
  \bold{24}(2), 123--140.
  
  Breiman, L. (2001). Random Forests. \emph{Machine Learning},
  \bold{45}(1), 5--32.

  Hastie, T., R. Tibshirani and J. Friedman (2001), The Elements of 
  Statistical Learning (2 ed.), Springer-Verlag.
}

\seealso{\command{\link{fvcm-methods}}, \command{\link{tvcm}}, 
  \command{\link{glm}}, \command{\link{olmm}}}

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
control <-
  fvcm_control(maxstep = 3, minsize = 10,
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