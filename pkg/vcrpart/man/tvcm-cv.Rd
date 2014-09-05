\name{tvcm-assessment}			

\alias{tvcm-assessment}
\alias{folds_control}
\alias{cvloss}
\alias{cvloss.tvcm}
\alias{loss}
\alias{loss.tvcm}
\alias{oobloss}
\alias{oobloss.tvcm}
\alias{print.cvloss.tvcm}
\alias{plot.cvloss.tvcm}
\alias{prune}
\alias{prune.tvcm}

\title{Model assessment and model selection for \code{\link{tvcm}} objects.}

\description{
  Out-of-bag loss, cross-validation and pruning for 
  \code{\link{tvcm}} objects. 
}

\usage{

folds_control(type = c("kfold", "subsampling", "bootstrap"),
      K = ifelse(type == "kfold", 5, 30),
      prob = 0.5, weights = c("case", "freq"),
      seed = NULL)

\method{cvloss}{tvcm}(object, folds = folds_control(),
       fun = NULL, dfpar = 2, direction = c("backward", "forward"),
       papply = mclapply, verbose = FALSE, ...)

\method{print}{cvloss.tvcm}(x, ...)

\method{plot}{cvloss.tvcm}(x, legend = TRUE, details = TRUE, ...)

\method{oobloss}{tvcm}(object, newdata = NULL, weights = NULL,
        fun = NULL, ...)

\method{prune}{tvcm}(tree, dfsplit = NULL, dfpar = 2.0,
      direction = c("backward", "forward"),
      alpha = NULL, maxstep = NULL, terminal = NULL,
      papply = mclapply, keeploss = FALSE, ...)
}

\arguments{
  \item{object, tree}{an object of class \code{\link{tvcm}}.}
  \item{x}{an object of class \code{cvloss.tvcm} as produced by
    \code{\link{cvloss}}.} 
  \item{type}{character string. The type of sampling scheme to be used
    to divide the data of the input model in a learning and a validation
    set.}    
  \item{K}{integer scalar. The number of folds.} 
  \item{prob}{numeric between 0 and 1. The probability for the
    \code{"subsampling"} cross-validation scheme.}
  \item{weights}{for \code{\link{folds_control}}, a character that defines
    whether the weights of \code{object} are case weights or frequencies
    of cases; for \code{\link{oobloss}}, a numeric vector
    of weights corresponding to the rows of \code{newdata}.}
  \item{seed}{an numeric scalar that defines the seed.}
  \item{folds}{a list with control arguments as produced by
    \code{\link{folds_control}}.}
  \item{fun}{the loss function for the validation sets. By default, the
    (possibly weighted) mean of the deviance residuals as defined by the
    \code{\link{family}} of the fitted \code{object} is applied.}
  \item{dfpar}{a numeric scalar larger than zero. The per parameter penalty
    to be applied. If the 2 log-likelihood prediction error is used,
    this value is typically set to 2.} 
  \item{direction}{either \code{"backward"} (the default) or
    \code{"forward"}. Indicates the pruning algorithm to be
    used. \code{"backward"} applies backward pruning where in each
    iteration the inner node that produces the smallest per-node
    increase is collapsed. \code{"forward"} applies premature stopping
    when no new split would improve the model.}
  \item{papply}{(parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. To run \code{cvloss} sequentially
    (i.e. not in parallel), use \code{\link{lapply}}.}
  \item{newdata}{a data.frame of out-of-bag data (including the response
    variable). See also \code{\link{predict.tvcm}}.}
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the validation.}
  \item{legend}{logical scalar. Whether a legend should be added.}
  \item{details}{logical scalar. Whether the foldwise validation errors
    and the in-sample prediction error should be shown.}
  \item{dfsplit}{numeric scalar. The per-split cost \code{dfsplit} with
    which the partitions are to be cross-validated. If no \code{dfsplit}
    is specified (default), the parameter is ignored for pruning.}
  \item{alpha}{numeric significance level. Represents the stopping
    parameter for \code{\link{tvcm}} objects grown with \code{sctest =
      TRUE}, see \code{\link{tvcm_control}}. A node is splitted when the 
    \eqn{p} value for any coefficient stability test in that node falls
    below \code{alpha}.}
  \item{maxstep}{integer. The maximum number of steps of the algorithm.}
  \item{terminal}{a list of integer vectors with the ids of the nodes
    the subnodes of which should be merged. }
  \item{keeploss}{logical scalar or numeric. If and how many times the
    statistics should be reused in following iterations. Specifically,
    the option activates approximating the AIC reduction per split based
    on AIC reduction of the last iteration adjusted by AIC reduction of
    the last collapse. Values different from \code{FALSE} or 0 will
    accelerate the algorithm at the expense of accuracy. In cases of
    multiple \code{\link{vc}} terms, only the statistics corresponding
    to the \code{\link{vc}} term in which a node was collapsed in the previous
    iteration are considered to be reused.}
  \item{...}{other arguments to be passed.} 
}

\details{As described in the help of \code{\link{tvcm}}, the TVCM is a two stage  
  procedure where first a model with overly fine partitions is grown 
  with \code{\link{tvcm}} and second the best-sized partitions are 
  selected with the \code{\link{prune}} method. An alternative is to specify
  \code{sctest = TRUE} in \code{\link{tvcm_control}} and select the 
  partition-sizes with the \code{alpha} and \code{bonferroni} parameters.

  The \code{\link{prune}} method collapses inner nodes of the overly large tree
  fitted with \code{\link{tvcm}} to minimize the in-sample prediction error. 
  The in-sample prediction error is, in what follows, defined as the mean of the
  in-sample loss plus \code{dfpar} times the number of coefficients plus 
  \code{dfsplit} times the number of splits. In the common likelihood setting, 
  the loss is equal - 2 times the maximum likelihood and \code{dfpar = 2}. The
  per-split penalty \code{dfsplit} generally unknown and estimated by using 
  cross-validation.

  \code{\link{folds_control}} and \code{\link{cvloss}} allow for estimating 
  \code{dfsplit} by cross-validation. The function \code{\link{folds_control}} 
  is used to specify the cross-validation scheme, where a random 5-fold 
  cross-validation scheme is set as the default. Alternatives are 
  \code{type = "subsampling"} (random draws without replacement) and 
  \code{type = "bootstrap"} (random draws with replacement). For 2-stage 
  models (with random-effects) fitted by \code{\link{olmm}}, the subsets 
  are based on subject-wise i.e. first stage sampling. For models where weights 
  represent frequencies of observation units with exactly the same values in 
  all variables (e.g., data from contingency tables), the option \code{weights =
    "freq"} should be used. 

  \code{\link{cvloss}} repeatedly fits \code{\link{tvcm}} objects based
  on the internally created folds and evaluates mean out-of-bag loss of
  the model at different levels of the tuning parameter \code{dfsplit}. 
  Out-of-bag loss refers here to the prediction error based on a loss function, 
  which is typically the -2 log-likelihood error (see the details for 
  \code{oobloss} below). Commonly, \code{dfsplit} is used for backward pruning 
  (\code{direction = "backward"}), but it is also possible to
  cross-validate \code{dfsplit} for premature stopping (\code{direction
  = "forward"}, see argument \code{dfsplit} in
  \code{\link{tvcm_control}}). \code{\link{cvloss}} returns an object
  for which a \code{print} and a \code{plot} generic is provided. The
  proposed estimate for \code{dfsplit} is the one that minimizes the
  validated loss and can be extracted from component \code{dfsplit.min}. 
  
  \code{\link{oobloss}} can be used for estimating the out-of-bag
  prediction error for out-of-bag data (the \code{newdata}
  argument). By default, the loss is defined as the sum of deviance
  residuals, see the return value \code{dev.resids} of
  \code{\link{family}} resp. \code{\link{family.olmm}}. Otherwise, the
  loss function can be defined manually by the argument \code{fun}, see 
  the examples below. In general the sum of deviance residual is equal
  the -2 log-likelihood. A special case is the gaussian family, where
  the deviance residuals are computed as \eqn{\sum_{i=1}^N w_i
  (y_i-\mu)^2} that is, the deviance residuals ignore the term
  \eqn{\log{2\pi\sigma^2}}. Therefore, the sum of deviance residuals for
  the gaussian model (and possibly others) is not exactly the -2
  log-likelihood prediction error but shifted by a constant. Another
  special case are models with random effects. For models based on
  \code{\link{olmm}}, the deviance residuals are based on the marginal
  predictions (where random effects are integrated out).

  The \code{\link{prune}} function is used to select a nested model
  of the current model, i.e. a model which collapses leaves to their
  inner nodes based on the estimated prediction error. The estimated
  prediction error is defined as the AIC of the model plus \code{dfsplit}
  times the number of splits. Pruning with \code{direction = "backward"}
  works as follows: In each iteration, all nested models of the current
  model are evaluated, i.e. models which collapse one of the inner nodes
  of the current model. The inner node that yields the smallest increase
  in the estimated prediction error is collapsed and the resulting model
  substitutes the current model. The algorithm is stopped as soon as all
  nested models have a higher estimated prediction error than the
  current model, which will be returned. 
}

\value{\code{\link{folds_control}} returns a list of parameters for building a
  cross-validation scheme. \code{\link{cvloss}} returns an
  \code{cvloss.tvcm} object with the following essential components: 

  \item{grid}{a list with two matrices \code{dfsplit} and
    \code{nsplit}. Specifies the grid of values at which the
    cross-validated loss was evaluated.}  
  \item{loss}{a list with two matrices \code{dfsplit} and
    \code{nsplit}. The cross-validated loss for each fold corresponding
    to the values in \code{grid}.}
  \item{dfsplit.min}{numeric scalar. The tuning parameter which minimizes
    the cross-validated loss.}
  \item{folds}{the used folds to extract the learning and the validation
    sets.}

  Further, \code{\link{oobloss}} returns the loss of the \code{newdata}
  validation set and \code{\link{prune}} returns a (possibly modified)
  \code{\link{tvcm}} object.  
}

\references{
  Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) 
  \emph{Classification and Regression Trees}. Wadsworth.

  T. Hastie, R. Tibshirani, J. Friedman (2001), The elements of
  statistical learning, Springer.
}

\seealso{
  \code{\link{tvcm}}
}

\examples{
## --------------------------------------------------------- #
## Dummy Example 1:
##
## Model selection for the 'vcrpart_2' data. The example is
## merely a syntax template.
## --------------------------------------------------------- #

## load the data
data(vcrpart_2)

## fit the model
control <- tvcm_control(maxstep = 2L, minsize = 5L)
model <- tvcglm(y ~ vc(z1, z2, by = x1) + vc(z1, by = x2), data = vcrpart_2,
                family = gaussian(), control = control, subset = 1:75)

## cross-validate 'dfsplit'
set.seed(1)
cv <- cvloss(model, folds_control(type = "kfold", K = 2))
cv
plot(cv)

## out-of-bag error
oobloss(model, newdata = vcrpart_2[76:100,])

## use an alternative loss function
rfun <- function(y, mu, wt) sum(abs(y - mu))
oobloss(model, newdata = vcrpart_2[76:100,], fun = rfun)
}

\author{Reto Buergin}

\keyword{validation}
