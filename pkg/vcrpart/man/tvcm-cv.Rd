\name{tvcm-assessment}			

\alias{tvcm-assessment}
\alias{prune}
\alias{prune.tvcm}
\alias{folds_control}
\alias{cvloss}
\alias{cvloss.tvcm}
\alias{loss}
\alias{loss.tvcm}
\alias{oobloss}
\alias{oobloss.tvcm}
\alias{print.cvloss.tvcm}
\alias{plot.cvloss.tvcm}

\title{Model selection utility functions for \command{\link{tvcm}} objects.} 

\description{
  Pruning, cross-validation to find the optimal pruning parameter and computing 
  validation set errors for \command{\link{tvcm}} objects. 
}

\usage{

\method{prune}{tvcm}(tree, cp = NULL, alpha = NULL, maxstep = NULL,
      terminal = NULL, original = FALSE, ...)

folds_control(type = c("kfold", "subsampling", "bootstrap"),
      K = ifelse(type == "kfold", 5, 30),
      prob = 0.5, weights = c("case", "freq"),
      seed = NULL)

\method{cvloss}{tvcm}(object, folds = folds_control(), ...)

\method{print}{cvloss.tvcm}(x, ...)

\method{plot}{cvloss.tvcm}(x, legend = TRUE, details = TRUE, ...)

\method{oobloss}{tvcm}(object, newdata = NULL, weights = NULL,
        fun = NULL, ...)
}

\arguments{
  \item{object, tree}{an object of class \command{\link{tvcm}}.}
  \item{x}{an object of class \code{cvloss.tvcm} as produced by
    \command{\link{cvloss}}.} 
  \item{type}{character string. The type of sampling scheme to be used
    to divide the data of the input model in a learning and a validation
    set.}    
  \item{K}{integer scalar. The number of folds.} 
  \item{prob}{numeric between 0 and 1. The probability for the
    \code{"subsampling"} cross-validation scheme.}
  \item{weights}{for \command{\link{folds_control}}, a character that
    defines whether the weights of \code{object} are case weights or
    frequencies of cases; for \command{\link{oobloss}}, a numeric vector
    of weights corresponding to the rows of \code{newdata}.} 
  \item{seed}{an numeric scalar that defines the seed.}
  \item{folds}{a list with control arguments as produced by
    \command{\link{folds_control}}.}
  \item{legend}{logical scalar. Whether a legend should be added.}
  \item{details}{logical scalar. Whether the foldwise validation errors 
    should be shown.} 
  \item{fun}{the loss function for the validation sets. By default, the
    (possibly weighted) mean of the deviance residuals as defined by the
    \command{\link{family}} of the fitted \code{object} is applied.}
  \item{newdata}{a data.frame of out-of-bag data (including the response
    variable). See also \command{\link{predict.tvcm}}.}
  \item{cp}{numeric scalar. The complexity parameter to be cross-validated 
    resp. the penalty with which the model should be pruned.}
  \item{alpha}{numeric significance level. Represents the stopping
    parameter for \command{\link{tvcm}} objects grown with
    \code{sctest = TRUE}, see \command{\link{tvcm_control}}. A node is
    splitted when the \eqn{p} value for any coefficient stability test
    in that node falls below \code{alpha}.} 
  \item{maxstep}{integer. The maximum number of steps of the algorithm.} 
  \item{terminal}{a list of integer vectors with the ids of the nodes
    the inner nodes to be set to terminal nodes. The length of the list must 
    be equal the number of partitions.}
  \item{original}{logical scalar. Whether pruning should be based on the
    trees from partitioning rather than on the current trees.}
  \item{...}{other arguments to be passed.} 
}

\details{By default, \command{\link{tvcm}} is a two stage procedure that 
  first grows overly large trees and second selects the best-sized trees by 
  pruning. The here presented functions may be interesting for advanced users 
  who want to process the model selection stage separately.

  In normal practice, the \command{\link{prune}} function is used to collapse 
  inner nodes of the tree structures by the tuning parameter \code{cp}. The 
  aim of pruning by \code{cp} is to collapse inner nodes to minimize the 
  cost-complexity criterion

  \deqn{error(cp) = error(tree) + cp * complexity(tree)}
 
  whereby, the training error \eqn{error(tree)} is defined by \code{lossfun} 
  and \eqn{complexity(tree)} is defined as the total number of coefficients times
  \code{dfpar} plus the total number of splits times \code{dfsplit}. The function
  \code{lossfun} and the parameters \code{dfpar} and \code{dfsplit} are defined 
  by the \code{control} argument of \command{\link{tvcm}}, see also 
  \command{\link{tvcm_control}}. By default, \eqn{error(tree)} is minus two 
  times the total likelihood of the model and \eqn{complexity(tree)} the number 
  of splits. The minimization of \eqn{error(cp)} is implemented by the 
  following iterative backward-stepwise algorithm 

  \enumerate{
  \item fit all \code{subtree} models that collapse one inner node of the 
    current \code{tree} model.
  \item compute the per-complexity increase in the training error
    \deqn{dev = (error(subtree) - error(tree)) / 
        (complexity(tree) - complexity(subtree))}
    for all fitted \code{subtree} models
  \item if any \code{dev} < \code{cp} then set as the \code{tree} model
    the \code{subtree} that minimizes \code{dev} and repeated 1 to 3, 
    otherwise stop.
  }

  The penalty \code{cp} is generally unknown and is estimated adaptively from 
  the data. \command{\link{cvloss}} implements the cross-validation 
  method to do this. \command{\link{cvloss}} repeats for each fold the following 
  steps
 
  \enumerate{
  \item fit a new model with \command{\link{tvcm}} based on the training data of 
    the fold. 
  \item prune the new model for increasing \code{cp}. Compute for each 
    \code{cp} the average validation error.
  }

  Doing so yields for each fold a  sequence of values for \code{cp} and a 
  sequence of average validation errors. The obtained sequences for \code{cp} 
  are combined to a fine grid and the average validation error is averaged 
  correspondingly. From these two sequences we choose the \code{cp} that 
  minimizes the validation error. Notice that the average validation error
  is computed as the total prediction error of the validation set divided  
  by the sum of validation set weights. See also the argument \code{ooblossfun} in 
  \command{\link{tvcm_control}} and the function \command{\link{oobloss}}.

  The function \command{\link{folds_control}} is used to specify the 
  cross-validation scheme, where a random 5-fold cross-validation scheme 
  is set as the default. Alternatives are \code{type = "subsampling"} 
  (random draws without replacement) and \code{type = "bootstrap"} (random 
  draws with replacement). For 2-stage models (with random-effects) 
  fitted by \command{\link{olmm}}, the subsets are based on subject-wise 
  i.e. first stage sampling. For models where weights represent frequencies
  of observation units (e.g., data from contingency tables), the option
  \code{weights = "freq"} should be considered. \command{\link{cvloss}} 
  returns an object for which a \code{print} and a \code{plot} generic is
  provided.  
  
  \command{\link{oobloss}} can be used to estimate the total prediction error 
  for validation data (the \code{newdata} argument). By default, the loss is 
  defined as the sum of deviance residuals, see the return value \code{dev.resids} 
  of \command{\link{family}} resp. \command{\link{family.olmm}}. Otherwise, the 
  loss function can be defined manually by the argument \code{fun}, see the 
  examples below. In general the sum of deviance residual is equal the sum of 
  the -2 log-likelihood errors. A special case is the gaussian family, where 
  the deviance residuals are computed as \eqn{\sum_{i=1}^N w_i (y_i-\mu)^2}, 
  that is, the deviance residuals ignore the term \eqn{log 2\pi\sigma^2}. 
  Therefore, the sum of deviance residuals for the gaussian model (and 
  possibly others) is not exactly the sum of -2 log-likelihood prediction 
  errors (but shifted by a constant). Another special case are models with
  random effects. For models based on \command{\link{olmm}}, the deviance 
  residuals are retrieved from marginal predictions (where random effects are
  integrated out). 
}

\value{\command{\link{prune}} returns a \command{\link{tvcm}} object, 
  \command{\link{folds_control}} returns a list of parameters for building a 
  cross-validation scheme. \command{\link{cvloss}} returns an \code{cvloss.tvcm} 
  object with at least the following components: 

  \item{grid}{a list with values for \code{cp}.} 
  \item{oobloss}{a matrix recording the validated loss for each value in 
    \code{grid} for each fold.}
  \item{cp.hat}{numeric scalar. The tuning parameter which
    minimizes the cross-validated error.}
  \item{folds}{the used folds to extract the learning and the validation
    sets.}

  \command{\link{oobloss}} returns a scalar representing the total prediction 
  error for \code{newdata}.
}

\references{
  Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) 
  \emph{Classification and Regression Trees}. Wadsworth.

  T. Hastie, R. Tibshirani, J. Friedman (2001), The elements of
  statistical learning, Springer.
}

\seealso{
  \command{\link{tvcm}}
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
control <- tvcm_control(maxstep = 2L, minsize = 5L, cv = FALSE)
model <- tvcglm(y ~ vc(z1, z2, by = x1) + vc(z1, by = x2),
                data = vcrpart_2, family = gaussian(),
                control = control, subset = 1:75)

## cross-validate 'dfsplit'
cv <- cvloss(model, folds = folds_control(type = "kfold", K = 2, seed = 1))
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
