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

\title{Model assessment and model selection for \command{\link{tvcm}} objects.}

\description{
  Out-of-bag loss, cross-validation and pruning for 
  \command{\link{tvcm}} objects. 
}

\usage{

folds_control(type = c("kfold", "subsampling", "bootstrap"),
      K = ifelse(type == "kfold", 5, 30),
      prob = 0.5, weights = c("case", "freq"),
      seed = NULL)

\method{cvloss}{tvcm}(object, folds = folds_control(), ...)

\method{print}{cvloss.tvcm}(x, ...)

\method{plot}{cvloss.tvcm}(x, legend = TRUE, details = TRUE, ...)

\method{oobloss}{tvcm}(object, newdata = NULL, weights = NULL,
        fun = NULL, ...)

\method{prune}{tvcm}(tree, cp = NULL, alpha = NULL, maxstep = NULL,
      terminal = NULL, original = FALSE, ...)
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
    and the in-sample prediction error should be shown.}
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
    the subnodes of which should be merged. }
  \item{original}{logical scalar. Whether pruning should be based on the
    trees from partitioning rather than on the current trees.}
  \item{...}{other arguments to be passed.} 
}

\details{As described in the help of \command{\link{tvcm}}, TVCM
  (combined with \code{sctest = FALSE}) is a two stage procedure that
  first grows overly fine partitions and second selects the best-sized
  partitions by pruning. Both steps can be carried out with a single
  \command{\link{tvcm}} and several parameters can be specified with
  \command{\link{tvcm_control}}. The here presented functions may be
  interesting for advanced users who want to process the two stages by
  separate calls.  

  The \command{\link{prune}} method collapses inner nodes of the overly
  large tree fitted with \command{\link{tvcm}} according to the tuning
  parameter \code{cp}, the goal being to minimize the cost-complexity
  in-sample error. The cost-complexity in-sample error is, in what 
  follows, defined as the in-sample error plus \code{cp} times the
  complexity of the model. Specifically, in the considered likelihood
  setting, the in-sample error is equal - 2 times the log-likelihood at
  the estimated parameters. Further, the complexity of the model is
  defined as \code{dfpar} times the number of coefficients times
  \code{dfsplit} times the number of splits. By default (see
  \command{\link{tvcm_control}}), the complexity is just the number of
  splits since \code{dfpar = 0} and \code{dfpar = 0} and \code{dfsplit =
    1}. The complexity penalty \code{cp} is generally unknown and
  estimated by using cross-validation.  

  \command{\link{folds_control}} and \command{\link{cvloss}} allow for
  estimating \code{dfsplit} by cross-validation. The function
  \command{\link{folds_control}} is used to specify the cross-validation
  scheme, where a random 5-fold cross-validation scheme is set as the
  default. Alternatives are \code{type = "subsampling"} (random draws
  without replacement) and \code{type = "bootstrap"} (random draws with
  replacement). For 2-stage models (with random-effects) fitted by
  \command{\link{olmm}}, the subsets are based on subject-wise i.e. first
  stage sampling. For models where weights represent frequencies of
  observation units (e.g., data from contingency tables), the option
  \code{weights = "freq"} should be used.  

  \command{\link{cvloss}} repeatedly fits \command{\link{tvcm}} objects
  based on the internally created folds and evaluates mean out-of-bag
  loss of the model at different levels of the tuning parameter
  \code{dfsplit}. Out-of-bag loss refers here to the prediction error
  based on a loss function, which is typically the -2 log-likelihood 
  error (see the details for \code{oobloss} below). Commonly,
  \code{dfsplit} is used for backward pruning (\code{direction =
  "backward"}), but it is also possible to cross-validate \code{dfsplit}
  for premature stopping (\code{direction = "forward"}, see argument
  \code{dfsplit} in
  \command{\link{tvcm_control}}). \command{\link{cvloss}} returns an
  object for which a \code{print} and a \code{plot} generic is
  provided. The proposed estimate for \code{dfsplit} is the one that
  minimizes the validated loss and can be extracted from component
  \code{dfsplit.min}.   
  
  \command{\link{oobloss}} can be used for estimating the out-of-bag
  prediction error for out-of-bag data (the \code{newdata}
  argument). By default, the loss is defined as the sum of deviance
  residuals, see the return value \code{dev.resids} of
  \command{\link{family}} resp. \command{\link{family.olmm}}. Otherwise,
  the loss function can be defined manually by the argument \code{fun},
  see the examples below. In general the sum of deviance residual is
  equal the sum of the -2 log-likelihood errors. A special case is the
  gaussian family, where the deviance residuals are computed as
  \eqn{\sum_{i=1}^N w_i (y_i-\mu)^2}, that is, the deviance residuals
  ignore the term \eqn{\log{2\pi\sigma^2}}. Therefore, the sum of
  deviance residuals for the gaussian model (and possibly others) is not
  exactly the sum of -2 log-likelihood prediction errors (but shifted by
  a constant). Another special case are models with random effects. For
  models based on \command{\link{olmm}}, the deviance residuals are
  retrieved from marginal predictions (where random effects are
  integrated out). 

  The \command{\link{prune}} function is used to select a nested model
  of the current model, i.e. a model which collapses terminal nodes to
  their inner nodes based by minimizing the in-sample error (see
  above). Pruning works as follows: In each iteration, all models that
  collapse one inner node of the current models are fitted. The inner
  node that yields the smallest increase in the in-sample error per
  complexity unit (e.g., by the number of splits) is collapsed and the
  resulting model substitutes the current model. The algorithm is
  stopped as soon as the increase of in-sample error per complexity
  unit of several nested models is larger than \code{cp}, 


  all nested models
  have a higher estimated prediction error than the 
  current model, which will be returned. 
}

\value{\command{\link{folds_control}} returns a list of parameters for
  building a cross-validation scheme. \command{\link{cvloss}} returns an
  \code{cvloss.tvcm} object with the following essential components: 

  \item{grid}{a list with two matrices \code{dfsplit} and
    \code{nsplit}. Specifies the grid of values at which the
    cross-validated loss was evaluated.}  
  \item{loss}{a list with two matrices \code{dfsplit} and
    \code{nsplit}. The cross-validated loss for each fold corresponding
    to the values in \code{grid}.} 
  \item{dfsplit.min}{numeric scalar. The tuning parameter which
    minimizes the cross-validated loss.}
  \item{folds}{the used folds to extract the learning and the validation
    sets.}

  Further, \command{\link{oobloss}} returns the loss of the \code{newdata}
  validation set and \command{\link{prune}} returns a (possibly modified)
  \command{\link{tvcm}} object.  
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
