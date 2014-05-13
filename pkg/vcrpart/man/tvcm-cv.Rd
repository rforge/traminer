\name{tvcm-cv}			

\alias{tvcm-cv}
\alias{AIC.tvcm}
\alias{BIC.tvcm}
\alias{cvfolds}
\alias{cvrisk}
\alias{cvrisk.tvcm}
\alias{risk}
\alias{risk.tvcm}
\alias{oobrisk}
\alias{oobrisk.tvcm}
\alias{stabpath}
\alias{stabpath.tvcm}

\title{Model assessment and selection for \code{\link{tvcm}} objects.}

\description{
  In-sample Akaike information criterion (AIC) and Bayesian information
  criterion (BIC), out-of-bag risk, cross-validated risk and variable selection 
  stability paths for \code{\link{tvcm}} objects.
}

\usage{

\method{AIC}{tvcm}(object, sub = FALSE, alpha = NULL, maxwidth = NULL, 
    k = 2,  verbose = FALSE, ...)

\method{BIC}{tvcm}(object, ...)

cvfolds(x, type = c("kfold", "subsampling", "bootstrap"),
        K = ifelse(type == "kfold", 10L, 50L),
        prob = 0.5)

\method{cvrisk}{tvcm}(object, folds = cvfolds(object, "kfold", 10),
       fun = NULL, alpha = NULL, maxwidth = NULL, verbose = FALSE, ...)

\method{oobrisk}{tvcm}(object, newdata = NULL, weights = NULL, 
     fun = NULL, ...)

\method{stabpath}{tvcm}(object, q, alpha = NULL, maxwidth = NULL, ...)

}

\arguments{
  \item{object}{an object of class \code{\link{tvcm}}.}
  \item{x}{either a \code{\link{tvcm}} object or a \code{factor}
    representing the stage-1 structure (subjects) of the data.}
  \item{sub}{logical scalar. Whether a list of AIC's of submodels should
    be computed.}
  \item{k}{numeric scalar. The penalty per parameter to be used; the
    default \code{k=2} is the classical AIC.}  
  \item{alpha}{numerical scalar. Maximal alpha parameter to be
    evaluated. See \code{\link{tvcm_control}}.} 
  \item{maxwidth}{integer scalar. Maximal width of the tree to be
    evaluated. See \code{\link{tvcm_control}}.} 
  \item{type}{character string. The type of sampling scheme to be
    used to extract subsets.}   
  \item{K}{integer. The number of folds to be evaluated.} 
  \item{prob}{numeric. The probability for the \code{"subsampling"} 
    cross-validation scheme.}
  \item{folds}{a integer matrix. User specified folds for
    cross-validation from which columns are used to extract subsets of
    the data.}
  \item{fun}{the risk function. By default, the mean of the deviance 
    residuals are used. \code{\link{family}}.} 
  \item{newdata}{a data.frame of out-of-bag data (including the response
    variable).}
  \item{weights}{a numeric vector of weights corresponding to entries in \code{newdata}.} 
  \item{q}{maximum number of variables to be selected in each
    evaluation.} 
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the validation.}
  \item{...}{other arguments to be passed.} 
}

\details{The \code{\link{tvcm}} algorithm provides with 
  \code{method = "mob"} an implementation that selects the size of 
  the tree automatically. An alternative to such premature stopping 
  is to evaluate the predictive performance at different levels of 
  the tuning parameter, which is \code{alpha} for \code{method = "mob"} 
  and \code{maxwidth} for \code{method = "partreg"}, selecting the 
  tuning parameters that performs best and re-fit the model using 
  the selected parameter. This is what you can do with \code{AIC}, 
  \code{BIC} and \code{vcrisk}.

  \code{AIC} resp. \code{BIC} with argument \code{sub = TRUE} show 
  the Akaike information criterion resp. Bayesian information criterion, 
  evaluated for different \code{alpha}s resp. \code{maxwidth}s. 
  If the assigned \code{alpha} resp. \code{maxwidth} exceeds the 
  corresponding parameter used for fitting, the model will be re-fitted. 
  
  \code{cvrisk} cross-validates the risk, which is by default the mean 
  of deviance residuals, for different levels of the principal tuning 
  parameter, i.e., for different levels of \code{alpha} for 
  \code{method = "mob"} and for different levels of\code{maxwidth} 
  for \code{method = "partreg"}. \code{cvrisk} requires an \eqn{n \times K} 
  integer matrix \code{folds} where \eqn{n} is equal \code{nobs(object)} 
  and \eqn{K} the number of desired evaluation steps. For the \eqn{k}-th 
  evaluation step, the \eqn{k}-th column is used to extract a subset of 
  the learning data to re-fit the model. Specifically, the integer 
  numbers of entries in \code{folds} define the number of times
  observations appear in the learning data for re-fitting the model. Zero 
  entries are used to estimate the out-of-bag error. We recommend to 
  produce the \code{folds} matrix by the utility function \code{cvfolds}.

  The default for \code{cvfolds} is to produce a random 10-fold 
  cross-validation scheme. Alternatives are \code{type = "subsampling"} 
  (random draws without replacement) and \code{type = "bootstrap"} 
  (random draws with replacement). For 2-stage models fitted by 
  \code{\link{olmm}}, the subsets are based on subject-wise i.e. 
  first stage sampling.
  
  \code{oobrisk} can be used for estimating the prediction error of the
  selected model based on \code{newdata}. By default, the risk is defined
  as the sum of deviance residuals. The risk function may be defined
  manually by the argument \code{fun}, see the examples below.
  
  The function \code{stabpath} implements the "stability paths" of
  Meinshausen and Buehlmann (2010). The function estimates the probability
  with which variables are selected by \code{\link{tvcm}}, by fitting
  the model on multiple subsamples. In each evaluation step, a tree
  is grown until the maximum tuning parameter is reached or the number
  of selected variables is equal \code{q}. For details, see 
  Meinshausen and Buehlmann (2010).
}

\value{An object of class \code{cvrisk.tvcm} or \code{stabpath.tvcm}.}

\references{
  Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) 
  \emph{Classification and Regression Trees}. Wadsworth.

  T. Hastie, R. Tibshirani, J. Friedman (2001), The elements of
  statistical learning, Springer.
  
  Meinshausen, N. and Buehlmann, P. (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}(4).
}

\seealso{
  \code{\link{tvcm}}
}

\examples{
## --------------------------------------------------------- #
## Dummy Example 1:
##
## Model selection for the 'vcrpart_2' data with 
## 'method = "partreg"'. The example is merely a syntax 
## template.
## --------------------------------------------------------- #

## load the data
data(vcrpart_2)

## fit the model
control <- tvcm_control(method = "partreg",  maxwidth = 2, minbucket = 1)

model <- tvcglm(y ~ vc(z1, z2, by = x1) + x2, data = vcrpart_2,
                family = gaussian(), control = control, subset = 1:75)

## AIC and BIC
AIC(model, sub = TRUE)
BIC(model, sub = TRUE)

## cross-validation
folds <- cvfolds(model, type = "kfold", K = 2)
cv <- cvrisk(model, folds, maxwidth = 3)
plot(cv) # use maxwidth = 2, as already specified

## out-of-bag error
oobrisk(model, newdata = vcrpart_2[76:100,])

## use an alternative risk function
rfun <- function(y, mu, wt) sum(abs(y - mu))
oobrisk(model, newdata = vcrpart_2[76:100,], fun = rfun)

## stability paths
folds <- cvfolds(model, type = "subsampling", K = 2)
sl <- stabpath(model, q = 1, folds = folds, maxwidth = 3)
plot(sl)
}

\author{Reto Buergin}

\keyword{validation}
