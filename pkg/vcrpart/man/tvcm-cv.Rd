\name{tvcm-cv}			

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
  criterion (BIC), out-of-bag sample risk, cross-validated risk and 
  variable selection stability paths of \code{\link{tvcm}} objects
  for varying tuning parameters.
}

\usage{

\method{AIC}{tvcm}(object, sub = FALSE, maxalpha = NULL, maxwidth = NULL, 
    k = 2,  verbose = FALSE, ...)

\method{BIC}{tvcm}(object, ...)

cvfolds(x, type = c("kfold", "subsampling", "bootstrap"),
        K = ifelse(type == "kfold", 10L, 50L),
        prob = 0.5)

\method{cvrisk}{tvcm}(object, folds = cvfolds(object, "kfold", 10),
       fun = NULL, maxalpha = NULL, maxwidth = NULL, verbose = FALSE, ...)

oobrisk(object, ...)

\method{oobrisk}{tvcm}(object, newdata = NULL, weights = NULL, 
     fun = NULL, ...)

stabpath(object, q, ...)

\method{stabpath}{tvcm}(object, q, maxalpha = NULL, maxwidth = NULL, ...)

}

\arguments{
  \item{object}{an object of class \code{\link{tvcm}}.}
  \item{x}{either a \code{\link{tvcm}} object or a \code{factor}
    representing the stage-1 structure (subjects) of the data.}
  \item{sub}{logical scalar. Whether a list of AIC's of submodels should
    be computed.}
  \item{k}{numeric scalar. The penalty per parameter to be used; the
    default \code{k=2} is the classical AIC.}  
  \item{maxalpha}{numerical scalar. Maximal alpha parameter to be
    evaluated.} 
  \item{maxwidth}{integer scalar. Maximal width of the tree to be
    evaluated.} 
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

\details{\code{AIC} and \code{BIC} with \code{sub = TRUE} may
  be used to select the principal tuning parameter, which is
  \code{alpha} (the significance threshold for coefficient constancy
  tests) for \code{method = "mob"} and \code{maxwidth} (the maximum
  number of nodes) for \code{method = "partreg"}, see
  \code{\link{tvcm_control}}. The function returns the values of
  these criteria for a increasing tuning parameters between
  \code{alpha=0} and \code{maxalpha} resp. between \code{width=1} and
  \code{maxwidth}. If \code{maxalpha} or \code{maxwidth} exceed the
  corresponding parameter used for fitting, the model will be
  re-fitted. 
  
  The \code{cvrisk} function aims at selecting the principal tuning
  parameters via cross-validation. Again, \code{maxalpha} and
  \code{maxwidth} are the maximum specifications
  considered. \code{cvrisk} requires a subset matrix \code{folds} of
  0's and 1's, which has the same number of rows than the learning
  data matrix and which may be produced by \code{cvfolds}. The default
  for \code{cvfolds} is to produce a random 10-fold cross-validation
  scheme. Alternatives are \code{type = "subsampling"} (random draws
  without replacement) and \code{type = "bootstrap"} (random draws with
  replacement). For 2-stage models fitted by \code{\link{olmm}}, the
  subsets are based on subject-wise i.e. first stage sampling.
  
  \code{oobrisk} can be used for estimating the prediction error of the
  selected model based on \code{newdata}.
  
  The function \code{stabpath} implements the "stability paths" of
  Meinshausen and Buehlmann (2010) for examining the variable selection.  
}

\value{An object of class \code{cvrisk.tvcm} or \code{stabpath.tvcm}.}

\references{
  Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) 
  \emph{Classification and Regression Trees}. Wadsworth.

  Meinshausen, N. and Buehlmann, P. (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}(4).
}

\seealso{
  \code{\link{tvcm}}
}

\examples{

\dontrun{
## Validate Example 2 in '?tvcm'

## load the data
data(tvcm_2)

## fit the model
mod <- tvcm(y ~ time + treat + (1|id) | z1 + z2 + z3, data = tvcm_2,
               control = tvcm_control(
                alpha = 0.05, bonferroni = FALSE),
                linear = "time")

## cross-validation
folds <- vcolmm_folds(mod, type = "kfold", K = 10)
cv <- cv.tvcm(mod, folds, control = list(alpha = 0.5))

## Overfitting is hardly a problem here.
plot(cv, log = "x")

## Stability paths
sl <- stabpath.tvcm(mod, q = 1, folds = folds,
                       control = list(alpha = 0.1))
plot(sl, log = "x")
}

}

\author{Reto Buergin}

\keyword{validation}
