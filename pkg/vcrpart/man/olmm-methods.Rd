\name{olmm-methods}

\alias{olmm-methods}
\alias{anova.olmm}
\alias{coefficients.olmm}
\alias{coef.olmm}
\alias{deviance.olmm}
\alias{formula.olmm}
\alias{fixef}
\alias{fixef.olmm}
\alias{getCall.olmm}
\alias{logLik.olmm}
\alias{model.frame.olmm}
\alias{model.matrix.olmm}
\alias{ranef}
\alias{ranef.olmm}
\alias{ranefCov}
\alias{ranefCov.olmm}
\alias{resid.olmm}
\alias{residuals.olmm}
\alias{simulate.olmm}
\alias{terms.olmm}
\alias{update.olmm}
\alias{VarCorr}
\alias{VarCorr.olmm}
\alias{print.VarCorr.olmm}
\alias{vcov.olmm}
\alias{weights.olmm}

\title{Methods for \code{olmm} objects}

\description{Standard methods for computing on \code{olmm} objects.}

\usage{
\method{anova}{olmm}(object, ...)

\method{coef}{olmm}(object, which = c("all", "fe"), ...)

\method{fixef}{olmm}(object, which = c("all", "ce", "ge"), ...)

\method{model.matrix}{olmm}(object, which = c("fe", "fe-ce", "fe-ge",
             "re", "re-ce", "re-ge"), ...)

\method{ranef}{olmm}(object, norm = FALSE, ...)

\method{ranefCov}{olmm}(object, ...) 

\method{simulate}{olmm}(object, nsim = 1, seed = NULL,
         newdata = NULL, ranef = TRUE, ...)

\method{terms}{olmm}(x, which = c("fe-ce", "fe-ge", "re-ce", "re-ge"), ...)

\method{VarCorr}{olmm}(x, sigma = 1., rdig = 3)

\method{weights}{olmm}(object, level = c("observation", "subject"), ...)
}

\arguments{
  \item{object, x}{an \code{\link{olmm}} object.}
  \item{which}{optional character string. For  \code{\link{coef}} and
    \code{\link{fixef}}, it indicates whether \code{"all"} coefficients,
    the fixed effects \code{"fe"}, the category-specific fixed effects
    \code{"ce"} (i.e. non-proportional odds) or the global fixed effects
    \code{"ge"} (i.e. proportional odds) should be extracted. For
    \code{model.matrix} it indicates whether the model matrix of the
    fixed- (\code{"fe"}) or the random effects (\code{"re"}) should be
    extracted.}
  \item{level}{character string. Whether the results should be on the
    observation level (\code{level = "observation"}) or on the subject
    level (\code{level = "subject"}).} 
  \item{norm}{logical. Whether residuals should be divided by their
    standard deviation.}
  \item{nsim}{number of response vectors to simulate.  Defaults to 1.}
  \item{seed}{an object specifying if and how the random number
    generator should be initialized. See \code{\link{simulate}}}
  \item{newdata}{a data frame with predictor variables.}
  \item{ranef}{whether the simulated responses should be based on
    conditional probabilities.}
  \item{sigma}{ignored but obligatory argument from original generic.} 
  \item{rdig}{ignored but obligatory argument from original generic.}
  \item{...}{potential further arguments passed to methods.}
}

\details{
  \code{\link{anova}} implements log-likelihood ratio tests for model 
  comparisons, based on the marginal likelihood. At the time being, 
  at least two models must be assigned.

  \code{ranefCov} extracts the variance-covariance matrix of the random
  effects. Similarly, \code{VarCorr} extracts the estimated variances,
  standard deviations and correlations of the random effects. 

  \code{resid} extracts the residuals of Li and Sheperd (2012). By
  default, the marginal outcome distribution is used to compute these residuals. 
  The conditional residuals can be computed by assigning \code{ranef = TRUE} as 
   a supplementary argument.  
  
  Further, undocumented methods are \code{\link{deviance}},
  \code{\link{extractAIC}}, \code{\link{fitted}}, \code{\link{formula}},
  \code{\link{getCall}}, \code{\link{logLik}},
  \code{\link{model.frame}}, \code{\link{nobs}}, \code{\link{update}},
  \code{\link{vcov}}. 
}

\author{Reto Buergin}

\references{
  Li, C. and Sheperd, B. E. (2012). A new residual for ordinal
  outcomes, \emph{Biometrika} \bold{99} (2): 437-480
}

\seealso{\code{\link{olmm}}, \code{\link{olmm-class}},
  \code{\link{predict.olmm}}, \code{\link{gefp.olmm}}}

\examples{
## --------------------------------------------------------- #
## Example 1: Schizophrenia (see also example of 'olmm')
## --------------------------------------------------------- #

data(schizo)

schizo <- schizo[1:181,]
schizo$id <- droplevels(schizo$id)

## fit two alternative models for the 'schizo' data
model.0 <- olmm(imps79o ~ tx + sqrt(week) + re(1|id), schizo)
model.1 <- olmm(imps79o ~ tx + sqrt(week)+tx*sqrt(week)+re(1|id),schizo)

anova(model.0, model.1)
coef(model.1)
fixef(model.1)
head(model.matrix(model.1, "fe-ge"))
head(weights(model.1))
ranefCov(model.1)
head(resid(model.1))
head(simulate(model.1))
terms(model.1)
VarCorr(model.1)
head(weights(model.1, "subject"))
}

\keyword{methods}