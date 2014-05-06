\docType{class}
\name{olmm-class}
\alias{olmm-class}
\alias{anova.olmm}
\alias{coefficients.olmm}
\alias{coef.olmm}
\alias{deviance.olmm}
\alias{fitted.olmm}
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
\alias{reweight.olmm}
\alias{show,olmm-method}
\alias{terms.olmm}
\alias{update.olmm}
\alias{VarCorr}
\alias{VarCorr.olmm}
\alias{print.VarCorr.olmm}
\alias{vcov.olmm}
\alias{weights.olmm}

\title{Class \code{olmm}}

\description{
The \code{olmm} class represents 2-stage ordinal linear mixed models
fit by the \code{\link{olmm}} function.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("olmm", ...)}. 
}

\section{Slots}{
  \describe{
    \item{\code{env}:}{environment in which the object was built (class
      \code{"environment"}).} 
    \item{\code{frame}:}{the model frame.}
    \item{\code{call}:}{the matched call to the function that created
      the object (class \code{"call"}).}
    \item{\code{control}:}{a control object produced by
      \code{\link{olmm_control}}.} 
    \item{\code{formula}:}{the formula of the call (class
      \code{"formula"}).}
    \item{\code{terms}}{a list of \code{terms} of the fitted 
      model. The first element corresponds to \code{X}, the second to
      \code{W}.}
    \item{\code{family}:}{an object of class \code{family.olmm}} 
    \item{\code{y}:}{ordered response vector.} 
    \item{\code{X}:}{model matrix for the fixed effects.}
    \item{\code{W}:}{model matrix for the random effects.}
    \item{\code{subject}:}{a factor vector with grouping levels.}
    \item{\code{subjectName}:}{variable name of the subject vector.}
    \item{\code{weights}:}{numeric weights vector of length n.}
    \item{\code{weights_sbj}:}{numeric weights vector of length N.}
    \item{\code{offset}:}{numeric offset matrix}
    \item{\code{xlevels}:}{(only where relevant) a list of levels of the
      factors used in fitting.} 
    \item{\code{contrasts}:}{(only where relevant) a list of contrasts
      used.} 
    \item{\code{dims}:}{a named integer of dimensions. Some of the
      dimensions are: \eqn{n} is the number of observations, \eqn{p} is
      the number of fixed effects per predictor and \eqn{q} is the total
      number of random effects.}
    \item{\code{fixef}:}{a matrix of fixed effects (one column for each
      predictor).}
    \item{\code{ranefCholFac}:}{a lower triangular matrix. The cholesky
      decomposition of the covariance matrix of the random effects.}
    \item{\code{coefficients}:}{a numeric vector of several fitted model
      parameters}
    \item{\code{restricted}:}{a logical vector indicating which elements
      of the \code{coefficients} slot are restricted to an initial value
      at the estimation.}
    \item{\code{eta}:}{a matrix of unconditional linear predictors of
      the fixed effects without random effects.}  
    \item{\code{u}:}{a matrix of orthogonal standardized random
      effects (one row for each subject level).} 
    \item{\code{logLik_obs}:}{a numeric vector of log likelihood value
      (one value for each observation).}
    \item{\code{logLik_sbj}:}{a numeric vector of log likelihood values
      (one value for each subject level).}
    \item{\code{logLik}:}{a numeric value. The log likelihood of the
      model.} 
    \item{\code{score_obs}:}{a matrix of partial derivates of
      the log-likelihood equation (one row for each observation).} 
    \item{\code{score_sbj}:}{a matrix of partial derivates of the
      log-likelihood equation (one row for each subject level).} 
    \item{\code{score}:}{a numeric vector of (total) partial derivates
      of the log-Likelihood function.}
    \item{\code{info}:}{the information matrix (default is the expected
      information).}
    \item{\code{ghx}:}{a matrix of quadrature points for the
      Gauss-Hermite quadrature integration.} 
    \item{\code{ghw}:}{a matrix of weights for the Gauss-Hermite
      quadrature integration.} 
    \item{\code{ranefElMat}:}{a transformation matrix}
    \item{\code{optim}:}{a list of arguments for calling the optimizer
      function.} 
    \item{\code{output}:}{the output of the optimizer (class
      \code{"list"}).} 
  }
}

\usage{

\method{fixef}{olmm}(object, which = c("all", "ce", "ge"), ...)

\method{model.matrix}{olmm}(object, which = c("fe", "fe-ce", "fe-ge",
             "re", "re-ce", "re-ge"), ...)

\method{weights}{olmm}(object, level = c("observation", "subject"), ...)

\method{resid}{olmm}(object, norm = FALSE, ...)

\method{terms}{olmm}(x, which = c("fe-ce", "fe-ge", "re-ce", "re-ge"), ...)

\method{weights}{olmm}(object, level = c("observation", "subject"), ...)
}

\arguments{
  \item{object, x}{an \code{\link{olmm}} object.}
  \item{level}{character string. Indicates whether the returned values
    should on the observation level (\code{level = "observation"}) or on 
    the subject level (\code{level = "subject"}).}
  \item{norm}{logical. Whether residuals should be divided by their
    standard deviation.}
  \item{which}{optional character string. For \code{\link{fixef}}, it
    indicates whether \code{"all"} coefficients, the category-specific
    coefficients \code{"ce"} (i.e. non-proportional odds) or the
    global coefficients \code{"po"} (i.e. proportional odds) should be
    returned. For \code{model.matrix} it indicates whether the model
    matrix of the fixed- (\code{"fe"}) or the random effects
    (\code{"re"}) should be extracted.}
  \item{...}{potential further arguments passed to methods.}
}

\section{Methods}{
  \describe{
    \item{\code{anova(object, ...)}:}{compute analysis of
      deviance table for one or more fitted \code{olmm} objects.} 
    \item{\code{coef(object, ...)}:}{ extracts model coefficients.} 
    \item{\code{deviance(object, ...)}:}{returns the deviance of a
      fitted model.} 
    \item{\code{estfun.olmm(x, level = c("observation", "subject"),
	decorrelate = FALSE, center = TRUE, silent = FALSE, ...)}:}{extracts the 
      empirical estimating functions (negative partial derivatives of the
      log-likelihood of the fitted model. See \code{\link{gefp.olmm}}.} 
    \item{\code{extractAIC(fit, scale, k, ...)}:}{computes the AIC. See
      \code{\link{extractAIC}}.}  
    \item{\code{fitted(object, type, ...)}:}{extracts the fitted values
      from the model.}
    \item{\code{fixef(object, which = c("all", "npo", "po"),
	...)}:}{extracts estimated fixed effects.}
    \item{\code{formula(x, ...)}:}{extracts the formula of the object.}
    \item{\code{gefp.olmm(object, scores = estfun.olmm(object,
	decorrelate = TRUE), order.by = NULL, terms = NULL, subset =
	NULL, center = TRUE, ...)}:}{extracts the empirical cumulative
      score process of a model. See 
      \code{\link{gefp.olmm}}.}    
    \item{\code{getCall(x, ...)}:}{extract the call of the fitted
      model.} 
    \item{\code{logLik(object, ...)}:}{extract the log-likelihood of the
      fitted object.}
    \item{\code{model.frame(formula, ...)}:}{extract the model frame
      from the fitted object.}
    \item{\code{model.matrix(object, which, ...)}:}{extract the model
      matrix of the fitted model.}
    \item{\code{predict(object, newdata, type, ranef, na.action,
	...)}}{prediction for fitted values or new data. See
      \code{\link{predict.olmm}}.}   
    \item{\code{print(x, digits, ...)}:}{prints a short model
      summary. See \code{\link{summary.olmm}}.} 
    \item{\code{ranef(x, norm, ...)}}{extract the fitted random
      effects.} 
    \item{\code{ranefCov(object, ...)}:}{extract the variance-covariance
      matrix of the random effects from the fitted object.}
    \item{\code{resid(object, norm, ...)}:}{extract the residuals propagated
      by Li and Sheperd (2012). By default, the residuals from the
      marginal outcome distribution are computed. The conditional
      residuals can be computed by assigning \code{ranef = TRUE}.}
    \item{\code{reweight(object, weights, ...)}:}{re-fit the model with
      new weights.}
    \item{\code{show(object)}:}{show the summary of the fitted model.}
    \item{\code{summary(object, silent, ...)}:}{extracts a summary of
      the fitted model. See \code{\link{summary.olmm}}.}
    \item{\code{terms(x, which, ...)}:}{extracts the terms of model
      matrices. Available are \code{"fixef-po"},
      \code{"fixef-npo"}, \code{"ranef-po"} and
      \code{"ranef-npo"}} 
    \item{\code{update(object, formula., evaluate, ...)}:}{update and
      re-fit the object.}
    \item{\code{VarCorr(x, sigma, rdig, ...)}:}{extract the estimated
      variances, standard deviations and correlations for the random
      effects.}  
    \item{\code{vcov(object, ...)}:}{computes the variance-covariance
      matrix of the main parameters of the fitted model.}
    \item{\code{weights(object, level, ...)}:}{extract model weights.}
  }
}
\references{
  Agresti, A. (2010). \emph{Analysis of Ordinal Categorical Data}, 10
  edn, Wiley.
  
  Hartzel, J., Agresti A. and Caffo, B. (2001). Multinomial Logit Random
  Effect Models, \emph{Statistical Modelling} \bold{1}: 81--102
  
  Hedeker, D. and Gibbons, R. (1994). A random-effects ordinal
  regression model for multilevel analysis. \emph{Biometrics} \bold{20}
  (4): 933--944
  
  Hedeker, D. and Gibbons, R. (2006). Longitudinal Data Analysis. Wiley,
  Palo Alto, CA.

  Li, C. and Sheperd, B. E. (2012). A new residual for ordinal
  outcomes, \emph{Biometrika} \bold{99} (2): 437-480
  
  Tutz, G. and Hennevogl W. (1996). Random effects in ordinal regression
  models, \emph{Computational Statistics & Data Analysis} \bold{22} (5):
  537--557

  Vermunt, J.K and Hagenaars, J.A. (2004). Ordinal longitudinal data
  analysis, \emph{Methods in Human Growth Research}, Cambridge
  University Press.
}

\author{Reto Buergin}

\seealso{\code{\link{olmm}}, \code{\link{predict.olmm}},
  \code{\link{gefp.olmm}}} 

\examples{
showClass("olmm")
}

\keyword{classes}
\keyword{methods}