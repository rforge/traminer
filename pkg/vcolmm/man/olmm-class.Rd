\name{olmm-class}
\Rdversion{1.1}
\docType{class}
\alias{olmm-class}
\alias{anova,olmm-method}
\alias{coefficients,olmm-method}
\alias{coef,olmm-method}
\alias{deviance,olmm-method}
\alias{drop1.olmm}
\alias{estfun.olmm}
\alias{fitted,olmm-method}
\alias{formula,olmm-method}
\alias{fixef,olmm-method}
\alias{getCall,olmm-method}
\alias{logLik,olmm-method}
\alias{model.frame,olmm-method}
\alias{model.matrix,olmm-method}
\alias{neglogLik,olmm-method}
\alias{predict,olmm-method}
\alias{print,olmm-method}
\alias{ranef,olmm-method}
\alias{ranefCov,olmm-method}
\alias{resid,olmm-method}
\alias{residuals,olmm-method}
\alias{reweight,olmm-method}
\alias{show,olmm-method}
\alias{summary,olmm-method}
\alias{terms,olmm-method}
\alias{terms.olmm}
\alias{update,olmm-method}
\alias{VarCorr,olmm-method}
\alias{vcov,olmm-method}
\alias{weights,olmm-method}

\title{Class \code{olmm}}

\description{
The \code{olmm} class represents 2-level ordinal linear mixed models
fitted by the \code{\link{olmm}} function.}

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
    \item{\code{formula}:}{the formula of the call (class
      \code{"formula"}).}
    \item{\code{terms}}{a list of \code{terms} of the fitted 
      model. The first element corresponds to \code{X}, the second to
      \code{W}.}
    \item{\code{y}:}{ordered response vector.} 
    \item{\code{X}:}{model matrix for the fixed effects.}
    \item{\code{W}:}{model matrix for the random effects.}
    \item{\code{subject}:}{a factor vector with grouping levels.}
    \item{\code{subjectName}:}{variable name of the subject vector.}
    \item{\code{weights}:}{numeric weights vector of length n.}
    \item{\code{weights_sbj}:}{numeric weights vector of length N.}
    \item{\code{offset}:}{numeric offset vector}
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
      the fixed effects (\code{object@fixef} times \code{object@X},
      random effects are not considered).}  
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
## Methods with "surprising" arguments

\S4method{fixef}{olmm}(object, which = c("all", "npo", "po"), ...)

\S4method{predict}{olmm}(object, newdata = NULL,
        type = c("prob", "class", "link"),
        na.action = na.pass, offset = NULL,
        ranef = FALSE, subset = NULL, ...)

\S4method{model.matrix}{olmm}(object, which = c("fixef", "ranef"), ...)

\S4method{weights}{olmm}(object, level = c("observation", "subject"),
...)

## S3 methods

\method{drop1}{olmm}(object, scope, scale = 0, test = c("none", "Chisq"),
       k = 2, trace = FALSE, ...)

\method{estfun}{olmm}(x, level = c("observation", "subject"), ...)

\method{terms}{olmm}(x, ...)
}

\arguments{
  \item{object, x}{an \code{\link{olmm}} object.}
  \item{k}{see \code{\link{drop1}}.}
  \item{level}{character string. Indicates whether the returned weights
    should be by observation or by subject.}
  \item{na.action, offset, subset}{further arguments to be passed to
    \code{\link{model.frame}} to extract the data.}
  \item{newdata}{an optional data frame in which to look for variables
    to predict. If omitted, the fitted values are used.}
  \item{ranef}{logical or a numeric matrix. The logical indicates
    whether the prediction is to be unconditional (FALSE) or conditional
    on random effects. If TRUE and \code{newdata = NULL} the predicted
    random effects will be used, otherwise a matrix with random effects
    is required.}  
  \item{scale}{see \code{\link{drop1}}.}
  \item{scope}{see \code{\link{drop1}}.}
  \item{test}{see \code{\link{drop1}}.}
  \item{trace}{see \code{\link{drop1}}.}
  \item{type}{character string. Denotes for \code{\link{predict}} the
    type of predicted value and for \code{\link{resid}} the type of
    residuals.}
  \item{which}{optional character string. For \code{\link{predict}} it
    indicates whether \code{all} coefficients, just the
    predictor-variable (\code{npo}) or the predictor-invariant (\code{po},
    i.e. proportional odds effects).. For \code{model.matrix} it indicates
    whether the model matrix of the fixed- or the random effects should
    be extracted.} 
  \item{...}{potential further arguments passed to methods.}
}

\section{Methods}{
  \describe{
    \item{\code{anova(object, ...)}:}{compute analysis of
      deviance table for one or more fitted \code{olmm} objects.} 
    \item{\code{coef(object, ...)}:}{ extracts model coefficients.} 
    \item{\code{deviance(object, ...)}:}{returns the deviance of a
      fitted model.}
    \item{\code{estfun(x, level, ...)}:}{extracts the
      empirical estimating functions i.e. negative scores i.e. negative
      partial derivatives of the log-likelihood of a fitted model.}
    \item{\code{fitted(object, type, ...)}:}{extracts the fitted values
      from the model.}
    \item{\code{fixef(object, which, ...)}:}{extract the fitted fixed
      effects.}
    \item{\code{formula(x, ...)}:}{extracts the formula of the object.} 
    \item{\code{getCall(x, ...)}:}{extract the call of the fitted
      model.} 
    \item{\code{logLik(object, ...)}:}{extract the log-likelihood of the
      fitted object.}
    \item{\code{model.frame(formula, ...)}:}{extract the model frame
      from the fitted object.}
    \item{\code{model.matrix(object, which, ...)}:}{extract the model
      matrix of the fitted model.}
    \item{\code{neglogLik(object, ...)}:}{extract negative
      log-likelihood.} 
    \item{\code{predict(object, newdata, type, na.action, offset,
	ranef, subset, ...)}:}{prediction for fitted values or new data.}
    \item{\code{print(x, digits, ...)}:}{print the model summary.}
    \item{\code{ranef(x, norm, ...)}}{extract the fitted random
      effects.} 
    \item{\code{ranefCov(object, ...)}:}{extract the variance-covariance
      matrix of the random effects from the fitted object.}
    \item{\code{resid(object, ...)}:}{extract the residuals propagated
      by Li and Sheperd (2012). By default, the residuals from the
      marginal outcome distribution are computed. The conditional
      residuals can be computed by assigning \code{ranef = TRUE}.}
    \item{\code{reweight(object, weights, ...)}:}{re-fit the model with
      new weights.}
    \item{\code{show(object)}:}{show the summary of the fitted model.}
    \item{\code{summary(object, ...)}:}{extract the summary of the
      fitted model.}
    \item{\code{terms(x, which, ...)}:}{extracts the terms of model
      matrices. Available are \code{"fixefEtaInv"},
      \code{"fixefEtaVar"}, \code{"ranefEtaInv"} and
      \code{"ranefEtaVar"}} 
    \item{\code{update(object, formula., evaluate, ...)}:}{update and
      re-fit the object.}
    \item{\code{VarCorr(x, ...)}:}{extract the estimated variances,
      standard deviations and correlations for the random effects.} 
    \item{\code{vcov(object, ...)}:}{calculate the variance-covariance
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

\seealso{\code{\link{olmm}}}

\examples{
showClass("olmm")
}

\keyword{classes}
\keyword{methods}