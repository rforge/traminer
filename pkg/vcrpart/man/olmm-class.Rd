\name{olmm-class}
\docType{class}
\Rdversion{1.1}
\alias{olmm-class}
\alias{show,olmm-method}

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
    \item{\code{score_obs}:}{a matrix of observation-wise partial
      derivates of the marginal log-likelihood equation.} 
    \item{\code{score_sbj}:}{a matrix of subject-wise partial derivates
      of the marginal log-likelihood equation.} 
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
    \item{\code{control}:}{a list of used control arguments produced by
      \code{\link{olmm_control}}.} 
    \item{\code{output}:}{the output of the optimizer (class
      \code{"list"}).} 
  }
}

\author{Reto Buergin}

\seealso{\code{\link{olmm}}, \code{\link{olmm-methods}},
  \code{\link{predict.olmm}}, \code{\link{gefp.olmm}}
}

\examples{
showClass("olmm")
}

\keyword{classes}
\keyword{classes}