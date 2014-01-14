\name{olmm-summary}
\Rdversion{1.1}
\docType{class}
\alias{summary.olmm-class}
\alias{print,summary.olmm-method}
\alias{show,summary.olmm-method}

\title{Class \code{summary.olmm}}

\description{The \code{summary.olmm} class represents summaries of
  \code{\link{olmm}} objects.    
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.olmm",
  ...)}. 
}

\section{Slots}{
  \describe{
    \item{\code{methTitle}:}{character string. The title to be returned
      at \code{\link{show}} or \code{\link{print}} calls.}
    \item{\code{formula}:}{character string. The called formula.}
    \item{\code{data}:}{character string. The name of the data set
      used.} 
    \item{\code{subset}:}{character string. A description of the used
      subset of \code{data}.}
    \item{\code{AICtab}:}{a data frame of values of \code{\link{AIC}},
      \code{\link{BIC}}, log-likelihood and \code{\link{deviance}}.}
    \item{\code{FEmatEtaInv}:}{the matrix of estimates, standard errors,
      etc. for the predictor-invariant fixed-effects coefficients}
    \item{\code{FEmatEtaVar}:}{the matrix of estimates, standard errors,
      etc. for the predictor-variable fixed-effects coefficients}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{dims}:}{}
  }  
}

\section{Methods}{
  \describe{
    \item{\code{print(x, digits, ...)}:}{print information about the
      fitted model.} 
    \item{\code{show(object)}:}{same as the \code{\link{print}} method.}
  }
}

\seealso{\code{\link{olmm}}, \code{\link{olmm-class}}}

\author{Reto Buergin}

\keyword{classes}