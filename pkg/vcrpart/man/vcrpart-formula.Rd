\name{tvcm-formula}

\alias{fe}
\alias{re}
\alias{vc}
\alias{ce}
\alias{ge}

\title{Special functions for formulas.}

\description{Functions to be used for the formula argument of
  \code{\link{tvcm}}, \code{\link{fvcm}} and \code{\link{olmm}}. 
}

\usage{
fe(formula, intercept = TRUE)
re(formula, intercept = TRUE)
vc(..., by, intercept = TRUE)
ce(formula)
ce(formula)
}

\arguments{
  \item{formula}{a symbolic description of the model. See examples.}
  \item{...}{the names of moderators i.e. effect modifiers or
    partitioning variables, separated by a comma.}
  \item{by}{a predictor equation of variables the effect of which is
    modified.} 
  \item{intercept}{logical or character string. \code{TRUE} indicates
    that a intercept is incorporated. Character strings may be used for
    categorical regression models using \code{\link{olmm}}.}
}

\value{
  a list used by \code{\link{tvcm}}, \code{\link{fvcm}} and
  \code{\link{olmm}} for constructing the model formulas.
}

\details{
  \code{fe} is used to specify the fixed effects of a model,
  e.g. \code{y ~ fe(x1 + x2)}. \code{re} is used to specify random
  effects. 
}

\seealso{\code{\link{tvcm}}, \code{\link{fvcm}}, \code{\link{olmm}}}

\author{Reto Buergin}