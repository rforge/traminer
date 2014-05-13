\name{vcrpart-formula}

\alias{vcrpart-formula}
\alias{fe}
\alias{re}
\alias{vc}
\alias{ce}
\alias{ge}

\title{Special terms for formulas.}

\description{Terms for the formula argument of \code{\link{tvcm}},
  \code{\link{fvcm}} and \code{\link{olmm}}.  
}

\usage{
fe(formula, intercept = TRUE)
re(formula, intercept = TRUE)
vc(..., by, intercept = TRUE)
ce(formula)
ce(formula)
}

\arguments{
  \item{formula}{a symbolic description for predictors in the model
    term. See examples.} 
  \item{...}{the names of moderators i.e. partitioning variables,
    separated by commas.} 
  \item{by}{a symbolic description for predictors the effect of which is
    moderated by the variables in \code{...}. See \code{\link{tvcm}} and
    the examples below.} 
  \item{intercept}{logical or character. \code{intercept = TRUE}
    indicates that an intercept is incorporated. For models based on
    \code{\link{olmm}}, \code{intercept = "ce"} incorporates category
    specific intercepts and \code{intercept = "ge"} incorporates global
    i.e. proportional odds intercepts. For \code{\link{vc}} terms, see
    the detail section in \code{\link{tvcm}}.}
}

\value{
  a list used by \code{\link{tvcm}}, \code{\link{fvcm}} and
  \code{\link{olmm}} for constructing the model formulas.
}

\details{
  The \pkg{vcrpart} package provides special formula terms for fixed effects 
  \code{\link{fe}}, random effects \code{\link{re}} and varying coefficients 
  \code{\link{vc}}. We recommend the use of these formula terms to ensure that
  the correct model is fitted. Some examples are given below.

  Variables which are not defined within one of these three special
  terms will be assigned to the fixed effect predictor equation. The
  deletion of intercepts must be defined using \code{intercept = "none"}
  in the corresponding special term. Terms like "\code{-1}" are
  ignored. The handling of intercept terms in connection with
  \code{\link{vc}} terms is discussed in detail in the details section
  of \code{\link{tvcm}}. The terms \code{\link{ce}} (category-specific
  effects) and \code{\link{ge}} (global effect) are mainly designed for
  categorical regression models, see \code{\link{olmm}}. Notice that the
  model model will not keep the exact order of the terms except for all
  \code{\link{ce}} terms appear in the formula before the
  \code{\link{ge}} terms. At the time being, the term "\code{.}", which
  is generally use to extract all variables in the data frame, is ignored.
}

\seealso{\code{\link{tvcm}}, \code{\link{fvcm}}, \code{\link{olmm}}}

\examples{
## Formula for a model with 2 fixed effects (x1 and x2) and a random
## intercept.

formula <- y ~ fe(x1 + x2) + re(1|id)

## Formula for a model with 1 fixed effect and a varying coefficient term
## with 2 moderators and 2 varying coefficient predictors.

formula <- y ~ x1 + vc(z1, z1, by = x2 + x3, intercept = TRUE)
}

\author{Reto Buergin}