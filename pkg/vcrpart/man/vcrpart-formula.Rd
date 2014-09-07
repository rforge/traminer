\name{vcrpart-formula}

\alias{vcrpart-formula}
\alias{fe}
\alias{vc}
\alias{re}
\alias{ce}
\alias{ge}

\title{Special terms for formulas.}

\description{Special terms for formulas assigned to \code{\link{tvcm}},
  \code{\link{fvcm}} and \code{\link{olmm}}.   
}

\usage{
fe(formula, intercept = TRUE)
re(formula, intercept = TRUE)
vc(..., by, intercept = missing(by), nuisance = character())
ce(formula)
ce(formula)
}

\arguments{
  \item{formula}{a symbolic description for the corresponding component
    of the formula component. See examples.}
  \item{intercept}{logical or character vector. \code{intercept = TRUE}
    (default) indicates that an intercept is incorporated.
    \code{intercept = FALSE} replaces the otherwise allowed \code{"-1"}
    term, that is ignored by \code{\link{fvcm}}, \code{\link{olmm}} and
    \code{\link{tvcm}}. Character strings may be used in connection with
    \code{\link{olmm}}. Intercepts have specific interpretations, see
    details.}
  \item{\ldots}{the names of moderators i.e. partitioning variables,
    separated by commas. It is also possibly to assign a character
    vector that includes all the variable names.} 
  \item{by}{a symbolic description for predictors the effect of which is
    moderated by the variables in \code{\ldots}. See \code{\link{tvcm}}
    and the examples below. Note that the by variable must be numeric
    and factor variables must be recoded to dummy variables by hand.}
  \item{nuisance}{character vector of variables in \code{by} which have
    to be estimated separately for each partition but the algorithm
    should not focus on this variable when searching for splits.}
}

\value{
  a list used by \code{\link{tvcm}}, \code{\link{fvcm}} and
  \code{\link{olmm}} for constructing the model formulas.
}

\details{
  Special formula terms to define fixed effects \code{\link{fe}}, 
  varying coefficients \code{\link{vc}} and random effects \code{\link{re}}. 
  The use of these formula terms ensures that \code{\link{fvcm}}, 
  \code{\link{tvcm}} and \code{\link{olmm}} fit the intended model. 
  Some examples are given below and in the corresponding documentation 
  pages. 

  Variables which are not defined within one of these three special
  terms will be assigned to the fixed effect predictor equation. The
  deletion of the intercept can be indicated by a \code{-1} or 
  \code{vc(intercept = FALSE)}. The terms \code{\link{ce}} (category-specific
  effects) and \code{\link{ge}} (global effect or proportional odds
  effect) are mainly designed for \code{\link{olmm}}. Notice that 
  \code{\link{tvcm}} may changes, for internal reasons, the order of the
  terms in the specified formula. At present, the term \code{"."},
  which is generally use to extract all variables of the data, is
  ignored. On the other hand, \code{\link{vc}} interprets character
  vectors, assigned as unnamed arguments, as lists of variables of 
  moderators to be extracted from \code{data}.
}

\seealso{\code{\link{tvcm}}, \code{\link{fvcm}}, \code{\link{olmm}}}

\examples{
## Formula for a model with 2 fixed effects (x1 and x2) and a random
## intercept.

formula <- y ~ fe(x1 + x2) + re(1|id)

## Formula for a model with 1 fixed effect and a varying coefficient term
## with 2 moderators and 2 varying coefficient predictors. 'tvcm' will
## fit one common partition for the two moderated predictors 'x2' and
## 'x3'.

formula <- y ~ x1 + vc(z1, z1, by = x2 + x3, intercept = TRUE)

## Similar formula as above, but the predictors 'x2' and 'x3' have
## separate 'vc' terms. 'tvcm' will fit a separate partition for each 
## 'vc' term

formula <- y ~ x1 + vc(z1, z1, by = x2 + x3, intercept = TRUE)
}

\author{Reto Buergin}