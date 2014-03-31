\name{predict.olmm}

\alias{predict.olmm}

\title{Predict outcome probabilities and responses for \code{olmm}
  objects} 

\description{\code{predict} method for \code{olmm} objects}

\usage{
\method{predict}{olmm}(object, newdata = NULL,
        type = c("prob", "class", "link"),
        ranef = FALSE, na.action = na.pass, ...)
}

\arguments{
  \item{object}{a fitted \code{olmm} object.}
  \item{newdata}{data frame for which to evaluate predictions.}
  \item{type}{character string. \code{type = "prob"} yields response
    probabilities, \code{type = "class"} the response category with
    highest probability and \code{type = "link"} the linear predictor
    matrix.}
  \item{ranef}{logical or a numeric matrix. For \code{type = "prob"} or
    \code{type = "class"}, \code{ranef = FALSE} yields the
    marginal outcome probabilities resp. outcomes with highest marginal
    probability. Marginal means that random effects are (numerically)
    integrated out. If \code{type = "link"}, \code{ranef = FALSE} yields
    the linear predictor with random effects set to zero. \code{ranef =
      TRUE} yields conditional prediction based on subject-specific
    effects and works only if \code{newdata = NULL}, i.e. for prediction
    on the learning data. Conditional prediction for new data can be
    accomplished by assigning a numeric random effects matrix. This
    requires the \code{newdata} to include a column indicating the
    subject and \code{ranef} to include a row for each subject in
    \code{newdata} with corresponding rowname labels.}
  \item{na.action}{function determining what should be done with missing
    values for fixed effects in \code{newdata}. The default is to
    predict \code{NA}: see \code{\link{na.pass}}.} 
  \item{...}{optional additional parameters. Includes \code{offset} and
    \code{subset}.} 
}

\value{
  A matrix or a vector of predicted values or response probabilities. 
}

\note{
  The method can not yet handle new categories in categorical predictors
  and will return an error. 
}

\author{Reto Buergin}

\seealso{\code{\link{olmm}}, \code{\link{olmm-class}}}

\keyword{methods}