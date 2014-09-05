\name{olmm-predict}

\alias{olmm-predict}
\alias{fitted.olmm}
\alias{predict.olmm}

\title{Predict outcome probabilities and responses for \code{olmm}
  objects} 

\description{\code{fitted} and \code{predict} method for \code{olmm} objects.}

\usage{
\method{fitted}{olmm}(object, ...)

\method{predict}{olmm}(object, newdata = NULL,
        type = c("link", "response", "prob", "class", "ranef"),
        ranef = FALSE, na.action = na.pass, ...)
}

\arguments{
  \item{object}{a fitted \code{olmm} object.}
  \item{newdata}{data frame for which to evaluate predictions.}
  \item{type}{character string. \code{type = "response"} and 
    \code{type = "prob"} yield response probabilities, 
    \code{type = "class"} the response category with highest 
    probability and \code{type = "link"} the linear predictor
    matrix. \code{type = "ranef"} yields the predicted random effects,
    see \code{\link{ranef.olmm}}.}
  \item{ranef}{logical or numeric matrix. For all of
    \code{type = "response"}, \code{type = "prob"} or
    \code{type = "class"}, the option \code{ranef = FALSE} yields the
    marginal outcome probabilities resp. outcomes with highest marginal
    probability (random effects are (numerically) integrated out). If
    \code{type = "link"}, \code{ranef = FALSE} yields the linear
    predictor with random effects set to zero. \code{ranef = TRUE} yields
    conditional prediction based on subject-specific effects. If
    \code{newdata = NULL}, the predicted random effects for the learning
    data are used, otherwise a random effects matrix must be
    assigned. This matrix must contain a row for each subject and a
    column for each random effect. The rownames of this matrix are used to
    link the \code{subject} vector in \code{newdata} with the random
    effects. The argument is ignored if the model has no random
    effects. \code{type = "ranef"} predict the random effects and yields
    the same result as \code{ranef(object)}.} 
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

\seealso{\code{\link{olmm}}, \code{\link{olmm-methods}}}

\examples{
## ------------------------------------------------------------------- #
## Example 1: Schizophrenia
## ------------------------------------------------------------------- #

data(schizo)

model <- olmm(imps79o ~ tx + sqrt(week) + tx * sqrt(week) + re(1|id), schizo, subset = 1:181)

head(fitted(model))
head(predict(model))

## now make predictions for subject "1103" which is in the learning sample
newdata <- data.frame(id = c("1103", "1103"), tx = c(0, 1), week = c(3, 3))
ranef <- predict(model, type = "ranef")["1103",,drop=FALSE]

## marginal prediction
predict(model, newdata = newdata, type = "response", ranef = FALSE)

## conditional prediction
predict(model, newdata = newdata, type = "response", ranef = ranef)

## conditional prediction with ranef = 0
predict(model, newdata = newdata, type = "response", ranef = ranef * 0)

## predict the response with highest probability
predict(model, newdata = newdata, type = "class")

}
\keyword{methods}