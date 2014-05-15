\name{fvcm-methods}

\alias{fvcm-methods}
\alias{fitted.fvcm}
\alias{print.fvcm}
\alias{oobrisk.fvcm}
\alias{plot.fvcm}
\alias{predict.fvcm}
\alias{print.fvcm}
\alias{ranef.fvcm}

\title{Methods for \code{fvcm} objects}

\description{Standard methods for computing on \code{\link{fvcm}} objects.}

\usage{
\method{oobrisk}{fvcm}(object, fun = NULL, ranef = FALSE, ...)

\method{plot}{fvcm}(x, type = c("default", "coef", 
       "simple", "terms"),
     which = 1L, ask = TRUE, ...)

\method{predict}{fvcm}(object, newdata = NULL,
        type = c("link", "response", "prob", "class", "coef", "ranef"),
        ranef = FALSE, na.action = na.pass, verbose = FALSE, ...)
}

\arguments{
  \item{object, x}{an object of class \code{\link{fvcm}}.}
  \item{fun}{the risk function. The default risk function is defined
    as the sum of the deviance residuals. For a user defined
    function \code{fun}, see the examples of
    \code{\link{oobrisk.tvcm}}.} 
  \item{newdata}{an optional data frame in which to look for variables
    with which to predict. If omitted, the training data are used.}
  \item{type}{character string indicating the type of plot or
    prediction. See \code{\link{plot.tvcm}} or
    \code{\link{predict.tvcm}}.}
  \item{which}{integer vector. Which trees should be plotted.}
  \item{ask}{logical scalar. Whether a user input should be asked before
    plotting the next tree.}
  \item{ranef}{logical scalar or matrix indicating whether predictions
    should be based on random effects. See \code{\link{predict.olmm}}.}
  \item{na.action}{function determining what should be done with missing
    values for fixed effects in \code{newdata}. The default is to
    predict \code{NA}: see \code{\link{na.pass}}.}
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the validation.}
  \item{...}{further arguments passed to other methods.}
}

\details{
  \code{\link{oobrisk.fvcm}} estimates the out-of-bag risk based on
  predictions on the model that aggregates only those trees in which the
  observation didn't appear (cf. Hastie et al, 2001, sec. 15). The prediction 
  error is computed as the sum of prediction errors obtained with \code{fun}, 
  which are the deviance residuals by default.

  The plot and the prediction methods are analogous to handle as
  \code{\link{plot.tvcm}} resp. \code{\link{predict.tvcm}}. Notice that
  the plot options \code{mean} and \code{conf.int} for
  \code{type = "coef"} resp are not available (and internally set to
  \code{FALSE}). 

  Further undocumented, available methods are \code{\link{fitted}},
  \code{\link{print}} and \code{\link{ranef}}. All these latter methods 
  have the same arguments as the corresponding default methods.  
}

\author{Reto Buergin}

\references{
  T. Hastie, R. Tibshirani, J. Friedman (2001), The elements of
  statistical learning, Springer.
}

\seealso{\code{\link{fvcm}}, \code{\link{tvcm-methods}}}

\examples{
## ------------------------------------------------------------------- #
## Dummy example 1:
##
## Fitting a random forest tvcm on artificially generated ordinal
## longitudinal data.
## ------------------------------------------------------------------- #

## load the data

data(vcrpart_1)

## fit and analyse the model

control <- fvcm_control(mtry = 2, maxwidth = 3,
                        minbucket = 15, maxevalsplit = 2)
folds <- cvfolds(vcrpart_1$id, type = "subsampling", K = 2, prob = 0.75)

model.1 <- fvcolmm(y ~ wave + vc(z1, z2, z3, z4, by=treat) + re(1|id),
                   family = cumulative(), subset = 1:100,
                   data = vcrpart_1, control = control, folds = folds)

suppressWarnings(oobrisk(model.1))

plot(model.1, "coef", 1:2)

## predicting responses and varying coefficients for subject '27'
subs <- vcrpart_1$id == "27"

## predict coefficients
predict(model.1, newdata = vcrpart_1[subs,], type = "coef")

## marginal response prediction
predict(model.1, vcrpart_1[subs,], "response", ranef = FALSE)

## conditional response prediction
re <- matrix(5, 1, 1, dimnames = list("27", "(Intercept)"))
predict(model.1, vcrpart_1[subs,], "response", ranef = re)
predict(model.1, vcrpart_1[subs,], "response", ranef = 0 * re)

## predicting in-sample random effects
head(predict(model.1, type = "ranef"))

## fitted responses (marginal and conditional prediction)
head(predict(model.1, type = "response", ranef = FALSE))
head(predict(model.1, type = "response", ranef = TRUE))

## ------------------------------------------------------------------- #
## Dummy example 2:
##
## Fitting a random forest tvcm on artificially generated normally
## distributed data.
## ------------------------------------------------------------------- #

data(vcrpart_2)

## fit and analyse the model

control <- fvcm_control(mtry = 1L, minbucket = 0, maxwidth = 3L)
folds <- cvfolds(factor(1:nrow(vcrpart_2)), "subsampling", 3L, 0.75)

model.2 <- fvcglm(y ~ vc(z1, z2, by = x1) + x2, data = vcrpart_2,
                 family = gaussian(), subset = 1:50,
                 folds = folds, control = control)

suppressWarnings(oobrisk(model.2))

plot(model.2, "coef", 1:3)

predict(model.2, vcrpart_2[91:100, ], "coef")
}

\keyword{methods}
\keyword{hplot}