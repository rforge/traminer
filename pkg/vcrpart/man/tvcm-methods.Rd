\name{tvcm-methods}

\alias{tvcm-methods}
\alias{coef.tvcm}
\alias{coefficients.tvcm}
\alias{extract}
\alias{extract.tvcm}
\alias{formula.tvcm}
\alias{fitted.tvcm}
\alias{getCall.tvcm}
\alias{logLik.tvcm}
\alias{model.frame.tvcm}
\alias{neglogLik2.tvcm}
\alias{nobs.tvcm}
\alias{predict.tvcm}
\alias{prunepath}
\alias{prunepath.tvcm}
\alias{print.tvcm}
\alias{ranef.tvcm}
\alias{resid.tvcm}
\alias{residuals.tvcm}
\alias{splitpath}
\alias{splitpath.tvcm}
\alias{summary.tvcm}
\alias{weights.tvcm}

\title{Methods for \code{tvcm} objects}

\description{Standard methods for computing on \code{tvcm} objects.}

\usage{
\method{coef}{tvcm}(object, ...)

\method{extract}{tvcm}(object, what = c(
              "control", "model", 
              "nodes", "sctest", "p.value",
              "lossgrid", "selected",
              "coef", "sd", "var"),
        steps = NULL, ...)

\method{neglogLik2}{tvcm}(object, ...)

\method{predict}{tvcm}(object, newdata = NULL,
        type = c("link", "response", "prob", "class",
          "node", "coef", "ranef"),
        ranef = FALSE, na.action = na.pass, ...)

\method{splitpath}{tvcm}(tree, steps = 1L,
         details = FALSE, ...)

\method{prunepath}{tvcm}(tree, steps = 1L, ...)

\method{summary}{tvcm}(object, ...)

}

\arguments{
  \item{object, tree}{an object of class \code{\link{tvcm}}.}
  \item{steps}{integer vector. The iteration steps from which information
    should be extracted.} 
  \item{newdata}{an optional data frame in which to look for variables
    with which to predict, if omitted, the fitted values are used.}
  \item{type}{character string. Denotes for \code{\link{predict}} the
    type of predicted value. See \code{\link{predict.glm}} or
    \code{\link{predict.olmm}}.}
  \item{na.action}{function determining what should be done with missing
    values for fixed effects in \code{newdata}. The default is to
    predict \code{NA}: see \code{\link{na.pass}}.}
  \item{ranef}{logical scalar or matrix indicating whether prediction
    should be based on random effects. See \code{\link{predict.olmm}}.}
  \item{what}{a character specifying the quantities to \code{extract}.}
  \item{details}{logical scalar. Whether detail results like
    coefficient constancy tests or loss minimizing grid search should be
    shown.}
  \item{...}{Additional arguments passed to the calls.}
}

\details{
  The \code{\link{predict}} function has two additional options for the
  \code{type} argument. The option \code{"node"} calls the node id and
  \code{"coef"} predicts the coefficients corresponding to an
  observation. 

  \code{\link{splitpath}} is a new method to trace the fitting
  procedure. It shows several information, such as coefficient constancy
  tests or the loss of different splits, computed at each step.
  
  Further undocumented, available methods are: \code{\link{fitted}}, 
  \code{\link{formula}}, \code{\link{getCall}}, \code{\link{logLik}},
  \code{\link{model.frame}}, \code{\link{nobs}}, \code{\link{print}},
  \code{\link{ranef}}, \code{\link{resid}}, and
  \code{\link{weights}}. All these methods have the same arguments as the
  corresponding default methods. 
}

\seealso{
  \code{\link{tvcm}}, \code{\link{tvcm-assessment}},
  \code{\link{tvcm-plot}}
}

\author{Reto Buergin}

\examples{
## ------------------------------------------------------------------- #
## Dummy example 1:
##
## Apply various methods on a fitted 'tvcm' object.
## ------------------------------------------------------------------- #

data(vcrpart_1)

control <- tvcm_control(sctest = TRUE, minsize = 10, verbose = TRUE)

model <- tvcolmm(y ~ -1 + wave + vc(z3, z4, by = treat, intercept = TRUE) + re(1|id), 
                 data = vcrpart_1, subset = 1:150, control = control)

coef(model)
extract(model, "selected")
extract(model, "model")
predict(model, newdata = vcrpart_1[191:200,], type = "node")
predict(model, newdata = vcrpart_1[191:200,], type = "response")
splitpath(model, 1)
summary(model, digits = 2)
}
\keyword{methods}