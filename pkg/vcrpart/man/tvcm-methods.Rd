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

\title{Methods for \command{\link{tvcm}} objects}

\description{Standard methods for computing on \command{\link{tvcm}}
  objects.} 

\usage{
\method{coef}{tvcm}(object, ...)

\method{extract}{tvcm}(object, what = c(
              "control", "model", 
              "nodes", "sctest", "p.value",
              "devgrid", "selected",
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
  \item{object, tree}{an object of class \command{\link{tvcm}}.}
  \item{steps}{integer vector. The iteration steps from which
    information should be extracted.} 
  \item{newdata}{an optional data frame in which to look for variables 
    with which to predict, if omitted, the fitted values are used.} 
  \item{type}{character string. Denotes for \command{\link{predict}} the
    type of predicted value. See \command{\link{predict.glm}} or
    \command{\link{predict.olmm}}.}
  \item{na.action}{function determining what should be done with missing
    values for fixed effects in \code{newdata}. The default is to
    predict \code{NA}: see \command{\link{na.pass}}.}
  \item{ranef}{logical scalar or matrix indicating whether prediction
    should be based on random effects. See
    \command{\link{predict.olmm}}. }
  \item{what}{a character specifying the quantities to \code{extract}.}
  \item{details}{logical scalar. Whether detail results like
    coefficient constancy tests or loss minimizing grid search should be 
    shown.}
  \item{...}{Additional arguments passed to the calls.}
}

\details{
  The \command{\link{predict}} function has two additional options for the
  \code{type} argument. The option \code{"node"} calls the node id and
  \code{"coef"} predicts the coefficients corresponding to an
  observation. In cases of multiple \command{\link{vc}} terms for the same
  predictor, the coefficients are summed up.

  \command{\link{splitpath}} and \command{\link{prunepath}} are new
  methods to trace the splitting resp. pruning procedures. They shows
  several information, such as the loss reduction of new splits during
  partitioning or the loss reduction of collapsing an inner node when
  pruning. 

  \command{\link{summary}} computes summary statistics of the fitted model, 
  including the estimated coefficients. The varying coefficient are printed
  by means of a printed decision tree. Notice that in cases there is no split
  for the varying coefficient, the average coefficient will be among the 
  fixed effects.	
  
  Further undocumented, available methods are: \command{\link{fitted}},
  \command{\link{formula}}, \command{\link{getCall}},
  \command{\link{logLik}}, \command{\link{model.frame}},
  \command{\link{nobs}}, \command{\link{print}}, \command{\link{ranef}},
  \command{\link{resid}}, and \command{\link{weights}}. All these
  methods have the same arguments as the corresponding default methods. 
}

\seealso{
  \command{\link{tvcm}}, \command{\link{tvcm-assessment}},
  \command{\link{tvcm-plot}}
}

\author{Reto Buergin}

\examples{
## ------------------------------------------------------------------- #
## Dummy example 1:
##
## Apply various methods on a 'tvcm' object fitted on the 'vcrpart_2'
## data. Cross-validation is omitted to accelerate the computations.
## ------------------------------------------------------------------- #

data(vcrpart_2)

model <- tvcm(y ~ -1 + vc(z1, z2) + vc(z1, z2, by = x1) + x2,
              data = vcrpart_2, family = gaussian(), subset = 1:90,
              control = tvcm_control(cv = FALSE))

coef(model)
extract(model, "selected")
extract(model, "model")
predict(model, newdata = vcrpart_2[91:100,], type = "node")
predict(model, newdata = vcrpart_2[91:100,], type = "response")
splitpath(model, steps = 1)
summary(model, digits = 2)
}
\keyword{methods}