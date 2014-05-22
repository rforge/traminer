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
\alias{neglogLik.tvcm}
\alias{nobs.tvcm}
\alias{predict.tvcm}
\alias{print.tvcm}
\alias{prune}
\alias{prune.tvcm}
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
              "sctest", "p.value",
              "riskgrid", "selected",
              "coef", "sd", "var"),
        step = NULL, ...)

\method{formula}{tvcm}(x, which = c("original", "root", "tree"), ...)

\method{neglogLik}{tvcm}(object, ...)

\method{predict}{tvcm}(object, newdata = NULL,
        type = c("link", "response", "prob", "class",
          "node", "coef", "ranef"),
        ranef = FALSE, na.action = na.pass, ...)

\method{prune}{tvcm}(tree, alpha = NULL, maxdepth = NULL,
      maxwidth = NULL, minsplit = NULL, minbucket = NULL,
      nselect = NULL, maxstep = NULL, terminal = NULL, ...)

\method{splitpath}{tvcm}(tree, step = 1L:tree$info$nstep, ...)

\method{summary}{tvcm}(object, ...)

}

\arguments{
  \item{object, x, tree}{an object of class \code{\link{tvcm}}.}
  \item{step}{integer vector. The iteration step from which information
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
  \item{which}{a character string specifying which quantities to extract.}
  \item{alpha}{numeric significance level. A node is splitted when the
    \eqn{p} value for any parameter stability test in that node falls
    below \code{alpha}.}
  \item{maxdepth}{integer. The maximum depth of the tree structure.}
  \item{maxwidth}{integer. The maximum width of the tree structure.}
  \item{minsplit}{numeric. The minimum sum of weights in a node in order
    to be considered for splitting.} 
  \item{minbucket}{numeric. The minimum sum of weights in a terminal
    node.}
  \item{nselect}{integer. The maximum number of paritioning variables in
    the tree structure.}
  \item{maxstep}{integer. The maximum number of steps of the algorithm.}
  \item{terminal}{integer vector. The inner nodes to be set to terminal nodes.}
  \item{...}{Additional arguments passed to the calls.}
}

\details{
  The \code{\link{predict}} function has two additional options for the
  \code{type} argument. The option \code{"node"} calls the node id and
  \code{"coef"} predicts the coefficients corresponding to an
  observation. 

  The \code{nselect} argument for \code{\link{prune}} reduces the tree
  structure to include \code{nselect} partitioning variables. If
  \code{nselect} is larger than the number of selected variables, the
  model remains the same. 

  \code{\link{splitpath}} is a new method to trace the fitting
  procedure. It shows several information, such as coefficient constancy
  tests or the risk of different splits, computed at each step.
  
  Further undocumented, available methods are: \code{\link{fitted}}, 
  \code{\link{getCall}}, \code{\link{logLik}},
  \code{\link{model.frame}}, \code{\link{nobs}}, \code{\link{print}},
  \code{\link{ranef}}, \code{\link{resid}}, and
  \code{\link{weights}}. All these methods have the same arguments as the
  corresponding default methods. 
}

\seealso{
  \code{\link{tvcm}}, \code{\link{tvcm-assessment}},
  \code{\link{tvcm-plot}} \code{\link{olmm-class}} 
}

\author{Reto Buergin}

\examples{
## ------------------------------------------------------------------- #
## Dummy example 1:
##
## Apply various methods on a fitted 'tvcm' object.
## ------------------------------------------------------------------- #

data(vcrpart_1)

control <- tvcm_control(minbucket = 10, verbose = TRUE)

model <- tvcolmm(y ~ wave + vc(z1, z2, z3, z4, by = treat) + re(1|id), 
                 data = vcrpart_1, subset = 1:150, control = control)

coef(model)
extract(model, "selected")
extract(model, "model")
predict(model, newdata = vcrpart_1[191:200,], type = "node")
predict(model, newdata = vcrpart_1[191:200,], type = "response")
prune(model, maxwidth = 1)
splitpath(model, 1:2)
summary(model, digits = 2)
}
\keyword{methods}