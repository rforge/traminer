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

\method{fitted}{tvcm}(object, ...)

\method{formula}{tvcm}(x, which = c("original", "root", "tree"), ...)

\method{getCall}{tvcm}(x, ...)

\method{logLik}{tvcm}(object, ...)

\method{model.frame}{tvcm}(formula, ...)

\method{nobs}{tvcm}(object, ...)

\method{predict}{tvcm}(object, newdata = NULL,
        type = c("link", "response", "prob", "class",
          "node", "coef"), ...)

\method{print}{tvcm}(x, ...)

\method{prune}{tvcm}(tree, alpha = NULL, maxdepth = NULL,
     maxwidth = NULL, minsplit = NULL, minbucket = NULL, nselect = NULL,
     maxstep = NULL, terminal = NULL, ...)

\method{ranef}{tvcm}(object, ...)

\method{resid}{tvcm}(object, ...)

\method{splitpath}{tvcm}(tree, step = 1L:tree$info$nstep, ...)

\method{summary}{tvcm}(object, ...)

\method{weights}{tvcm}(object, ...)
}

\arguments{
  \item{object, x, formula, tree}{an object of class \code{\link{tvcm}}.}
  \item{step}{integer vector. The iteration step from which information
    should be extracted.} 
  \item{newdata}{an optional data frame in which to look for variables
    with which to predict, if omitted, the fitted values are used.}
  \item{type}{character string. Denotes for \code{\link{predict}} the
    type of predicted value.}
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
  The \code{predict} function has two additional options for the
  \code{type} argument. The option \code{"node"} calls the node id and
  \code{"coef"} predicts the coefficients corresponding to an
  observation. 

  The \code{nselect} argument for pruning reduces the tree to consist of
  the assigned number of variables. If \code{nselect} is larger than the
  number of selected variables, the model remains the same.
}

\seealso{
  \code{\link{cvrisk.tvcm}}, \code{\link{prune.tvcm}},
  \code{\link{tvcm}}, \code{\link{tvcm-plot}},
  \code{\link{olmm-class}} 
}

\author{Reto Buergin}

\keyword{methods}