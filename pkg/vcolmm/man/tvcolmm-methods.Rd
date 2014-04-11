\name{tvcolmm-methods}

\alias{tvcolmm-methods}
\alias{coef.tvcolmm}
\alias{coefficients.tvcolmm}
\alias{extract.tvcolmm}
\alias{fitted.tvcolmm}
\alias{getCall.tvcolmm}
\alias{logLik.tvcolmm}
\alias{model.frame.tvcolmm}
\alias{neglogLik.tvcolmm}
\alias{predict.tvcolmm}
\alias{print.tvcolmm}
\alias{prune.tvcolmm}
\alias{ranef.tvcolmm}
\alias{resid.tvcolmm}
\alias{residuals.tvcolmm}

\title{Methods for \code{tvcolmm} objects}
\description{Standard methods for computing on \code{tvcolmm} objects.}

\usage{
\method{coef}{tvcolmm}(object, ...)

\method{extract}{tvcolmm}(object, what = c("control", "sctest",
              "model", "selected", "p.value"), ids = nodeids(object), ...)

\method{fitted}{tvcolmm}(object, ...)

\method{getCall}{tvcolmm}(x, ...)

\method{logLik}{tvcolmm}(object, cv = FALSE,
               folds = tvcolmm_folds(object, "kfold"), ...)

\method{model.frame}{tvcolmm}(formula, ...)

\method{neglogLik}{tvcolmm}(object, ...)

\method{predict}{tvcolmm}(object, newdata = NULL,
        type = c("prob", "class", "node", "link", "terms"), ...)

\method{print}{tvcolmm}(x, ...)

\method{prune}{tvcolmm}(tree, alpha = NULL, depth = NULL,
            minsplit = NULL, nselect = NULL, ...)

\method{ranef}{tvcolmm}(object, ...)

\method{resid}{tvcolmm}(object, ...)

}
\arguments{
  \item{object, x, formula, tree}{an object of class \code{\link{tvcolmm}}.}
  \item{ids}{integer. A vector of node identifiers to apply over.}
  \item{newdata}{an optional data frame in which to look for variables
    with which to predict, if omitted, the fitted values are used.}
  \item{type}{character string. Denotes for \code{\link{predict}} the
    type of predicted value.}
  \item{what}{a character specifying the quantities to \code{extract}.}
  \item{alpha}{numeric significance level. A node is splitted when the
    \eqn{p} value for any parameter stability test in that node falls
    below \code{alpha}.}
  \item{depth}{integer. The maximum depth of the tree structure.}
  \item{minsplit}{numeric. The minimum number of observations in a
    node.} 
  \item{nselect}{integer. The maximum number of paritioning variables in
    the tree structure.}
  \item{folds}{folds for cross-validation. See
    \code{\link{tvcolmm_folds}}.}
  \item{cv}{logical scalar. Whether the log-likelihood should be
    computed via cross validation.}
  \item{...}{Additional arguments passed to the calls.}
}

\details{
  Except for \code{print}, \code{predict} and \code{prune}, these
  methods extract the underneath \code{olmm} object and assign that
  object to the equally named methods. The detailed documentation with
  descriptions of arguments is in \code{\link{olmm-class}}.
  
  The \code{predict} function has two additional options for the
  \code{type} argument. The option \code{"node"} calls the node id of an
  observation in the tree structure. The option \code{"terms"} predicts
  the coefficients corresponding to an observation.

  For pruning via \code{\link{prune}}, the \code{nvar} argument
  reduces \code{alpha} such that \code{nvar} the number of effect
  modifiers remain in the tree structure. If the number of effect
  modifiers in \code{object} is smaller than \code{nvar}, the tree
  structure remains the same.
}

\seealso{
  \code{\link{cv.tvcolmm}}, \code{\link{prune.tvcolmm}},
  \code{\link{tvcolmm}}, \code{\link{tvcolmm-plot}},
  \code{\link{olmm-class}} 
}

\author{Reto Buergin}

\keyword{methods}