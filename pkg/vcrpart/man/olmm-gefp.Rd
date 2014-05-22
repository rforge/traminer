\name{olmm-gefp}

\alias{gefp.olmm}
\alias{predecor_control}
\alias{estfun.olmm}

\title{Methods for empirical processes of \code{\link{olmm}} objects}

\description{Methods to extract and pre-decorrelate the estimating
  equations (the negative maximum likelihood scores) and compute the
  empirical fluctuation process (the decorrelated, cumulative score
  process) from a fitted \code{\link{olmm}} object.
} 

\usage{
estfun.olmm(x, predecor = FALSE, control = predecor_control(),
            nuisance = NULL, ...)

predecor_control(impute = TRUE, seed = NULL, nit = 1L,
                 symmetric = TRUE, drop = TRUE,
                 reltol = 1e-6, maxit = 100L,
                 verbose = FALSE, silent = FALSE)

gefp.olmm(object, scores = NULL, order.by = NULL, subset = NULL,
          predecor = TRUE, parm = NULL, center = TRUE, drop = TRUE,
          silent = FALSE, ...)
}

\arguments{
  \item{x, object}{a fitted \code{\link{olmm}} object.}
  \item{predecor}{logical scalar. Indicates whether the within-subject
    correlation of the estimating equations should be removed by a linear
    transformation. See details.}
  \item{control}{a list of control parameter as produced by
    \code{\link{predecor_control}}.}
  \item{nuisance}{integer vector. Defines the coefficients which are
    regarded as nuisance and therefore omitted from the transformation.}
  \item{impute}{logical scalar. Whether missing values should be
    replaced using imputation.}
  \item{seed}{an integer scalar. Specifies the random number used for
    the \code{set.seed} call before the imputation. If set to
    \code{NULL}, \code{set.seed} is not processed.}
  \item{nit}{integer scalar. Specifies the number of score matrices to
    be imputed. Higher numbers will reduce the randomness of the
    imputations.}
  \item{symmetric}{logical scalar. Whether the transformation matrix
    should be symmetric.} 
  \item{drop}{logical. Whether singularities should be handled
    automatically (otherwise singularities yield an error).}
  \item{reltol}{convergence tolerance used to compute the transformation
    matrix.} 
  \item{maxit}{the maximum number of iterations used to compute the
    transformation matrix.}
  \item{silent}{logical scalar. Should the report of warnings be suppressed?}
  \item{verbose}{logical scalar. Produces messages.}
  \item{scores}{a function or a matrix. Function to extract the
    estimating equations from \code{object}, e.g.
    \code{\link{estfun.olmm}}, or a matrix representing the estimating
    equations.} 
  \item{order.by}{a numeric or factor vector. The explanatory variable to 
    be used to order the entries in the estimating equations. If set to 
    \code{NULL} (the default) the observations are assumed to be ordered.}
  \item{subset}{logical vector. For \code{\link{gefp.olmm}}, this argument
    extracts the subset of the estimating equations to be used. For 
    \code{estfun.olmm}, the argument specifies which observations of an 
    individual with more than \code{Nbal} observations should be kept / dropped. 
    See examples.}
  \item{parm}{integer, logical or a character vector. Extracts the
    columns of the estimating equations.}
  \item{center}{logical scalar. \code{TRUE} subtracts, if necessary, the
    column means of the estimating equations.} 
  \item{...}{arguments passed to other
    functions. \code{\link{gefp.olmm}} passes these arguments to
    \code{scores} if \code{scores} is a function.}
}

\value{
  \code{\link{predecor_control}} returns a list of control parameters
  for computing the pre-decorrelation transformation
  matrix. \code{\link{estfun.olmm}} returns a \code{\link{matrix}} with
  the estimating equations and \code{\link{gefp.olmm}} a list of class class
  \code{"gefp"}.  
}

\details{
  Complements the \code{estfun} method of the package \code{sandwich}
  and the \code{gefp} method of the package \code{strucchange} for
  \code{\link{olmm}} objects. \code{\link{estfun.olmm}} allows to
  pre-decorrelate the intra-individual correlation of observation
  scores, see argument \code{predecor}. The value returned by
  \code{gefp.olmm} may be used for testing coefficient constancy
  regarding an explanatory variable \code{order.by} by the \code{sctest}
  function of package \code{strucchange}, see the examples below. 

  If \code{predecor = TRUE} in \code{estfun.olmm}, a linear
  within-subject transformation is applied that removes (approximately) 
  the intra-subject correlation from the scores. Specifically,
  \eqn{u_{it}}, the ML score of the \eqn{t}'th observation of subject
  \eqn{i}, is transformed to \eqn{u^*_{it} = u_{it} + T \sum_{t'=1, t'
  \neq t}^{T_i} u_{it'}} such that \eqn{\mathrm{Cov}(u*_{it},
  u*_{it'}) =  \mathrm{Cov}(u*_{it}, u*_{i't''})} for all \eqn{i \neq
  i'} and \eqn{t \neq t'}.
  
  The pre-decorrelation approach above is prinicipally limited to
  balanced data. If the data of \code{object} are not balanced, a
  multiple immputation algorithm is applied to estimate values the
  scores would take if the data were balanced. Specifically, in each
  iteration of the imputation, (i) a set of predictors is drawn from the 
  empirical distribution of the model matrices in \code{object}, (ii) 
  responses are drawn based on predictors of (i) and \code{object}
  (cf. \code{\link{simulate.olmm}}) and (iii) the predictors of (i) and
  the responses of (ii) are added \code{object} and the scores are
  recomputed. The obtained sets of scores are averaged and then used to
  for transforming the scores instead of the original scores. The
  imputation algorithm can be controlled with \code{maxit.impute} (the
  number of imputation iterations), \code{abstol.impute} the maximum
  permitted quantile of differences between scores of the previous and
  the current step and \code{ratio.impute} specifies the quantile to be
  taken.
  
  Given a score matrix produced by \code{\link{estfun.olmm}}, the
  empirical fluctuation process can be computed by
  \code{\link{gefp.olmm}}. See Zeileis and Hornik
  (2007). \code{\link{gefp.olmm}} provides with \code{subset} and
  \code{parm} arguments specifically designed for nodewise tests in the
  \code{\link{tvcm}} algorithm. Using \code{subset} extracts the partial
  fluctuation process of the selected subset. \code{center = TRUE} makes 
  sure that the partial fluctuation process ends with zero. 
}

\references{
Zeileis A., Hornik K. (2007), Generalized M-Fluctuation Tests for Parameter
Instability, \emph{Statistica Neerlandica}, \bold{61}, 488--508.
doi:10.1111/j.1467-9574.2007.00371.x.
}

\author{Reto Buergin}

\seealso{\code{\link{olmm}}, \code{\link{olmm-class}}}

\examples{
## ------------------------------------------------------------------- #
## Dummy example 1:
##
## Testing coefficient constancy on 'z4' of the 'vcrpart_1' data.
## ------------------------------------------------------------------- #

data(vcrpart_1)

## extract a unbalanced subset to show to the full functionality of estfun
vcrpart_1 <- vcrpart_1[-seq(1, 100, 4),]
subset <- vcrpart_1$wave != 1L ## obs. to keep for fluctuation tests
table(table(vcrpart_1$id))

## fit the model
model <- olmm(y ~ treat + re(1|id), data = vcrpart_1)

## extract and pre-decorrelate the scores
scores <- estfun.olmm(model, predecor = TRUE,
                      control = predecor_control(verbose = TRUE,
                                                 nit = 5))
attr(scores, "T") # transformation matrix

## compute the empirical fluctuation process
fp <- gefp.olmm(model, scores, order.by = vcrpart_1$z4)

## process a fluctuation test
library(strucchange)
sctest(fp, functional = catL2BB(fp))
}

\keyword{methods}