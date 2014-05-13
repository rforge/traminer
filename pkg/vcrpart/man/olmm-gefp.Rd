\name{olmm-gefp}

\alias{estfun.olmm}
\alias{gefp.olmm}
\alias{decormat.olmm}

\title{Methods for empirical processes of \code{\link{olmm}} objects}

\description{Methods to extract and pre-decorrelate the estimating
  equations (the negative maximum likelihood scores) and compute the
  empirical fluctuation process (the decorrelated, cumulative score
  process) from a fitted \code{\link{olmm}} object.
} 

\usage{

estfun.olmm(x, level = c("observation", "subject"),
            predecor = FALSE, complete = predecor,
            Nbal = NULL, subset = NULL,
            nuisance = NULL, control = list(),
            verbose = FALSE, silent = FALSE, ...)

gefp.olmm(object, scores = NULL, predecor = TRUE,
          order.by = NULL, parm = NULL, subset = NULL,
          center = TRUE, drop = TRUE, silent = FALSE, ...)

decormat.olmm(object, method = c("symmetric", "unconstraint"),
              Nbal = NULL, parm = NULL, control = list(),
              verbose = FALSE, drop = TRUE, silent = FALSE)
}

\arguments{
  \item{object, x}{a fitted \code{olmm} object.}
  \item{level}{character string. Indicates whether the returned values
    should on the observation level (\code{level = "observation"}) or on 
    the subject level (\code{level = "subject"}).}
  \item{predecor}{logical scalar. Indicates whether the within-subject
    correlation of scores should be removed by a linear
    transformation. See details.}
  \item{complete}{logical scalar. If \code{TRUE} and the data are
    unbalanced, the estimation equations are adjusted. See details.}
  \item{Nbal}{scalar integer. The number of observations per subject for
    which the transformation should be performed when data are unbalanced.}
  \item{subset}{a vector to extract the subset of scores to be used. For
    \code{estfun.olmm}, this argument can be used to specify which
    observations of an individual with more than Nbal observations
    should be kept / dropped. See examples.}
  \item{nuisance}{integer vector. Defines the coefficients which are
    regarded as nuisance and therefore omitted from the transformation.}
  \item{control}{a list. Control parameter for optimization algorithms. To 
    control the algorithm to handle unbalanced data, you can specify 
   \code{abstol} and \code{Rmax} and to control the algorithm computing
   the transformation matrix you can specify \code{reltol}, 
   \code{maxreltol}, \code{maxit} and \code{stopreltol}.}
  \item{verbose}{logical scalar. Produces messages.}	 
  \item{silent}{logical scalar. Should the report of warnings be suppressed?}
  \item{scores}{function or matrix. Function to extract the estimating
    equations from a \code{\link{olmm}} object or a matrix representing
    the estimating equations. Commonly this is a matrix produced by
   \code{\link{estfun}}.}	       
  \item{center}{logical scalar. \code{TRUE} subtracts, if necessary, the
    column means of the estimating equations.} 
  \item{parm}{integer, logical or a character vector. Extracts the
    columns of the estimating equations.} 
  \item{order.by}{a vector of a single explanatory variable to be used
    to order the estimating equations. If set to \code{NULL} (the default) the
    observations are assumed to be ordered.}
   \item{method}{character. \code{"symmetric"} (default) forces the transformation 
    matrix to be symmetric.}
  \item{drop}{logical. Whether singularities should be handled
    automatically (otherwise singularities yield an error).}
  \item{...}{arguments passed to other functions. Specifically,
    \code{\link{estfun.olmm}} passes the arguments to
    \code{\link{decormat.olmm}} and \code{\link{gefp.olmm}} to
    \code{\link{estfun.olmm}}. }
}

\value{
  \code{\link{decormat.olmm}} returns a matrix for a linear 
  within-subject transformation of the estimating equations (i.e., the
  negative maximum likelihood scores), \code{\link{estfun.olmm}} a
  \code{\link{matrix}} with the estimating equations and
  \code{\link{gefp.olmm}} a list of class class \code{"gefp"}. 
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
  i'} and \eqn{t \neq t'}. The transformation matrix \eqn{T} is computed
  with the \code{decormat.olmm} function and involves a nonlinear
  optimization. 
  
  The pre-decorrelation approach above is prinicipally limited to
  balanced data. Setting \code{complete = TRUE}, an algorithm is applied
  to tackle this problem: First, the transformation matrix \eqn{T} is
  computed based on the scores of subjects with equal or more than
  \code{Nbal} observations. Second, the scores of subjects with fewer
  observations are computed by repeating the following steps: (i)
  simulate responses based \code{object} (including predicted random
  effects) to have \code{Nbal} responses for each of these subjects,
  (ii) recompute the score function (iii) transform the scores (iv)
  delete the simulated observations. The resulting scores are averaged. 
  
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
vcrpart_1 <- vcrpart_1[-c(seq(1, 100, 4), 2, 6),]
subset <- vcrpart_1$wave != 1L ## obs. to keep for fluctuation tests
table(table(vcrpart_1$id))

## fit the model
model <- olmm(y ~ treat + re(1|id), data = vcrpart_1)

## extract and pre-decorrelate the scores
scores <- estfun.olmm(model, predecor = TRUE, complete = TRUE,
                      Nbal = 3, subset = subset, verbose = TRUE)
attr(scores, "T") # transformation matrix

## compute the empirical fluctuation process
fp <- gefp.olmm(model, scores, order.by = vcrpart_1$z4)

## process a fluctuation test
library(strucchange)
sctest(fp, functional = catL2BB(fp))
}

\keyword{methods}