\name{gefp.olmm}

\alias{estfun.olmm}
\alias{gefp.olmm}
\alias{olmm_scoreTransfMat}

\title{Methods for empirical processes of \code{\link{olmm}} objects}

\description{Methods to extract the estimating equations i.e. the
  negative maximum likelihood scores and the empirical fluctuation
  process i.e. the coefficient-wise decorrelated, cumulative score
  process from a fitted \code{\link{olmm}} object.  
} 

\usage{

estfun.olmm(x, level = c("observation", "subject"),
            prewhite = FALSE, prewhite.fail = c("stop", "ignore"),
            complete = FALSE, nuisance = NULL, control = list(),
            verbose = FALSE, silent = FALSE, ...)

gefp.olmm(object, scores = NULL,
          order.by = NULL, terms = NULL, subset = NULL,
          center = TRUE, omit.terms = TRUE,
          silent = FALSE, ...)

olmm_scoreTransfMat(object, terms = NULL, 
                    method = c("symmetric", "unconstraint"),
		    Nmax = NULL, control = list(),
                    verbose = FALSE, omit.terms = TRUE,
                    silent = FALSE)
}

\arguments{
  \item{object, x}{a fitted \code{olmm} object.}
  \item{level}{character string. Indicates whether the returned values
    should on the observation level (\code{level = "observation"}) or on 
    the subject level (\code{level = "subject"}).}
  \item{prewhite}{logical scalar. Indicates whether the within-subject
    correlation of scores should be removed by a linear
    transformation. See details.}
  \item{prewhite.fail}{character. Instructs whether an error should
    be reported (\code{"stop"}) or the raw estimating equations
    should be returned (\code{"ignore"}) in cases the prewhitening failed.}
  \item{complete}{logical scalar. If \code{TRUE} and the data are
    unbalanced, the estimation equations are adjusted. See details.}
  \item{nuisance}{integer vector. Defines the coefficients which are
    regarded as nuisance and therefore omitted from the transformation.}
  \item{control}{list. Control parameter for optimization algorithms. To 
    control the algorithm to handle unbalanced data, you can specify 
   \code{abstol} and \code{Rmax} and to control the algorithm computing
   the transformation matrix you can specify \code{reltol}, 
   \code{maxreltol}, \code{maxit} and \code{stopreltol}.}
  \item{verbose}{logical. Produces messages.}	 
  \item{silent}{logical scalar. Should the report of warnings be suppressed?}
  \item{scores}{function or matrix. Function to extract the estimating
    equations from a \code{\link{olmm}} object or a matrix representing
    the estimating equations. Commonly this is a matrix produced by
   \code{\link{estfun}}.}	       
  \item{center}{logical scalar. \code{TRUE} subtracts, if necessary, the
    column means of the estimating equations.} 
  \item{terms}{integer, logical or a character vector. Extracts the
    columns of the estimating equations.} 
  \item{subset}{subset of estimating function to be used.}
  \item{order.by}{a vector of a single explanatory variable to be used
    to order the estimating equations. If set to \code{NULL} (the default) the
    observations are assumed to be ordered.}
   \item{method}{character. \code{"symmetric"} (default) forces the transformation 
    matrix to be symmetric.}
  \item{Nmax}{integer. The number of observation per subject for which the 
    transformation matrix should be computed.}
  \item{omit.terms}{logical. Whether singularities should be handled
    automatically.}
  \item{...}{arguments passed to other functions. Specifically,
    \code{\link{estfun.olmm}} passes the arguments to
    \code{\link{olmm_scoreTransfMat}} and \code{\link{gefp.olmm}} to
    \code{\link{estfun.olmm}}. }
}

\value{
  \code{\link{olmm_scoreTransfMat}} returns a matrix for a linear 
  within-subject transformation of the estimating equations (i.e., the
  negative maximum likelihood scores), \code{\link{estfun.olmm}} a
  \code{\link{matrix}} with the estimating equations and
  \code{\link{gefp.olmm}} a list of class class \code{"gefp"}. 
}

\details{
  Complements the \code{estfun} method of the package \code{sandwich}
  and the \code{gefp} method of the package \code{strucchange} for
  \code{\link{olmm}} objects. The implementation is different from
  originals to allow the score function to be prewhitened with respect
  to intra-subject correlation. The value returned by \code{gefp.olmm}
  may be used for testing coefficient constancy regarding an explanatory
  variable \code{order.by} by the \code{sctest} function of package
  \code{strucchange}.

  If \code{prewhite = TRUE} in \code{estfun.olmm}, a linear
  within-subject transformation is applied that removes (approximately) 
  the intra-subject correlation from the scores. Specifically, u_it, 
  the ML score of observation t of subject i, is transformed to 
  u*_it = u_it + T sum(u_it'', t'' != t) such that Cov(u*_it, u*_it') = 
  Cov(u*_it, u*_i't'') for all i != i' and t != t'. The transformation 
  matrix T is computed with the \code{\link{olmm_scoreTransfMat}}
  function and involves a nonlinear optimization. According to our
  simulations, the method works stable for a model for a balanced
  dataset with about 100 subject and between 3 and 10 observations per
  individual. 

  If \code{complete = TRUE} and \code{prewhite = TRUE} and \code{level =
  "observation"}, an algorithm is applied to handle unbalanced data:
  First, the transformation matrix is estimated based on the scores of
  subjects that share the same maximal number of observations. Then, for
  the transformation of scores of subjects with fewer observations, we
  repeat the following steps: (i) simulate responses based on the model
  (including predicted random effects), (ii) recompute the score
  function (iii) transform the scores (iv) delete the simulated
  observations. The resulting scores are averaged.
  
  \code{gefp.olmm} proceeds in the following sequence:

  \enumerate{
    \item Extract the estimating equations (EE) using the \code{scores}
    argument (i.e. the negative ML scores).
    \item If \code{subset} is not \code{NULL}, extract the specified
    subset from the EE matrix and the \code{order.by} vector.
    \item If \code{center = TRUE}, subtract the column means from the
    EE.  
    \item Divide the EE by the number of observations.
    \item Multiply the EE by the inverse of the square root of their
    crossproduct. 
    \item Order the EE by the \code{order.by} vector and cumulate them.
    \item If \code{terms} is not \code{NULL}, extract the specified
    columns of the EE. 
  }

  The resulting estimating equations are stored in the \code{process}
  slot, see above.
  
  The \code{subset} and the \code{terms} arguments are  designed for the
  \code{\link{tvcolmm}} algorithm. Centering the scores is necessary if
  the column sums of the estimating equations is unequal 0, which may
  occurs if the scores are previously prewhitened.
}

\references{
Zeileis A., Hornik K. (2007), Generalized M-Fluctuation Tests for Parameter
Instability, \emph{Statistica Neerlandica}, \bold{61}, 488--508.
doi:10.1111/j.1467-9574.2007.00371.x.
}

\author{Reto Buergin}

\seealso{\code{\link{olmm}}, \code{\link{olmm-class}}}

\keyword{methods}