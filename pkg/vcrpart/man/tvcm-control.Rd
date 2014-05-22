\name{tvcm-control}

\alias{tvcm_control}

\title{Control parameters for \code{\link{tvcm}} and \code{\link{fvcm}}.} 

\description{
  Various parameters that control aspects for \code{\link{tvcm}}.  
}

\usage{
tvcm_control(method = c("mob", "partreg"),
             alpha = 0.05, bonferroni = TRUE,
             maxwidth = ifelse(method == "partreg", 10L, Inf),
             minbucket = 50L, mtry = Inf,
             minsplit = 2 * minbucket, trim = 0.1,
             maxdepth = Inf, maxstep = maxwidth - 1L, 
             nselect = Inf, estfun = list(),  
             maxevalsplit = 20, riskfun = neglogLik,
             fast = 0L, verbose = FALSE, ...)
}

\arguments{
  \item{method}{single character. Defines the algorithm to be used.}
  \item{alpha}{numeric significance threshold between 0 and 1. A node is
    splitted when the smallest (possibly Bonferroni-corrected) \eqn{p}
    value for any coefficient coefficient test in the current step falls
    below \code{alpha}.}  
  \item{bonferroni}{logical. Should \eqn{p} values be Bonferroni
    corrected?}
  \item{maxwidth}{integer. The maximum width of the tree.}  
  \item{minbucket}{numeric. The minimum sum of weights in a terminal
    node. Specifies also the trimming in parameter coefficient tests for
    numeric variables.}
  \item{mtry}{integer. The number of input variables randomly sampled as
    candidates at each node for random forest like algorithms.}
  \item{minsplit}{numeric. The minimum sum of weights in a node in order
    to be considered for splitting.}
  \item{trim}{numeric between 0 and 1. Specifies the trimming parameter in 
    coefficient constancy tests for continuous partitioning variables.}	
  \item{maxdepth}{integer. The maximum depth of the tree.}
  \item{maxstep}{integer. The maximum number of iterations. Equals to
    \code{"maxwidth"} minus 1.} 
  \item{nselect}{integer. The maximum number of distinct variables to be
    selected.}
  \item{estfun}{list of arguments to be passed to \code{\link{gefp.olmm}}.}
  \item{maxevalsplit}{integer. The maximum number of evaluation for
    splits on numeric partitioning variables.}
  \item{riskfun}{function. A function that extracts a risk measure from a
    fitted object, e.g., the negative log likelihood (default).}
  \item{fast}{integer. Whether the algorithm should be
    accelerated. Currently, for values larger than 1, optimizations use
    coefficients from the most actual estimation as initial values.}
  \item{verbose}{logical. Should information about the fitting process
    (such as test statistics, \eqn{p} values, 
    selected splitting variables and split points) be printed to the
    screen?}
  \item{...}{further arguments to be passed. For \code{tvcm_control}, these 
    are undocumented, see details. For \code{fvcm_control}, the arguments 
    are passed to \code{tvcm_control}.} 
}

\value{
  A list of class \code{tvcm_control} containing
  the control parameters.
}

\author{Reto Buergin}

\seealso{\code{\link{tvcm}}, \code{\link{fvcm}}}

\examples{
tvcm_control(method = "partreg")
}