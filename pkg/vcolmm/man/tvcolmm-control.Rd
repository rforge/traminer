\name{tvcolmm-control}

\alias{tvcolmm_control}
\alias{fvcolmm_control}

\title{Control parameters for \code{\link{tvcolmm}} and \code{\link{fvcolmm}}.} 

\description{
  Various parameters that control aspects of recursive partitioning for
  the function \code{\link{tvcolmm}}. 
}

\usage{
tvcolmm_control(alpha = 0.05, bonferroni = TRUE,
                minsplit = 50L, minbucket = 25L,
                maxdepth = Inf, maxwidth = Inf,
                mtry = Inf, nselect = Inf,
                estfun = list(), sctest = TRUE, 
                maxevalsplit = 20, lossfun = neglogLik,
                fast = 0L, verbose = FALSE,...)

fvcolmm_control(alpha = 1.0, bonferroni = TRUE,
                minsplit = 50L, minbucket = 25L,
                maxdepth = Inf, maxwidth = 10L,
                mtry = 5L, nselect = Inf,
                estfun = list(), sctest = TRUE, 
                maxevalsplit = 20, lossfun = neglogLik,
                fast = 0L, verbose = FALSE,...)
}

\arguments{
  \item{alpha}{numeric significance level. A node is splitted when the
    smallest (possibly Bonferroni-corrected) \eqn{p} value for any
    coefficient coefficient test in a iteration falls below \code{alpha}.} 
  \item{bonferroni}{logical. Should \eqn{p} values be Bonferroni
    corrected?}
  \item{minsplit}{numeric. The minimum sum of weights in a node in order
    to be considered for splitting.} 
  \item{minbucket}{numeric. The minimum sum of weights in a terminal
    node. Specifies also the trimming in parameter coefficient tests for
    numeric variables.} 
  \item{maxdepth}{integer. The maximum depth of the tree.}
  \item{maxwidth}{integer. The maximum width of the tree.}
  \item{mtry}{integer. The number of input variables randomly sampled as
    candidates at each node for random forest like algorithms.}
  \item{nselect}{integer. The maximum number of distinct variables to be
    selected.}
  \item{estfun}{list of arguments to be passed to \code{\link{gefp.olmm}}.}
  \item{sctest}{logical. Whether coefficient constancy tests are used
    for variable- and node selection.}
  \item{maxevalsplit}{integer. The maximum number of evaluation for
    splits on numeric partitioning variables.}
  \item{lossfun}{function. A function that extracts a loss measure from a
    fitted \code{\link{olmm}} object.}
  \item{fast}{integer. Whether the algorithm should be
    accelerated. Currently, for values larger than 1, optimizations use
    coefficients from the most actual estimation as initial values.}
  \item{verbose}{logical. Should information about the fitting process
    (such as test statistics, \eqn{p} values, 
    selected splitting variables and split points) be printed to the
    screen?}
  \item{...}{further undocumented arguments to be passed. See details.}
}

\details{\code{lossfun} is the function to be minimized by the
  algorithm. By default, this is the likelihood function times minus one
  (\code{\link{neglogLik}}). User-specified must be applicable on
  \code{\link{olmm}} objects and yield a numeric scalar.
}

\value{
  A list of class \code{tvcolmm_control} containing
  the control parameters.
}

\author{Reto Buergin}

\seealso{\code{\link{tvcolmm}}, \code{\link{fvcolmm}}}
