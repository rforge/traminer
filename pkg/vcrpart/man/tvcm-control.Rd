\name{tvcm-control}

\alias{tvcm_control}

\title{Control parameters for \code{\link{tvcm}}.} 

\description{
  Various parameters that control aspects for \code{\link{tvcm}}.  
}

\usage{
tvcm_control(lossfun = neglogLik2, 
             maxstep = Inf, maxwidth = Inf,
             minsize = 30, maxdepth = Inf,
             dfpar = 2.0, dfsplit = 0.0,
             maxoverstep = ifelse(sctest, Inf, 0),
             sctest = FALSE, alpha = 0.05, bonferroni = TRUE,
             trim = 0.1, estfun = list(), ninpute = 5L,
             maxfacsplit = 5L, maxordsplit = 10, maxnumsplit = 10,
             cv = !sctest, folds = folds_control("kfold", 5),
             prune = cv, keeploss = FALSE, papply = mclapply,
             verbose = FALSE, ...)

}

\arguments{
  \item{lossfun}{a function that extracts a loss measure from a
    fitted object, e.g., two times the negative log likelihood
    (default).}
  \item{maxstep}{integer. The maximum number of iterations i.e. the
    total number of splits processed.} 
  \item{maxwidth}{integer (vector). The maximum width of the tree(s).}  
  \item{minsize}{numeric. The minimum sum of weights in a terminal
    node. The default is the number of varying coefficients times 10.
    The parameter specifies also the trimming in parameter coefficient
    tests for numeric variables (if \code{sctest = TRUE}).}
  \item{maxdepth}{integer (vector). The maximum depth of the tree(s).}
  \item{dfpar}{a numeric scalar larger than zero. The per-parameter penalty to
    be applied for stopping. See also argument \code{maxoverstep}.}
  \item{dfsplit}{a numeric scalar larger than zero. The per-split penalty to
    be applied for stopping. See also argument \code{maxoverstep}.}
  \item{maxoverstep}{integer scalar. The maximum number of consecutive
    times the penalized reduction statistic is allowed to be smaller than 
    \code{dfsplit} before stopping. Specifically, the penalized loss 
    is computed as the loss (see argument \code{lossfun}) of the model plus 
    \code{dfpar} times the number of coefficients  (as extracted from
    \code{\link{extractAIC}}), and the penalized loss reduction statistic is the
    difference between the penalized loss of the current model and the 
    penalized loss of the selected (best) update model. The argument is 
    disabled if \code{sctest = TRUE}.}  
  \item{sctest}{logical scalar. Defines whether coefficient constancy
    tests should be used for variable and node selection.}
  \item{alpha}{numeric significance threshold between 0 and 1. A node is
    splitted when the smallest (possibly Bonferroni-corrected) \eqn{p}
    value for any coefficient constancy test in the current step falls
    below \code{alpha}.}  
  \item{bonferroni}{logical. Indicates if and how \eqn{p}-values of
    coefficient constancy tests must be Bonferroni
    corrected. See details.}
  \item{trim}{numeric between 0 and 1. Specifies the trimming parameter in 
    coefficient constancy tests for continuous partitioning variables.}	
  \item{estfun}{list of arguments to be passed to \code{\link{gefp.olmm}}.}
  \item{ninpute}{the number of times the coefficient constancy tests are
    repeated in each node in case of two-stage data that are
    unbalanced.}
  \item{maxfacsplit}{integer.}
  \item{maxordsplit}{integer.}
  \item{maxnumsplit}{integer. The maximum number of evaluation for
    splits on numeric partitioning variables.}
  \item{cv}{logical scalar. Whether or not the \code{dfsplit} parameter should
    be cross-validated.}
  \item{folds}{a list of parameters to create folds as produced by 
   \code{\link{folds_control}}.}
  \item{prune}{logical scalar. Whether or not the overly large partitions should
    be pruned by the estimated \code{dfsplit} parameter from cross-validation. 
    Note that \code{prune = TRUE} conflicts with \code{cv = FALSE}}
  \item{keeploss}{logical scalar or a numeric equal or larger than
    0. Indicates if and how many times the computed penalized loss reduction
    statistics should be reused in the following
    iterations. Specifically, the option activates approximating the penalized loss
    reduction statistic of a split based on the penalized loss reduction statistic
    of the same split in the previous iteration. Values different from
    \code{FALSE} or 0 will accelerate the algorithm at the expense of
    accuracy. In cases of multiple \code{\link{vc}} terms, only the
    statistics corresponding to the \code{\link{vc}} term which was
    splitted in the previous iteration are reused.} 
   \item{papply}{(parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. The function will parallelize the 
    partition stage and the evaluation of the cross-validation folds as
    well as the final pruning stage.}  
  \item{verbose}{logical. Should information about the fitting process
    (such as test statistics, \eqn{p} values, 
    selected splitting variables and split points) be printed to the
    screen?}
  \item{\ldots}{further, undocumented arguments to be passed. These can
    include arguments for the \code{papply} function.} 
}

\value{
  A list of class \code{tvcm_control} containing
  the control parameters.
}

\author{Reto Buergin}

\seealso{\code{\link{tvcm}}, \code{\link{fvcm}}}

\examples{
tvcm_control()
}