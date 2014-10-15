\name{tvcm-control}

\alias{tvcm_control}

\title{Control parameters for \command{\link{tvcm}}.} 

\description{
  Various parameters that control aspects for \command{\link{tvcm}}.  
}

\usage{
tvcm_control(minsize = 30, mindev = 2.0,
             sctest = FALSE, alpha = 0.05, bonferroni = TRUE,
             trim = 0.1, estfun.args = list(), nimpute = 5, 
             maxnomsplit = 5, maxordsplit = 9, maxnumsplit = 9,
             maxstep = 1e3, maxwidth = 1e9, maxdepth = 1e9,
             lossfun = neglogLik2, ooblossfun = NULL,
             cp = 0.0, dfpar = 0.0, dfsplit = 1.0,
             cv = !sctest, folds = folds_control("kfold", 5),
             prune = cv, papply = mclapply, papply.args = list(),
             center = TRUE, verbose = FALSE, ...)
}

\arguments{
  \item{minsize}{numeric (vector). The minimum sum of weights in
    terminal nodes.}
  \item{mindev}{numeric scalar. The minimum permitted training error 
    reduction a split must exhibit to be considered of a new split.
    The main role of this parameter is to save computing time by early
    stopping. May be set lower for very few partitioning variables
    resp. higher for many partitioning variables. }
  \item{sctest}{logical scalar. Defines whether coefficient constancy
    tests should be used for variable and node selection.}
  \item{alpha}{numeric significance threshold between 0 and 1. A node is
    splitted when the smallest (possibly Bonferroni-corrected) \eqn{p}
    value for any coefficient constancy test in the current step falls
    below \code{alpha}.}  
  \item{bonferroni}{logical. Indicates if and how \eqn{p}-values of
    coefficient constancy tests must be Bonferroni
    corrected. See details.}
  \item{trim}{numeric between 0 and 1. Specifies the trimming parameter
    in coefficient constancy tests for continuous partitioning
    variables.} 	
  \item{estfun.args}{list of arguments to be passed to
    \command{\link{gefp.olmm}}. The argument is specific for models
    fitted with \command{\link{olmm}} with random effects.}
  \item{nimpute}{the number of times coefficient constancy tests
    should be repeated in each iteration. See details.}
  \item{maxnomsplit}{integer. For nominal partitioning variables with
    more the \code{maxnomsplit} the categories are ordered an treated as
    ordinal.}
  \item{maxordsplit}{integer. The maximum number of splits of ordered
    partitioning variables to be evaluated.}
  \item{maxnumsplit}{integer. The maximum number of splits of numeric
    partitioning variables to be evaluated.}
  \item{maxstep}{integer. The maximum number of iterations i.e. number
    of splits to be processed.}  
  \item{maxwidth}{integer (vector). The maximum width of the partition(s).}  
  \item{maxdepth}{integer (vector). The maximum depth of the partition(s).} 
  \item{lossfun}{a function to extract the training error, typically
    minus two times the negative log likelihood of the fitted model (see
    \command{\link{neglogLik2}}).}  
  \item{ooblossfun}{a loss function that defines how to compute the
    validation error during cross-validation. The function will be
    assigned to the \code{fun} argument of \command{\link{oobloss}}.}
  \item{cp}{numeric scalar. The penalty to be multiplied with the
    complexity of the model during partitioning. The complexity of the
    model is defined as the number of coefficients times \code{dfpar}
    plus the number of splits times \code{dfsplit}. By default, \code{cp
      = 0} (no penalization during partitioning) and \code{dfpar = 0} and
    \code{dfsplit = 1} (the complexity is measured as the total number
    of splits). \code{cp} also presents the minimum evaluated value at
    cross-validation.}
  \item{dfpar}{numeric scalar. The degree of freedom per model
    coefficient. Is used to compute the complexity of the model, see
    \code{cp}.}  
  \item{dfsplit}{a numeric scalar. The degree of freedom per split. Is
    used to compute the complexity of the model, see \code{cp}.}   
  \item{cv}{logical scalar. Whether or not the \code{cp} parameter
    should be cross-validated. If \code{TRUE} \command{\link{cvloss}} is
    called.} 
  \item{folds}{a list of parameters to create folds as produced by 
    \command{\link{folds_control}}. Is used for cross-validation.}
  \item{prune}{logical scalar. Whether or not the initial tree should be
    pruned by the estimated \code{cp} parameter from
    cross-validation. Cannot be \code{TRUE} if \code{cv = FALSE}.} 
  \item{papply}{(parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. The function will parallelize the 
    partition stage and the evaluation of the cross-validation folds as
    well as the final pruning stage.}
  \item{papply.args}{a list of arguments to be passed to \code{papply}.}
  \item{center}{logical integer. Whether the predictor variables of
    update models during the grid search should be centered. Note that
    \code{TRUE} will not modify the predictors of the fitted model.} 
  \item{verbose}{logical. Should information about the fitting process
    be printed to the screen?}
  \item{\ldots}{further, undocumented arguments to be passed.} 
}

\details{
  The arguments \code{alpha}, \code{bonferroni}, \code{trim},
  \code{estfun.args}, \code{nimpute} are relevant only if \code{sctest =
  TRUE}. Moreover, the argument \code{nimpute} is relevant only for the
  special case where the model is fitted by \command{\link{olmm}}
  including random effects and the data are unbalanced (individuals have
  different numbers of observations). In these cases, the internally
  used \command{\link{estfun.olmm}} function to extract the
  (decorrelated) likelihood scores will balance out the data by
  generating observations from the model under the \eqn{H_0} of the
  test. To account for the randomness of generating obervations, the
  test procedure is \code{nimpute} times repeated. The final p-values
  are then the averages of the repetations. The procedure is explained
  in detail in Buergin and Ritschard (2014a). 
}

\value{
  A list of class \code{tvcm_control} containing
  the control parameters for \command{\link{tvcm}}.
}

\author{Reto Buergin}

\seealso{\command{\link{tvcm}}, \command{\link{fvcm}}}

\examples{
tvcm_control()
}