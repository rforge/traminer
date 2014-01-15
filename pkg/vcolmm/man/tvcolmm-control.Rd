\name{tvcolmm-control}

\alias{tvcolmm_control}

\title{Control parameters for \code{\link{tvcolmm}}.} 

\description{
  Various parameters that control aspects of recursive partitioning for
  the function \code{\link{tvcolmm}}. 
}

\usage{
tvcolmm_control(alpha = 0.05, bonferroni = TRUE,
                minsplit = 50L, trim = 0.1,
                lossfun = neglogLik, breakties = FALSE,
                terms = NULL, verbose = FALSE, ...)
}

\arguments{
  \item{alpha}{numeric significance level. A node is splitted when the
    (possibly Bonferroni-corrected) \eqn{p} value for any parameter
    stability test in that node falls below \code{alpha}.} 
  \item{bonferroni}{logical. Should \eqn{p} values be Bonferroni
    corrected?}
  \item{minsplit}{numeric. The minimum number of subjects represented in
    a node. If smaller than 1, it is interpreted as the fraction relative
    to the current node size.} 
  \item{trim}{numeric. This specifies the trimming in the parameter
    instability test for the numerical variables. If smaller than 1, it is
    interpreted as the fraction relative to the current node size.}
  \item{lossfun}{function. A function for extracting the minimized value
    of the loss function (a.k.a. objective- or risk functin) from a
    fitted model in a node.} 
  \item{breakties}{logical. Should ties in numeric variables be broken
    randomly for computing the associated parameter instability test?}
  \item{terms}{character vector. The names of the predictors to be
    tested on coefficient constancy.}
  \item{verbose}{logical. Should information about the fitting process
    (such as test statistics, \eqn{p} values, 
    selected splitting variables and split points) be printed to the
    screen?}
  \item{...}{further undocumented arguments to be passed. See details.}
}

\details{There are a few hidden argument allowing further control:
  \code{maxevalsplit} is the maximum of evaluated cutpoints in numeric
  effect modifiers, which is 20 by default. \code{depth} regularizes the
  maximum depth of the tree, \code{nselect} the maximum number of
  variables to select, \code{mtry} the number of randomly selected
  effect modifiers i.e. partitioning variables candidates in a node
  (i.e. the random feature selection mechanism in random forests). 

  \code{lossfun} is the function to be minimized by the algorithm. By
  default, this is the likelihood function times minus one
  (\code{\link{neglogLik}}). User-specified must be applicable on
  \code{\link{olmm}} objects and yield a numeric scalar.
}

\value{
  A list of class \code{tvcolmm_control} containing
  the control parameters.
}

\author{Reto Buergin}

\seealso{\code{\link{tvcolmm}}}
