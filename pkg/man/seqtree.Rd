\name{seqtree}
\alias{seqtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Tree structured analysis of a state sequence object.
}
\description{
	Facility to produce a tree structured analysis of a sequence object. The function provides a simplified interface for applying \code{\link{disstree}} on state sequence objects.
}
\usage{
seqtree(formula, data = NULL, weighted = TRUE, minSize = 0.05, 
     maxdepth = 5, R = 1000, pval = 0.01, weight.permutation = "replicate", 
     seqdist_arg = list(method = "LCS", norm = TRUE), 
     diss = NULL, squared = FALSE, first = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula where the left hand side is a state sequence object (see \code{\link{seqdef}})  and the right hand specifies the candidate variables for partitioning the set of sequences.}
  \item{weighted}{Logical. If \code{TRUE}, use the weights of the state sequence object.}
  \item{data}{a data frame where variable in formula will be searched}
  \item{minSize}{minimum number of cases in a node, in percentage if less than 1.}
  \item{maxdepth}{maximum depth of the tree.}
  \item{R}{Number of permutations used to assess the significance of the split.}
  \item{pval}{Maximum p-value, in percent.}
  \item{weight.permutation}{Weights permutation method: "diss" (attach weigths to the dissimilarity matrix), "replicate" (replicate case according to the \code{weights} arguments), "rounded-replicate" (replicate case according to the rounded \code{weights} arguments), "random-sampling" (random assignment of covariate profiles to the objects using distributions defined by the weights.)}
  \item{seqdist_arg}{list of arguments directly passed to \code{\link{seqdist}}}
  \item{diss}{An optional dissimilarity matrix.}
  \item{squared}{Logical. If \code{TRUE}, the dissimilarity matrix is squared}
  \item{first}{Character. An optional variable name to force the first split.}
}
\value{
	a seqtree object.
}

\seealso{
\code{\link{disstree}}
}
