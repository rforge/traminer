\docType{methods}
\name{prune}

\alias{prune}

\title{
  Pruning of regression tree objects 
}
\description{
  Determines a subtree of a tree object by recursively \code{snipping}
  off the least important splits, based on a criteria. 
}
\usage{
prune(tree, \dots)
}

\arguments{
  \item{tree}{fitted tree object}
  \item{\dots}{further arguments passed to or from other methods.}
}

\seealso{
  \code{\link{tvcolmm}}
}

\keyword{methods}