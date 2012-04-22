\encoding{latin1}
\name{dissvar.grp}
\alias{dissvar.grp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Discrepancy by group.
}
\description{
This function computes the dissimilarity-based discrepancy measure of the groups defined by the group variable. The function is a wrapper for the TraMineR \code{\link{dissvar}} function.
}
\usage{
dissvar.grp(mdis, group=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdis}{
        a dissimilarity matrix or a \code{dist} object.
        }
  \item{group}{
        group variable. If \code{NULL} a single group is assumed.
        }
  \item{...}{
        additional arguments passed to \code{\link{dissvar}} .
        }
}

\details{
The function is a wrapper for running \code{\link{dissvar}} on the different groups defined by the \code{group} variable.
}
\value{
A vector with the group discrepancy measures.
}

\note{This function is a pre-release and further testing is still needed, please report any problems.
}
\author{
Gilbert Ritschard
}

\seealso{
\code{\link{dissvar}}
}
\examples{
data(mvad)
mvad.shortlab <- c("EM","FE","HE","JL","SC","TR")
mvad.seq <- seqdef(mvad[, 17:86], states=mvad.shortlab, labels=mvad.shortlab)
dist <- seqdist(mvad.seq, method="LCS")

## discrepancy based on non-squared dissimilarities
dissvar.grp(dist, mvad$gcse5eq) 
## square root of discrepancy based on squared dissimilarities
sqrt(dissvar.grp(dist, mvad$gcse5eq, squared=TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{state sequences}

