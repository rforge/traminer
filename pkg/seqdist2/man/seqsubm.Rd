\encoding{latin1}

\name{seqsubm}

\alias{seqsubm}

\title{Create a substitution-cost matrix}

\description{
The function proposes different ways to generate substitution costs
(supposed to reflect state dissimilarities). Proposed methods are:
\code{"CONSTANT"} (same cost for all substitutions), \code{"TRATE"} (derived from the observed transition rates), \code{"FUTURE"} (Chi-squared distance between conditional state distributions \code{lag} positions ahead), \code{"FEATURES"} (Gower distance between state features), \code{"INDELS"}, \code{"INDELSLOG"} (based on estimated indel costs)
The substitution-cost matrix is intended to serve as \code{sm} argument in the \code{\link{seqdist}} function that computes distances between sequences.
 }

\details{
Actually this function calls \code{\link{seqcost}} and extracts the substitution cost matrix from the object returned by \code{seqcost}. See \code{\link{seqcost}} help page for details about the available methods.  
}


\usage{
seqsubm(...)
}

 \arguments{
  \item{...}{arguments passed to \code{\link{seqcost}}.}

}

\value{
 A matrix of size \eqn{ns*ns}, where
\eqn{ns} is the number of states in the alphabet of the sequence object.
}

\seealso{
 \code{\link{seqcost}}.
}

\examples{
## Defining a sequence object with columns 10 to 25
## in the 'biofam' example data set
data(biofam)
biofam.seq <- seqdef(biofam,10:25)

## Optimal matching using transition rates based substitution-cost matrix
## and insertion/deletion costs of 3
trcost <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq,method="OM",indel=3,sm=trcost)

## Optimal matching using constant value (2) substitution-cost matrix
## and insertion/deletion costs of 3
ccost <- seqsubm(biofam.seq, method="CONSTANT", cval=2)
biofam.om.c2 <- seqdist(biofam.seq, method="OM",indel=3,sm=ccost)

## Displaying the distance matrix for the first 10 sequences
biofam.om.c2[1:10,1:10]

## =================================
## Example with weights and missings
## =================================
data(ex1)
ex1.seq <- seqdef(ex1,1:13, weights=ex1$weights)

## Unweighted
subm <- seqsubm(ex1.seq, method="TRATE", with.missing=TRUE, weighted=FALSE)
ex1.om <- seqdist(ex1.seq, method="OM", sm=subm, with.missing=TRUE)

## Weighted
subm.w <- seqcost(ex1.seq, method="TRATE", with.missing=TRUE, weighted=TRUE)
ex1.omw <- seqdist(ex1.seq, method="OM", sm=subm.w, with.missing=TRUE)

ex1.om == ex1.omw
}

\author{Matthias Studer and Alexis Gabadinho (first version) (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity measures}
