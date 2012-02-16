\name{ex2}
\docType{data}
\alias{ex2}
\title{Example data sets with weighted and unweighted sequence data}
\description{
Example data set used to demonstrate the handling of weights. The 'ex2.weighted' data set contains 6 sequences whith weights inflating to 100 sequences (sum of weights is 100), the same as in the 'ex2.unweighted' data set. 

The sequences are in the 'seq' column, weights (in the 'ex2.weighted' data set) are in the 'weight' column.

The alphabet is made of four possible states: A, B, C and D.
\cr

These data sets are intended mainly to test and illustrate the handling of weights in TraMineR's functions. Weighted results obtained with 'ex2.weighted' data set should be exactly the same as unweighted results obtained with the 'ex2.unweighted' data set.
}
\usage{data(ex2)}
\format{ex2.weighted: a data frame with 6 rows, 1 variable containing sequences as character strings, 1 weight variable.
\cr
ex2.unweighted: a data frame with 100 rows, 1 variable containing sequences as character strings.

}

\examples{
data(ex2)

ex2w.seq <- seqdef(ex2.weighted, 1, weights=ex2.weighted$weight)
ex2u.seq <- seqdef(ex2.unweighted)
}


\source{The brain of the TraMineR package maintainer.}

\keyword{datasets}
