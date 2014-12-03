\name{contr.wsum}

\alias{contr.wsum}

\title{Contrast matrices}

\description{
  Returns a category-weighted contrast matrix
}

\usage{
contr.wsum(x, weights = rep.int(1.0, length(x)))
}

\arguments{
  \item{x}{a factor vector}
  \item{weights}{a vector of weights with the same length as \code{x}.}
}

\details{
  Computes a contrast matrix similar to \command{\link{contr.sum}}. The
  reference category is however weighted by the sum of weights of the
  other categories.
}

\value{A matrix with \code{nlevels(x)} rows and \code{nlevels(x)- 1}
  columns.}
\author{Reto Buergin}

\seealso{\command{\link{contr.sum}}}

\examples{
x <- factor(rep(LETTERS[1:3], c(10, 20, 30)))
contr.wsum(x)
}