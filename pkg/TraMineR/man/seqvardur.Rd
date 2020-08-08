\name{seqvardur}
\alias{seqvardur}
\alias{print.seqvardur}
%
\author{Gilbert Ritschard}
%
\title{Variance of spell durations in individual sequences}
%
\description{
  Variance of spell durations of individual state sequences.
}
\usage{
  seqvardur(seqdata, type=1, with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a state sequence object (\code{stslist}) as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{type}{either 1 or 2. The default \code{type=1} ignores non visited states. Type 2 accounts for the 0-time spent in non-visited states (see Details).}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}
}
\details{
  The returned variance is a population variance, i.e. obtained by dividing the sum of squares by the number of terms. 
  
  When \code{type=1}, this is the variance of the observed spells in the sequence. When \code{type=2}, in addition to the observed spells one spell of length 0 is considered for each non-visited state.
}

\value{
  A numeric vector with the spell duration variance of each sequence. 
  
  The returned object has a \code{vmax} attribute with the maximum value that the variance of each sequence could reach given its number of spells.
}

\seealso{
\code{\link{seqdur}}, \code{\link{seqindic}}
}
\references{
Ritschard, Gilbert (2020) Measuring the nature of state sequences, University of Geneva, preprint.
}


\examples{
sq.dat <- c('AAAA','AAAC','ABC','ABAA','AC')
sq <- seqdef(seqdecomp(sq.dat, sep=''), right=NA)

seqvardur(sq, type=1)
seqvardur(sq, type=1, with.missing=TRUE)
(vd2 <- seqvardur(sq, type=2))
attr(vd2,'vmax')
(vd2m <- seqvardur(sq, type=2, with.missing=TRUE))
attr(vd2m,'vmax')

}

\keyword{Longitudinal characteristics}
