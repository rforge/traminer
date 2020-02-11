\name{seqivolatility}
\alias{seqivolatility}
%
\author{Gilbert Ritschard}
%
\title{Volatility of individual state sequences}
%
\description{
  Returns Brzinsky-Fay's objective volatility of each sequence.
}
\usage{
  seqivolatility(seqdata, type=1, w=.5, with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a state sequence object (\code{stslist}) as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{type}{either 1 or 2. The default \code{type=1} is the indicator adjusted such that it can reach its bounds 0 and 1. Type 2 is the unadjusted indicator (see Details).}

  \item{w}{real in the range [0, 1]. Weight given to the proportion of states visited (see Details).}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}
}
\details{
  The (objective) volatility is the weighted average between the proportion \eqn{pvisited} of states visited and the frequency \eqn{ftrans} of transitions (state changes). Formally,

  \deqn{volatility = w \cdot pvisited + (1-w) \cdot ftrans}{volatility = w * pvisited + (1-w) * ftrans}

  The proportion of states visited is computed as \eqn{(visited - 1)/(|a| - 1}) when \code{type=1} and as \eqn{visited / |a|} when \code{type=2}. Here, \eqn{visited} is the number of states visited and \eqn{|a|} the size of the alphabet.

  The frequency of transition is \eqn{ftrans = \frac{transn}{max.transn}} where
  \eqn{transn} is the number of transitions (state changes) within the sequence, and \eqn{max.transn} the maximum possible transitions in the sequence.
  
  For the normative volatility computed for a set of positive states, see \code{\link{seqipos}}.
}

\value{
  A numeric vector with the volatility of each sequence.
}

\seealso{
\code{\link{seqintegration}}, \code{\link{seqipos}}, \code{\link{seqindic}}
}
\references{
Brzinsky-Fay, C. Unused Resources: Sequence and Trajectory Indicators. International Symposium on Sequence Analysis and Related Methods, Monte Verita, TI, Switzerland, Oct 10-11, 2018
}


\examples{
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")

seqivolatility(sx)
seqivolatility(sx, type=2)
seqivolatility(sx, with.missing=TRUE)
seqivolatility(sx, w=.7, with.missing=TRUE)

}

\keyword{Longitudinal characteristics}