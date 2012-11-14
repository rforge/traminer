\name{seqentrans}
\alias{seqentrans}
%
\author{Nicolas Müller and Gilbert
Ritschard}
%
\title{Event sequence length and number of events}
%
\description{
  Adds the sequence length (number of transitions) and total number of events of event sequences to the data attribute of a \code{subseqelist} event sequence object.
}
\usage{
  seqentrans(fsubseq)
}
\arguments{
  \item{fsubseq}{ A \code{subseqelist} object as returned by \code{\link[TraMineR]{seqefsub}}.}
}
\details{
  An event sequence object is an ordered list of transitions, with each transition a non-ordered list of events occurring at a same position.
}
\value{
  The object \code{fsubseq} updated with the additional information.
}
%\seealso{
%}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse[1:500,])

##Searching for frequent subsequences appearing at least 30 times
fsubseq <- seqefsub(actcal.seqe, minSupport=10)
fsubseq <- seqentrans(fsubseq)
## dispaying only those with at least 3 transitions
fsubseq[fsubseq$data$ntrans>2]
## dispaying only those with at least 3 events
fsubseq[fsubseq$data$nevent>2]
}
\keyword{util}
