\name{seqeid}
\alias{seqeid}
\title{Retrieve id of an event sequence object.}
\description{
  Retrieve id of an event sequence or a list of event sequence object.
}
\usage{
seqeid(s)
}
\arguments{
\item{s}{A sequence or a list of sequence}
}
%\seealso{}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)
seqeid(actcal.seqe)
}
\keyword{attribute}

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Event Sequences}