\name{STS_to_SPELL}
\alias{STS_to_SPELL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data conversion from STS to SPELL format.
}
\description{
Convert data from STS to SPELL format.
}
\usage{
STS_to_SPELL(seqdata, id=NULL, pdata=NULL, birthdate=NULL, with.missing=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{
a state sequence object of type \code{"stslist"}.
}
  \item{id}{Either a
vector of sequence id's or name of the id column in \code{pdata}.
}
  \item{pdata}{
Data frame with id and birth dates.
}
  \item{birthdate}{
Either a
vector of sequence birth dates or name of the birth date column in \code{pdata}.
}
    \item{with.missing}{Logical. Should explicit spells of missing states be created?}  
}
%\details{
%}
\value{
A \code{data.frame} with the sequence in SPELL format.
}
\author{
Matthias Studer
}

\seealso{
See Also \code{\link[TraMineR]{seqformat}}.
}
\note{This function is a pre-release and further testing is still needed, please report any problems.}
%\examples{
%hspell <- data.frame(begin1=rep(1, 5), end1=c(2:5, NA), status1=1:5,
%                     begin2=c(3:6, NA), end2=rep(NA, 5), status2=5:1)
%sts <- HSPELL_to_STS(hspell, begin=c("begin1", "begin2"), end=c("end1", "end2"),
%                     status=c("status1", "status2"))
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data format}
