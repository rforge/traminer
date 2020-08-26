\name{seqdur}
\alias{seqdur}
\title{Extract state durations from a sequence object.}
\description{
Extracts states durations from a sequence object. Returns a matrix containing the states durations for the sequences. The states durations in \code{'D-D-D-D-A-A-A-A-A-A-A-D'} are 4,7,1. Distinct states can be extracted with the \code{\link{seqdss}} function.
 }
\usage{
seqdur(seqdata, with.missing=FALSE)
}
\arguments{
   \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
	\item{with.missing}{Should non-void missing values be considered as regular states? See Details.}

}
\details{
When \code{with.missing=FALSE} (default) missing values are ignored and a substring \code{AA***A} will be considered as a spell \code{AAA} of duration 3. When \code{with.missing=TRUE}, durations are also computed for spells of missing values (gaps in sequences).

See \link{seqdef} on options for handling missing values when creating sequence objects.
}

\value{a matrix containing the states durations for each distinct state in each sequence.}

\seealso{
\code{\link{seqdss}}.
}

\examples{
## Creating a sequence object with the columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Retrieving the DSS
actcal.dur <- seqdur(actcal.seq)

## Displaying the durations for the first 10 sequences
actcal.dur[1:10,]
}

\author{Alexis Gabadinho and Gilbert Ritschard}

\keyword{Longitudinal characteristics}
