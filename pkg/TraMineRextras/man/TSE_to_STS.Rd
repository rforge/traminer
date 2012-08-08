\name{TSE_to_STS}
\alias{TSE_to_STS}

\title{
Conversion from TSE data format (for event sequences) to STS (state sequence) data format.
}
\description{
Conversion from TSE data format (for event sequences) to STS (state sequence) data format.
}
\usage{
TSE_to_STS(seqdata, id = 1, timestamp = 2, event = 3, stm = NULL, tmin = 1, 
    tmax = NULL, firstState = "None")
}

\arguments{
  \item{seqdata}{a data frame or matrix containing sequence data.}
  \item{id}{Name or index of the column containing the identification numbers for the sequences.}
  \item{timestamp}{Name or index of the column containing the timestamp of the events.}
  \item{event}{Name or index of the column containing the events.}
  \item{stm}{An event to state transition matrix (See \code{\link{seqe2stm}}).}
  \item{tmin}{Integer. Starting time of the state sequence.}
  \item{tmax}{Integer. Ending time of the state sequence.}
  \item{firstState}{Character. The name of the state before any events has occured.}
}
\details{
Convert data in TSE format (for event sequences) to STS (state sequence). By default, the states are defined has the combination of events that already has occured.
Different scheme may be specified using the function \code{\link{seqe2stm}} and the \code{stm} argument.
}
\value{
A \code{data.frame} with the sequence in STS format.
}
\author{
Matthias Studer
}

\note{This function is a pre-release and further testing is still needed, please report any problems.}

\seealso{
See Also \code{\link{seqe2stm}}, \code{\link[TraMineR]{seqformat}}.
}
\references{
Ritschard, G., Gabadinho, A., Studer, M. & Müller, N.S. (2009), "Converting between various sequence representations",
In Ras, Z. & Dardzinska, A. (eds) \emph{Advances in Data Management}. Series: \emph{Studies in Computational Intelligence}. Volume 223, pp. 155-175. Berlin: Springer.
}
\examples{
data(actcal.tse)
events <- c("PartTime", "NoActivity", "FullTime", "LowPartTime")
## Dropping all previous events.
stm <- seqe2stm(events, dropList=list("PartTime"=events[-1], 
   NoActivity=events[-2], FullTime=events[-3], LowPartTime=events[-4]))
mysts <- TSE_to_STS(actcal.tse, id=1, timestamp=2, event=3, 
   stm=stm, tmin=1, tmax=12, firstState="None")
}

\keyword{data format}
