\name{seqe2stm}
\alias{seqe2stm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Definition of an events to states matrix.
}
\description{
This function create a matrix specifying for each state (given in row) the new state if a given event (given in column) happen.
}
\usage{
seqe2stm(events, dropMatrix = NULL, dropList = NULL, firstState = "None")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{events}{Character. The vector of all possible events.}
  \item{dropMatrix}{Logical matrix. Specifying the event to forget once a given event has occured.}
  \item{dropList}{List. Same as \code{dropMatrix} but using a list (often more convenient).}
  \item{firstState}{Character. Name of the first state, before any event has occured.}
}

\details{
This function create a matrix specifying for each state (given in row) the new state if a given event (given in column) happen. This matrix is needed by \code{\link{TSE_to_STS}}.
By default, a new state is created for each combination of events that already has occured.

\code{dropMatrix} and \code{dropList} allow to specify which events should be "forgotten" once an event has occured.
For instance, we may want to forget the "marriage" event once the event "divorce" has occured.

\code{dropMatrix} specifies for each event given in row, the previous events, given in column that should be forgotten.
\code{dropList} uses a list to specify the same things. The form is list(event1=c(..., events to forgets), event2=c(..., events to forgets)).
See example below.
}
\value{
A matrix.
}
\note{This function is a pre-release and further testing is still needed, please report any problems.}
\references{
Ritschard, G., Gabadinho, A., Studer, M. & Müller, N.S. (2009), "Converting between various sequence representations", 
In Ras, Z. & Dardzinska, A. (eds) \emph{Advances in Data Management}. Series: \emph{Studies in Computational Intelligence}. Volume 223, pp. 155-175. Berlin: Springer. 
}
\author{
Matthias Studer
}

\seealso{
\code{\link{TSE_to_STS}}
}
\examples{
## Achive the same using dropMatrix or dropList.
## List of possible events.
events <- c("marr", "child", "div")
dm <- matrix(FALSE, 3,3, dimnames=list(events, events))
dm[3, ] <- c(TRUE, TRUE, FALSE)
dm[1, 3] <- TRUE
## Using the matrix, we "forget" marriage and child events once "divorce" occured.
## We also forget "divorce" once "marriage" occured.
print(dm)
stm <- seqe2stm(events, dropMatrix=dm)

## We achive the same using dropList argument.
stmList <- seqe2stm(events, dropList=list("div"=c("marr", "child"), "marr"="div"))

## test that the results are the same
all.equal(stm, stmList)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data format}

