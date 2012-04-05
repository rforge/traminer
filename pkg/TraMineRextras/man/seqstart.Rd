\name{seqstart}
\alias{seqstart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aligning sequence data to a new start time.
}
\description{
Aligning sequence data to a new start time.
}
\usage{
seqstart(seqdata, data.start, new.start, tmin = NULL, tmax = NULL, missing = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a data frame or matrix containing sequence data.}
  \item{data.start}{Integer. The actual starting date of each sequence.}
  \item{new.start}{Integer. The new starting date to use.}
  \item{tmin}{Integer. Minimum time in the new sequence data. If \code{NULL}, it is guessed from the data.}
  \item{tmax}{Integer. Minimum time in the new sequence data. If \code{NULL}, it is guessed from the data.}
  \item{missing}{Character. Code used to fill missing data in the new time axis.}
}
\value{A matrix.}

\author{Matthias Studer}

\note{Warning: This function needs further testing.}

\examples{
#An example data set
surveydata <- matrix(c("A" ,"A" , "B" , "B" , "B",
"A" , "A" , "B" , "B" , "B", 
"A" , "A",  "B" , "B" , "B" ,
"A"  ,"A" , "A" , "B"  ,"B" ,
"A"  ,"A" , "A" , "A" , "B"), byrow=TRUE, ncol=5)

print(surveydata)

## Usefull to change from survey time to age
surveytime <- 2000
birthdate <- 1995:1999
agedata <- seqstart(surveydata, data.start=surveytime,  new.start=birthdate)
print(agedata)
## Setting minimum and maximum age explicitly.
seqstart(surveydata, data.start=surveytime,  new.start=birthdate, tmax=8, tmin=4, missing="*")

## Usefull to change from age to survey time
surveyage <- surveytime - birthdate 
seqstart(agedata, data.start=1,  new.start=surveyage)

}
\keyword{data format}
