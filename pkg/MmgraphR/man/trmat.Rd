\name{trmat}
\alias{trmat}
\alias{trmat.depmix.fitted}
\title{Extract a Probability Transition Matrix}
\description{
Extracts a probability transition matrix.
}
\usage{
## S3 method for class 'depmix.fitted'

trmat ( d )

}
\arguments{
	\item{d}{Object of class \command{\link[depmixS4]{depmix.fitted}}.}
}
\value{
	\code{trmat} returns an the probability transition matrix as an object of class \command{\link[base]{matrix}}.
}

\details{
The function \code{trmat} extracts the probability transition matrix from a \command{\link[depmixS4]{depmix.fitted}} object and returns it as a \command{\link[base]{matrix}}. \code{trmat} is used internally by the function \command{\link{trmatplot.depmix.fitted}} in order to extract the probability transition matrix to be plotted.
}

\author{
Pauline Poulcheria Adamopoulou
}
\seealso{
	\command{\link{trmatplot}},
	\command{\link{trmatplot.depmix.fitted}},  
	\command{\link[depmixS4]{depmix.fitted}}, 
}

