\name{STS_to_SPELL}
\alias{STS_to_SPELL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data conversion from STS to SPELL.
}
\description{
Convert data from STS to SPELL.
}
\usage{
STS_to_SPELL(seqdata, id = NULL, pdata = NULL, birthdate = NULL, with.missing = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{
a data frame or matrix containing sequence data in STS format
}
  \item{id}{
Either a column name (or number) of \code{pdata} or a vector of individual identification number or character. If \code{id=NULL}, \code{id} is set to \code{1:nrow(seqdata)} or \code{rownames} of \code{seqdata} (if not \code{NULL}). 
}
  \item{pdata}{
A data frame containing personal \code{id} and \code{birthdate}.
}
  \item{birthdate}{
Either a column name (or number) of \code{pdata} or a vector of individual process start date.
}
  \item{with.missing}{If \code{TRUE}, missing episodes are included. Otherwise, these episodes are ignored.}

}
\details{
Spell data format has the following caracteristics:
	- One row per episode (spell)
	- Each spells are specified using four variables: an id, a begin date, a end date, and the states.
}
\value{A \code{data.frame} with the sequence data in the SPELL format.
}
\author{
Matthias Studer
}

\note{This function is a pre-release and further testing is still needed, please report any problems.}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link[TraMineR]{seqformat}}.
}
\examples{
data(biofam)

## Create the sequence object
bfstates <- c("Parent", "Left", "Married", "Left/Married",  "Child", "Left/Child", "Left/Married/Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
bf.seq <- seqdef(biofam[ ,10:25], alphabet=0:7, states=bf.shortlab, labels=bfstates)

spell <- STS_to_SPELL(bf.seq, birthdate=biofam$birthyr)
head(spell)

}
