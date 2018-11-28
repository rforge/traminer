\encoding{utf8}

\name{seqprecstart}
\alias{seqprecstart}

\title{State precarity}

\description{
Determines the initial state precarity level from the state ordering. (See details.)
}

\usage{
seqprecstart(seqdata, state.order=alphabet(seqdata), state.equiv=NULL)
}

\arguments{
  \item{seqdata}{a state sequence object (class \code{stslist}) as returned by the \code{\link[TraMineR]{seqdef}} function.}

  \item{state.order}{Vector of short state labels defining the order of the states. First the less precarious (most positive) state and then the other states in increasing precariousness order. States of the alphabet that are not included here define the non-comparable states.}

  \item{state.equiv}{List of state equivalence classes. The classes in the list are each set as the vector of the short labels of the states forming the class. }

}


\details{
The initial state precarity levels are set as equivalently spaced values between 0 and 1 assigned to the states in the specified order.

Equivalent states get the mean value of the states in the equivalence class. 

Incomparable states (those not listed on the state.order list nor included in an equivalent class with at least one state listed on the state order list) receive a value of .5, i.e., the average of all state costs.

Only one state per class needs to be on the state order list. If more than one member is on the state order list, then they should be listed consecutively. Note that currently no check is performed.

}



\value{

The vector of assigned precarity levels sorted according to the original order of states in the alphabet.


}

\seealso{
 \code{\link{seqprecarity}}.
}

\references{

Ritschard, G., Bussi, M., and O'Reilly, J. (2018), "An index of precarity for
measuring early employment insecurity", in G. Ritschard, and M. Studer, \emph{Sequence Analysis and Related Approaches: Innovative Methods and Applications}, Series Life Course Research and Social
Policies, Vol. 10, pp 279-295. Cham: Springer.


}

\examples{
## Defining a sequence object with columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal <- actcal[1:200,] ## Here, only a subset
actcal.seq <- seqdef(actcal[,13:24])

## state precarity levels using the original state order
seqprecstart(actcal.seq)

## Assuming A and B as equivalent regarding precarity
seqprecstart(actcal.seq, state.equiv=list(c('A','B')))

## and letting C be non-comparable
seqprecstart(actcal.seq, state.order=c("A","B","D"), state.equiv=list(c('A','B')))

}

\author{Gilbert Ritschard}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sequence analysis}
