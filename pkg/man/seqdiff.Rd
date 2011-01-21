\encoding{latin1}
\name{seqdiff}
\alias{seqdiff}
\alias{print.seqdiff}
\alias{plot.seqdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Decompose the difference between groups of sequences}
\description{
  Decompose the difference between groups of sequences
}
\usage{
seqdiff(seqdata, group, cmprange = c(0, 1),
        seqdist_arg=list(method="LCS",norm=TRUE),
		with.missing = FALSE, weighted = TRUE, squared = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a state sequence object created with the \code{\link{seqdef}} function.}
  \item{group}{The group variable.}
  \item{cmprange}{The time range of the sliding window on which subsequences are compared.}
  \item{seqdist_arg}{List of arguments passed directly to seqdist.}
  \item{with.missing}{Logical. If \code{TRUE}, missing values are considered as an additional state. If \code{FALSE} subsequences with missing values are removed from the analysis.}
  \item{weighted}{Logical. If \code{TRUE}, \code{seqdiff} uses the weights specified in \code{seqdata}.}
  \item{squared}{Logical. If \code{TRUE} the dissimilarities are squared for computing the discrepancy.}
}
\details{
  The function analyses
  the sequence discrepancy explained by the \code{group} variable within a sliding
  time-window (of range defined by \code{cmprange}). At
  each timestamp, the method uses \code{\link{seqdist}} to compute a distance matrix over the time-window and then derives the explained discrepancy on that window with \code{\link{dissassoc}}.

  There are print and plot methods for the result returned.
}

\value{
  A \code{seqdiff} object, with the following items:

  \item{stat}{A \code{data.frame} with three statistics (PseudoF, PseudoR2 and PseudoT)
  for each timestamp of the sequence,
      see \code{\link{dissassoc}}}

  \item{discrepancy}{A \code{data.frame} with, at each time stamp, the discrepancy within each group defined by the
      \code{group} variable and for the whole population.}
}

\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2010)
  Discrepancy	analysis of complex objects using dissimilarities.
  In F. Guillet, G. Ritschard, D. A. Zighed and H. Briand (Eds.),
  \emph{Advances in Knowledge Discovery and Management},
  Studies in Computational Intelligence, Volume 292, pp. 3-19. Berlin: Springer.
	
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009)
  Analyse de dissimilarités par arbre d'induction. In EGC 2009,
  \emph{Revue des Nouvelles Technologies de l'Information}, Vol. E-15, pp. 7-18.

 }
\seealso{\code{\link{dissassoc}} to analyse the association with the whole sequence}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.diff <- seqdiff(mvad.seq, group=mvad$gcse5eq)
print(mvad.diff)
plot(mvad.diff)
plot(mvad.diff, stat="discrepancy")
}
