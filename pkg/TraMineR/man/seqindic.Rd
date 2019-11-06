\name{seqindic}
\alias{seqindic}
%
\author{Gilbert Ritschard}
%
\title{Sequence indicators}
%
\description{
  Returns a table with per sequence values of selected indicators.
}
\usage{
seqindic(seqdata, indic=c("visited","trans","ient","cplx"), with.missing=FALSE,
    ipos.args=list(), prec.args=list())
}
\arguments{
  \item{seqdata}{a state sequence object (\code{stslist}) as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{indic}{vector of character strings. List of selected indicators among \code{"lgth"} (sequence length), \code{"nonm"} (number of non-missing elements), \code{"dlgth"}, (number of successive distinct states), \code{"visited"} (number of visited states), \code{"trans"} (number of transitions), "ptrans" (number of transitions as a proportion of maximum possible transitions), \code{"ient"} (longitudinal normalized entropy), \code{"cplx"} (complexity index), \code{"turb"} (turbulence), \code{"nturb"} (normalized turbulence), \code{"ipos"} (proportion of positive states), \code{"prec"} (precarity index). Setting \code{indic="all"} selects all indicators but \code{"ipos"} and \code{"prec"} that require an additional argument list.}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}

  \item{ipos.args}{list: when \code{"ipos"} requested, the arguments passed to \code{seqipos} together with \code{seqdata}. See \code{\link{seqipos}}}
      
  \item{prec.args}{list: when \code{"prec"} requested, the arguments passed to \code{seqprecarity} together with \code{seqdata}. See \code{\link{seqprecarity}}}
      
}
\details{
  The number of visited states is the number of different elements in the sequence, e.g. 2 for aababba.
  
  The sequence length, number of transitions, longitudinal entropy, complexity, and turbulence are computed respectively with functions \code{\link{seqlength}}, \code{\link{seqtransn}}, \code{\link{seqient}}, \code{\link{seqici}}, and \code{\link{seqST}}. See corresponding help pages for details.
  
  Number of distinct successive states is one plus the number of transitions.
  
  The \code{with.missing} argument applies to all indicators but the length. \code{lgth} returns the length obtained with \code{with.missing=TRUE}, and \code{nonm} the length obtained with \code{with.missing=FALSE}.
  
  The precarity (indicator \code{"prec"}) is computed with function \code{seqprecarity} of TraMineRextras.
}
\value{
  A data frame with the selected indicators. Names are:
  \cr Lght: Length of the sequence
  \cr NonM: Number of non missing elements
  \cr Dlgth: Length of the sequence of distinct successive states (DSS)
  \cr Nvis: Number of visited states
  \cr Ntrn: Number of transitions
  \cr Ptrn: Number of transitions as a proportion of maximum number of transitions
  \cr Entr: Longitudinal entropy
  \cr Cplx: Complexity
  \cr Turbn: Normalized turbulence
  \cr Turb: Turbulence

}
\seealso{
\code{\link{seqlength}}, \code{\link{seqtransn}}, \code{\link{seqient}}, \code{\link{seqici}}, and \code{\link{seqST}}
}
\examples{
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
seqindic(sx, indic=c("lgth","nonm","visited","nturb","cplx"), with.missing=TRUE)
seqindic(sx, indic=c("lgth","nonm","visited","nturb","cplx"), with.missing=FALSE)
seqindic(sx, indic=c("lgth","dlgth","ipos"),
          ipos.args=list(pos.states=c("A","B")),with.missing=FALSE)
}
\keyword{util}
