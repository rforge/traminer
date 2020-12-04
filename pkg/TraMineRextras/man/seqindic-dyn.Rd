\name{seqindic.dyn}
\alias{seqindic.dyn}
\alias{print.dynin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic index
%%  ~~function to do ... ~~
}
\description{Dynamic (i.e. successive) values of an individual index. For each sequence, the values of the selected index is computed on sliding windows.
}
\usage{
seqindic.dyn(seqdata, indic="cplx", window.size = .2,
      with.missing=TRUE, ...)
}
\arguments{
  \item{seqdata}{state sequence object (\code{stslist}) as produced by \code{seqdef}
  }
  \item{indic}{character string: the individual index. Can be any value supported by \code{\link{seqindic}} except index group names.
  }
  \item{window.size}{Integer or real. If an integer > 1, window size in number of positions. If real number in the range ]0,1), the window size is set as that proportion of the length of the longest sequence.
  }
  \item{with.missing}{logical. Should the missing state be treated as a state of the alphabet?
  }
  \item{...}{additional arguments passed to \code{\link{seqindic}}
  }
}

\details{
The values of the \code{indic} index are computed for each sequence on sliding windows of size \code{window.size}.
}

\value{
A matrix of class \code{dynin} with attributes \code{xtstep}, \code{tick.last}, \code{window.size}, and \code{indic}. The first two as well as the row and column names are taken from \code{seqdata}. Column names indicate the middle position of the window.

There are \code{print} and \code{plot} methods for \code{dynin} objects. See \code{\link{plot.dynin}}.
}

\seealso{\code{\link{seqindic}}, \code{\link{plot.dynin}}
}

\references{
Pelletier, D.; Bignami-Van Assche, S. & Simard-Gendron, A. (2020) Measuring Life Course Complexity with Dynamic Sequence Analysis, \emph{Social Indicators Research} \doi{10.1007/s11205-020-02464-y}.

%% ~put references to the literature/web site here ~
}

\author{
Gilbert Ritschard
}

\examples{
data(actcal)
cases <- 1:100
actcal <- actcal[cases,] ## Here, only a subset
actcal.seq <- seqdef(actcal[,13:24], alphabet=c('A','B','C','D'))

a.dyn <- seqindic.dyn(actcal.seq, indic='cplx', with.missing=FALSE)

plot(a.dyn, group=actcal[cases,'sex'], fstat=mean)

## Trimmed mean
plot(a.dyn, group=actcal[cases,'sex'], fstat=function(x)mean(x, trim=.02))



}