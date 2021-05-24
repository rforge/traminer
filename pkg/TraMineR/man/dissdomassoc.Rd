\name{dissdomassoc}
\alias{dissdomassoc}
\title{Domain association measures}

\description{
Measures of association between domains are computed as the association between the pairwise dissimilarities in the domains. Measures are: Pearson correlation, Spearman correlation, global Cronbach alpha, and Cronbach alpha for each subset of the domains. The function can also return the share of variance (R-square) of the dissimilarities in one domain that can be reproduced from the dissimilarities in the other domains.
}

\usage{
  dissdomassoc(domdiss, jointdiss = NULL, what = c("pearson","R2"),
        names=NULL, weights=NULL)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{domdiss}{List of symmetrical matrices or dist objects: the pairwise dissimilarities per domain.}
  \item{jointdiss}{\code{NULL} (default), matrix or dist object: pairwise dissimilarities for joint dimensions.}
  %\item{weights}{optional numerical vector containing weights.}
  \item{what}{String or vector of strings: requested association measures among \code{'pearson'}, \code{'spearman'}, \code{'R2'} \code{'cronbach'}, \code{'cron.subsets'}. Can also be \code{'all'} for all measures but \code{'spearman'}.}
  \item{names}{\code{NULL} (default) or string vector of length equal to number of domains. Names of domains. If \code{NULL}, set as \code{names(domdiss)}.}
  \item{weights}{vector of non-negative weights. If \code{NULL}, no weights are applied. (Currently Cronbach measures ignore weights!)}
}

\details{ The \code{dissdomassoc} function computes the domain association measures proposed by \cite{Piccarreta (2017)}.
These are for each pair of domains the correlation (Pearson or Spearman) between the domain specific pairwise dissimilarities. When \code{'R2'} is requested, the function computes the share of variance of the dissimilarities in one domain that is reproduced by the dissimilarities in all other domains. The Cronbach alpha measures the coherence between the domains. With \code{'cronbach.subsets'}, the Cronbach alpha is computed for the entire set of domains as well as for each possible subset of two or more domains.

When a \code{jointdiss} distance matrix or object is provided, correlations of each individual with this joint domain are also computed. The R2 of the joint domain takes account of all domains. However, R2's for the domains remain unchanged, i.e. they ignore the joint domain. Likewise, \code{jointdiss} does not affect the Cronbach alpha.

Spearman correlations are based on rank values of the distances. The computation of these ranks for weighted data is done with the \code{\link[cNORM]{weighted.rank}} function of the \code{cNORM} package and may take a while for large number (\code{> 1000}) of sequences.
}

\value{A list with the requested association measures
}
\references{
    Piccarreta (2017). Joint Sequence Analysis: Association and Clustering. \emph{Sociological Methods \& Research}, 46(2), 252–287. \doi{10.1177/0049124115591013}.
}
%\note{ ~~further notes~~ }


\examples{
data(biofam)

## Building one channel per type of event left, children or married
## Using only first 200 sequences
bf <- as.matrix(biofam[1:200, 10:25])
children <-  bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6
## weights
weights <- biofam[1:200,"wp00tbgs"]

## Building sequence objects
child.seq <- seqdef(children)
marr.seq <- seqdef(married)
left.seq <- seqdef(left)
## distances by channel
dchild <- seqdist(child.seq, method="OM", sm="INDELSLOG")
dmarr <- seqdist(marr.seq, method="OM", sm="INDELSLOG")
dleft <- seqdist(left.seq, method="OM", sm="INDELSLOG")
dbiofam <- list(dchild,dmarr,dleft)
names(dbiofam) <- c("child","marr","left")

## Association between domains and R2 by domain
rass <- dissdomassoc(dbiofam, weights=weights)

## Joint distance matrix based on substitution costs
## derived from indels based on log of state frequencies
mcdist <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
                    method="OM", sm =list("INDELSLOG", "INDELSLOG", "INDELSLOG"))

rassj <- dissdomassoc(dbiofam, jointdiss=mcdist, what=c("all"), weights=weights)
rassj[["Pearson.Rsquare"]]
rass[["Pearson.Rsquare"]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\author{Gilbert Ritschard}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity-based analysis}
