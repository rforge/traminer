\encoding{latin1}
\name{seqdistOO}
\alias{seqdistOO}
\title{Distances (dissimilarities) between sequences}

\description{
Computes pairwise dissimilarities between sequences or
dissimilarities with a reference sequence. Several dissimilarities measures can be chosen including optimal matching (OM) and many of its variants, distance based on the count of common attributes, and distances between within sequence state distributions.
 }

\usage{
seqdistOO(seqdata, method, refseq=NULL, norm=FALSE,
        indel=1, sm=NA, with.missing=FALSE, full.matrix=TRUE,
        kweights=rep(1.0, ncol(seqdata)), tpow=1,
        expcost=0.5, context=1-2*expcost,
        link="mean", h=0.5, nu=0,
        transindel="constant", otto, previous=FALSE,
        addcolumn=TRUE, numcpu=1,
        breaks=NULL, step=1, overlap=FALSE, weighted=TRUE)
}

\arguments{
  \item{seqdata}{State sequence object defined with the \code{\link{seqdef}} function.}

  \item{method}{Character string indicating the metric to be used. One of
  \code{"OM", "OMloc", "OMslen", "OMspell", "OMstran", "LCS", "HAM", "DHD", "TWED", "SVRspell", "NMS", "NMSMST", "LCP", "RLCP", "CHI2", "EUCLID"} (see details)
}

  \item{refseq}{Optional baseline sequence to compute the distances from. Can be the index of a sequence in the state sequence object, 0 for the most frequent sequence, or an external sequence passed as a sequence object with a single row. %The latter works only with a subset of methods (see details). 
  }

  \item{norm}{String or logical. One of \code{"none", "maxlength", "gmean", "maxdist", "YujianBo"}. If \code{TRUE} the normalization method is automatically selected (see details).  Default is \code{FALSE}.}

  \item{indel}{Scalar or numeric vector. The insertion/deletion cost (for all OM variants including TWED). Default is 1. Ignored with non OM metrics. If it is a vector, the indel costs are state dependant (the vector length should be equal to the number of states).}

  \item{sm}{substitution-cost matrix (for OM variants, TWED, SVRspell, HAM and DHD method). Can also be one of the \code{\link{seqcost}} methods \code{"TRATE"} or \code{"CONSTANT"} for building a \code{sm} matrix. Default is \code{NA}.}

  \item{with.missing}{Logical. Must be \code{TRUE} when sequences contain non deleted gaps (missing values). See details.}

  \item{full.matrix}{Logical. If \code{TRUE} (default), the full distance matrix is returned. This is for compatibility with  earlier versions of the \code{seqdist} function. If \code{FALSE}, an object of class \code{\link{dist}} is returned, that is, a vector containing only values from the upper triangle of the distance matrix. Objects of class \code{dist} are smaller and can be passed directly as arguments to most clustering functions. Ignored when \code{refseq} is set.}

  \item{kweights}{Numerical vector with at position \eqn{k} the weight applied to subsequences of length \eqn{k} (for all NMS distances).}

  \item{expcost}{Scalar. Cost of spell length transformation (the exact interpretation is distance dependant). For \code{OMspell} and \code{OMloc}.}

  \item{context}{Scalar. Cost of local insertion (for \code{OMloc} only).}

  \item{link}{Function name. Function used to compute some internal costs.\cr
  For \code{OMloc}, the function used in the dynamic computation of insertion/deletion costs; one of \code{"mean"} (default), \code{"min"} (minimum) or \code{"previous"} (previous state only). \cr
   For \code{"OMslen"}, the function used to compute substitution costs; one of \code{"mean"} (arithmetic average) or \code{"gmean"} (geometric mean as in the original proposition of Halpin 2010).}
   %\cr
  %For \code{NMS, NMSMST, SVRspell}, if \code{link=="log"}, the logarithm of common attributes is computed before computing the final distance.}

  \item{h}{Scalar. Exponential weight of spell length (for \code{OMslen}).}

  \item{tpow}{Scalar. Exponential weight of spell length (for \code{OMspell} and \code{SVRspell}).}

  \item{transindel}{Character string. Method for computing transition indel costs. One of \code{"constant"} (single indel of 1), \code{"subcost"} (based on substitution costs), \code{"prob"} (based on transition probabilities) (for \code{OMstran}).}

  \item{otto}{Scalar. Origin-transition trade-off weight  (for \code{OMstran}).}

  \item{previous}{if \code{TRUE}, add transition from previous state (for \code{OMstran}).}

  \item{addcolumn}{if \code{TRUE}, repeat the first (if \code{previous=TRUE}) and the last column (for \code{OMstran}).}

  \item{breaks}{An optional list containing the spells used to build the intervals (for \code{CHI2}).}

  \item{step}{The length of the interval (only used if \code{breaks} is not \code{NULL}) (for \code{CHI2}).}

  \item{weighted}{Logical. Should distributions account for the sequence weights in the state sequence object? (see \code{\link{seqdef}}.) Set as \code{FALSE} if you want to ignore the weights (for \code{CHI2, OMstran}).}

  \item{overlap}{Logical. Are overlapping spells allowed (for \code{CHI2, EUCLID}).}

  \item{nu}{Scalar. Stiffness (for \code{TWED}). See \cite{Halpin (2014), p 88.}}

  \item{numcpu}{Scalar. Number of CPU to be used.}
}

\details{
The \code{seqdistOO} function returns a matrix of distances between sequences or a vector of distances to the reference sequence when \code{refseq} is set. The available metrics (see 'method' option) include

\itemize{
 \item{\emph{Edit distances}: optimal matching (\code{OM}), localized OM (\code{OMloc}), spell length sensitive OM (\code{OMslen}), OM of spells (\code{OMspell}), OM of sequences of transitions (\code{OMstran}), Hamming (\code{HAM}), dynamic Hamming (\code{DHD}), and the time warp edit distance (\code{TWED}).}

 \item{\emph{Metrics based on counts of common attributes}: distance based on the longest common subsequence (\code{LCS}), on the longest common prefix (\code{LCP}), on the longest common suffix (\code{RLCP}), on the number of matching subsequences (\code{NMS}), and the subsequence vectorial representation metric (\code{SVRspell}).}

 \item{\emph{Distances between state distributions}: Euclidean (\code{EUCLID}), Chi-squared (\code{CHI2}).
 %, and the Chi-squared distance between vectors of transition rates to the future states (\code{CHI2fut}).
 }
}

See \cite{Studer and Ritschard (2014)} for a description and the comparison of the above dissimilarity measures except \code{TWED} for which we refer to \cite{Marteau (2009)} and \cite{Halpin (2014)}.

Except for a few dissimilarity measures (\code{LCS, LCP, RLCP, NMS}) that do not require parameters, each method can be controlled with some parameters as described in the following table.

\tabular{ll}{
  method \tab parameters \cr
  ------------------ \tab ---------------------------------\cr
  \verb{OM} \tab \verb{sm, indel} \cr
  \verb{OMloc} \tab \verb{sm, expcost, context} \cr
  \verb{OMslen} \tab \verb{sm, indel, h} \cr
  \verb{OMspell} \tab \verb{sm, indel, expcost} \cr
  \verb{OMstran} \tab \verb{sm, indel, otto, previous, addcolumn, transindel} \cr
  \verb{TWED} \tab \verb{sm, indel, h, nu} \cr
  \verb{HAM, DHD} \tab  \verb{sm}\cr
  \verb{EUCLID, CHI2} \tab  \verb{step, breaks, overlap}\cr
  \verb{NMS, NMSMST} \tab \verb{kweights}\cr
  \verb{SVRspell} \tab \verb{sm, kweights, tpow}\cr
------------------ \tab ---------------------------------
  }

In addition, on multiprocessor computers, one can specify the number \code{numcpu} of cpu's to be used.

 %longest common prefix (\code{"LCP"}), longest common suffix (\code{"RLCP"}), longest common subsequence (\code{"LCS"}), Hamming distance (\code{"HAM"}) and Dynamic Hamming Distance (\code{"DHD"}). The Hamming distance is OM without indels and the
%Dynamic Hamming Distance is HAM with specific substitution costs at each position as proposed by \cite{Lesnard (2006)}. Note that HAM and DHD apply only to sequences of equal length.


%For OM, HAM and DHD, a user specified substitution cost matrix can be provided with the \code{sm} argument.

For DHD, when used \code{sm} should be a series of matrices grouped in a 3-dimensional matrix with the third index referring to the position in the sequence.

When \code{sm} is not specified, a single cost of 2 is used for all substitutions, except for HAM for which it is 1, and for DHD for which the costs are derived from the transition rates at the successive positions.

Distances can optionally be normalized by means of the \code{norm} argument. If set to \code{TRUE}, Elzinga's normalization (similarity divided by geometrical mean of the two sequence lengths) is applied to LCP, RLCP and LCS distances, while Abbott's normalization (distance divided by length of the longer sequence) is used for OM, HAM and DHD. Elzinga's method can be forced with \code{"gmean"} and Abbott's rule with \code{"maxlength"}. With \code{"maxdist"} the distance is normalized by its maximal possible value. For more details, see \cite{Gabadinho et al. (2009, 2011)}. Finally, \code{"YujianBo"} is the normalization proposed by \cite{Yujian and Bo (2007)} that preserves the triangle inequality.

%Distances to a provided reference sequence (specified with \code{refseq}) can currently not be computed with \code{CHI2}.

When sequences contain gaps and the \code{gaps=NA} option was passed to \code{\link{seqdef}}---i.e. when there are non deleted missing values---the \code{with.missing} argument should be set as \code{TRUE}. If left as \code{FALSE} the function stops when it encounters a gap. This is to make the user aware that there are gaps in his sequences. For methods that need an \code{sm} value, \code{seqdist} expects a substitution cost matrix with a row and a column entry for the missing state (symbol defined with the \code{nr} option of \code{\link{seqdef}}). Substitution cost matrices returned by \code{\link{seqcost}} and \code{\link{seqsubm}} will include these additional entries when the functions are called with \code{with.missing=FALSE}. More details on how to compute distances with sequences containing gaps can be found in \cite{Gabadinho et al. (2009)}.
}

\value{
 When \code{refseq} is \code{NULL} (default), the whole matrix (or a \code{dist} object when \code{full.matrix=FALSE}) of pairwise distances between sequences is returned. Otherwise a vector with distances from the sequences in the data sequence object to the reference sequence specified with \code{refseq} is returned.
 }

\seealso{
 \code{\link{seqsubm}}, \code{\link{seqdef}}, and for multichannel distances \code{\link{seqdistmc}}, \code{\link{seqcost}}.
}

\references{
Studer, M. and G. Ritschard (2014). "A Comparative Review of Sequence Dissimilarity Measures". \emph{LIVES Working Papers}, \bold{33}. NCCR LIVES, Switzerland. DOI: \url{http://dx.doi.org/10.12682/lives.2296-1658.2014.33}

%Elzinga, C. H. (2008). Sequence analysis: Metric representations of categorical time
%series. \emph{Technical Report}, Department of Social Science Research Methods, Vrije Universiteit, Amsterdam.

%Elzinga, C. H. and Studer, M. (2014). Spell sequences, state proximities and distance
%metrics. \emph{Sociological Methods and Research}. First published online. DOI: http://dx.doi.org/10.1177/0049124114540707

Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1--37.

Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide. Department of Econometrics and Laboratory of Demography, University of Geneva

Halpin, B. (2014). Three Narratives of Sequence Analysis, in Blanchard, P., Bühlmann, F. and Gauthier, J.-A. (Eds.) \emph{Advances in Sequence Analysis: Theory, Method, Applications}, Vol 2 of Series \emph{Life Course Research and Social Policies}, pages 75--103, Heidelberg: Springer. DOI: \url{http://dx.doi.org/10.1007/978-3-319-04969-4_5}

%Lesnard, L. (2006). Optimal Matching and Social Sciences. \emph{Série des Documents de Travail du CREST},  Institut National de la Statistique et des Etudes Economiques, Paris.

Marteau, P.-F. (2009). Time Warp Edit Distances with Stiffness Adjustment for Time Series Matching. \emph{IEEE Transactions on Pattern Analysis and Machine Learning Intelligence}, \bold{31}(2), 306--318. DOI: \url{http://dx.doi.org/10.1109/TPAMI.2008.76}

Yujian, L. and Bo, L. (2007). A normalized Levenshtein distance metric. \emph{IEEE Transactions
On Pattern Analysis And Machine Intelligence}, \bold{29}(6), 1091--1095. DOI: \url{http://dx.doi.org/10.1109/TPAMI.2007.1078}


%Rousset, P., Giret, J.-F., and Grelet, Y. (2012). Typologies de parcours et dynamique longitudinale.
%\emph{Bulletin of Sociological Methodology/Bulletin de Méthodologie Sociologique},
%114(1):5–34.

See also all references in \cite{Studer and Ritschard (2014)}
}

\examples{
## optimal matching distances with substitution cost matrix
## derived from transition rates
data(biofam)
biofam.seq <- seqdef(biofam, 10:25)
costs <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdistOO(biofam.seq, method="OM", indel=3, sm=costs)

## =====================
## Example with missings
## =====================
data(ex1)
ex1.seq <- seqdef(ex1,1:13)

subm <- seqsubm(ex1.seq, method="TRATE", with.missing=TRUE)
ex1.om <- seqdistOO(ex1.seq, method="OM", sm=subm, with.missing=TRUE)
}
\keyword{misc}
