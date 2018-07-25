\name{trmatplot}
\alias{trmatplot}
\alias{trmatplot.default}
\alias{trmatplot.array}
\alias{trmatplot.depmix.fitted}
\alias{trmatplot.dthmm}
\alias{trmatplot.hmm}
\alias{trmatplot.msm}
\title{Transition Matrix Plot}
\description{
A parallel coordinate plot which maps each element in the probability transition matrix (of any order) as a line, where each line is weighted by probability. Users can apply filters to emphasize the most (or least) probable state sequences overall, or by initial state. Various color palettes using the Hue-Chroma-Luminance color scheme can be easily selected by the user. Input is either an object of class \command{\link[base]{matrix}}, \command{\link[base]{array}}, \command{\link[depmixS4]{depmix.fitted}}, or  a \command{\link[base]{list}} of class \command{\link[msm]{msm}}, \code{hmm} or of class \code{dthmm}.
}
\usage{
## S3 method for classes 'matrix', 'array', 
## 'depmix.fitted', 'dthmm', 'hmm', and 'msm'

trmatplot ( d, seed = NULL, rowconstraint = TRUE, morder = 1,
cspal = NULL, cpal = NULL, main = NULL, 
xlab = NULL, ylab = NULL, ylim = NULL, xtlab = NULL, ytlab = NULL,
pfilter = NULL, shade.col = "grey80", num = NULL, 
hide.col = NULL, lorder = NULL, plot = TRUE, 
verbose = FALSE, ...)
}
\arguments{

	\item{d}{Object to be plotted. A transition matrix of probabilities. Can be of any order. An object of class \command{\link[base]{matrix}}, \command{\link[base]{array}}, \command{\link[depmixS4]{depmix.fitted}}, \code{dthmm}, \code{hmm}, or \command{\link[msm]{msm}}. See Details.}

	\item{seed}{A single value, interpreted as an integer, or \code{NULL} (default). See Details.}

	\item{rowconstraint}{Logical. Checks if the row constraint is satisfied, i.e. the sum of each row of the transition probability matrix equals one. Default is \code{TRUE}.}

	\item{morder}{Numeric. The order of the probability transition matrix, \code{d}, to be plotted. Default is 1. See Details.}

	\item{cspal}{A color palette that can be specified as one of: \code{"dynamic"}, \code{"harmonic"}, \code{"cold"}, \code{"warm"}, \code{"heat"}, or \code{"terrain"}. The \command{\link[colorspace]{rainbow_hcl}}, \command{\link[colorspace]{heat_hcl}}, and \command{\link[colorspace]{terrain_hcl}} commands are used to generate color palettes. The \command{\link[colorspace]{rainbow_hcl}} command is used to generate  a color palette if none is specified.}

	\item{cpal}{Color palette vector when coloring probability sequences. Cannot specify both \code{cspal} and \code{cpal}.}

	\item{main}{Title for the graphic. Default is \emph{Probability Transition Matrix}.}

	\item{xlab}{Label for the x axis. Default is \emph{Time}.}

	\item{ylab}{Label for the y axis. Default is \emph{States}.}

	\item{ylim}{Numeric vector of length 2 giving the y coordinates range.}

	\item{xtlab}{Labels for the x axis ticks. Default is time (t, t+1,...).}
	
	\item{ytlab}{Labels for the y axis ticks.}

	\item{pfilter}{Probability filter. Can be specified as one of \code{"tmax"}, \code{"tmin"}, \code{"smax"}, \code{"smin"}. See Details.}

	\item{shade.col}{The color for sequences shaded out using the \code{pfilter} argument. Default is \code{"grey80"}. See Details.}

	\item{num}{Numeric. The number of sequences to be highlighted, a whole number from 1 to \code{length(d) - 1}. Must be specified only  when using \code{pfilter="tmax"} or \code{pfilter="tmin"}. Default is NULL.}

	\item{hide.col}{The color for sequences shaded out using the \code{filter} argument. Default is \code{"grey80"}. See Details.}

	\item{lorder}{Line order. Either \code{"background"} or \code{"foreground"}. When \code{pfilter} is used \code{lorder} is set by default.}

	\item{plot}{Logical. Should the object be plotted. Default is \code{TRUE}.}

	\item{verbose}{Logical. Reports extra information on progress. Default is \code{FALSE}.}

	\item{...}{Additional arguments, such as graphical parameters, to be passed on. See \command{\link[graphics]{par}} and \command{\link[TraMineR]{seqpcplot}}.}
}
\value{
	\code{trmatplot} returns an object of class \command{\link{trmatplot}}. Some of the arguments are inherent in \command{\link[graphics]{par}} and \command{\link[TraMineR]{seqpcplot}}.
}

\details{
The object \code{d} to be plotted, is a probability transition matrix that can be of class \command{\link[base]{matrix}}, \command{\link[base]{array}}, \command{\link[depmixS4]{depmix.fitted}}, \code{dthmm}, \code{hmm}, \code{msm}.
If the probability transition matrix is the output of packages implemented using \code{S3} methods and classes, such as \pkg{msm}, \pkg{HiddenMarkov},  \pkg{HMM}, or \pkg{seqHMM} it is a \code{matrix} which is element of a \command{\link[base]{list}}. 
More specifically, using the package \pkg{HMM}, the function \command{\link[HMM]{initHMM}} returns a \command{\link[base]{list}} containing the element \code{transProbs} with the probability transition matrix.
Within the \pkg{HiddenMarkov}, \command{\link[HiddenMarkov]{dthmm}} returns a \command{\link[base]{list}} of class \code{dthmm} with the element \code{Pi} which is a probability transition matrix, either user-defined or estimated. 
Similarly, the \pkg{seqHMM} functions \command{\link[seqHMM]{build_mm}} and \command{\link[seqHMM]{build_hmm}} return a \command{\link[base]{list}} of class \code{hmm} containing the element \code{transition_probs} with the estimated probability transition matrix.
The package \pkg{msm} proposes the \command{\link[msm]{pmatrix.msm}} and \command{\link[msm]{pmatrix.piecewise.msm}} to extract the probability transition matrix from a fitted multi-state model (a \command{\link[base]{list}} of class \code{msm}), as returned by \command{\link[msm]{msm}}.
In the case of \pkg{depmixS4}, which uses \code{S4} classes, a probability transition matrix may be obtained depending on the type of model computed using \command{\link[depmixS4]{depmix}} followed by \command{\link[depmixS4]{fit}} and is embedded in an object of class \command{\link[depmixS4]{depmix.fitted}}.

Setting a \code{seed} allows the graphic to be replicated.

If \code{morder > 1}, in other words, the order of the probability transition matrix \code{d} is of order greater than one, then \code{d} must be specified in reduced form. Structural zeroes are not accepted.

The \code{pfilter} argument serves to highlight probability sequences that are either most probable while shading out those that are less probable in \code{shade.col} and vice-versa.
The four options for \code{pfilter} are described below, and are illustrated in Examples.
\describe{
	\item{\code{"smax"}}{For each initial state the most probable next state is highlighed.}
	\item{\code{"smin"}}{For each initial state the least probable next state is highlighed.}
	\item{\code{"tmax"}}{The sequence of states with the highest probability overall is highlighed. To highlight the \emph{n} most probable sequences of states, set \code{num = n}.}
	\item{\code{"tmin"}}{The sequence of states with the lowest probability overall is highlighed. To highlight the \emph{n} least probable sequences of states, set \code{num = n}.}
}

 The \code{filter} and \code{hide.col} arguments are inherent in and may be passed on to \command{\link[TraMineR]{seqpcplot}}.
 The \code{filter} argument serves to specify filters to gray less interesting patterns.
 The filtered-out patterns are displayed in the \code{hide.col} color.
 The filter argument expects a list with at least elements \code{type} and \code{value}.
 Most relevant within the context of probabilities is \code{type = "sequence"}, which highlights the specific pattern. 
 % for example assign \preformatted{filter = list(type = "sequence", value = "(1)-(4)")}
 See Examples.
}


\examples{
##########################################
# Plotting a probability transition matrix
##########################################

trMat<-matrix(c(0.1, 0.05, 0.05, 0.80,
  0.06, 0.02, 0.03, 0.89,
  0.03, 0.01, 0.01, 0.95,
  0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

trmatplot(trMat)

#--- Setting a seed so that the graphic can be replicated

trmatplot(trMat, seed = 2)

#--- Adjusting line width

trmatplot(trMat, seed = 2, lwd = 0.8)

#--- Defining a second order probability transition matrix as an array

trMatArray <- array(c(
  0.30, 0.70,
  0.65, 0.35, 
  0.05, 0.95,
  0.99, 0.01), dim = c(1, 2, 4))

#--- Plotting with user-defined colors

trmatplot(trMatArray, seed = 3, morder = 2, cpal = c("grey40", "grey70"))

###########################################################
# cspal: ready-to-use color palettes using colorspace
###########################################################

#--- Color palette "dynamic"

trmatplot(trMat, seed = 2, cspal = "dynamic")

#--- Color palette "harmonic"

trmatplot(trMat, seed = 2, cspal = "harmonic")

#--- Color palette "cold"

trmatplot(trMat, seed = 2, cspal = "cold")

#--- Color palette "warm"

trmatplot(trMat, seed = 2, cspal = "warm")

#--- Color palette "heat"

trmatplot(trMat, seed = 2, cspal = "heat")

#--- Color palette "terrain"

trmatplot(trMat, seed = 2, cspal = "terrain")

###########################################################
# pfilter: Filtering out most (or least) probable sequences
###########################################################

#--- The most probable sequence given a state

trmatplot(trMat, seed = 2, pfilter = "smax")

#--- The least probable sequence given a state

trmatplot(trMat, seed = 2, pfilter = "smin")

#--- The two most probable sequnces

trmatplot(trMat, seed = 2, pfilter = "tmax", num = 2 )

#--- The ten least probable sequences

trmatplot(trMat, seed = 2, pfilter = "tmin", num = 10 )

######################################################
# filter: Highlighting a specific sequence of interest
######################################################

#--- Highlight the probability that a sequence is initially in state 1 and then in state 4

trmatplot(trMat, seed = 2, filter = list(type = "sequence", value = "(1)-(4)"))
}
\author{
Pauline Adamopoulou
}
\references{
	Buergin, R. and G. Ritschard (2014), "A decorated parallel coordinate plot for categorical longitudinal data", \emph{The American Statistician}. Vol. 68(2), pp. 98-103.

	Zeileis, A.; Hornik, K. and P. Murrell (2009), "Escaping RGBland: Selecting Colors for Statistical Graphics", \emph{Computational Statistics & Data Analysis}. Vol. 53, pp. 3259-3270.
}
\seealso{
	\command{\link[TraMineR]{seqpcplot}}, 
	\command{\link[graphics]{par}}.
}

