\name{trmatplot}
\alias{trmatplot}
\alias{trmatplot.default}
\alias{trmatplot.array}
\alias{trmatplot.depmix.fitted}
\alias{trmatplot.hmm}
\title{Transition Matrix Plot}
\description{
A coordinate plot which maps each element in the probability transition matrix as a line, where each line is weighted by probability. Users can apply filters to emphasize the most (or least) probable state sequences overall, or by initial state. Various color palettes using the Hue-Chroma-Luminance color scheme can be easily selected by the user. Input is either an object of class \command{\link[base]{matrix}}, \command{\link[base]{array}}, \command{\link[depmixS4]{depmix.fitted}}, or \code{hmm} which is the output of either \command{\link[seqHMM]{build_hmm}} or \command{\link[seqHMM]{build_mm}}.
}
\usage{
## S3 method for classes 'matrix', 'array', 'depmix.fitted', and 'hmm'

trmatplot ( d, seed = NULL, rowconstraint = TRUE, morder = 1,
cspal = NULL, cpal = NULL, main = NULL, 
xlab = NULL, ylab = NULL, ylim = NULL, xtlab = NULL, ytlab = NULL,
pfilter = NULL, shade.col = "grey80", num = 1, 
hide.col = NULL, lorder = NULL, plot = TRUE, 
verbose = FALSE, title, ...)
}
\arguments{

	\item{d}{Object to be plotted of class \command{\link[base]{matrix}}, \command{\link[base]{array}}, \command{\link[depmixS4]{depmix.fitted}}, or \code{hmm}. A transition matrix of probabilities. See Details.}

	\item{seed}{A single value, interpreted as an integer, or \code{NULL} (default). See Details.}

	\item{rowconstraint}{Logical. Checks if the row constraint is satisfied, i.e. the sum of each row of the transition probability matrix equals one. Default is \code{TRUE}.}

	\item{morder}{Numeric. The order of the probability transition matrix, \code{d}, to be plotted. Default is 1.}

	\item{cspal}{A color palette that can be specified as one of: \code{"dynamic", "harmonic", "cold", "warm", "heat", "terrain"}. The \command{\link[colorspace]{rainbow_hcl}}, \command{\link[colorspace]{heat_hcl}}, and \command{\link[colorspace]{terrain_hcl}} commands are used to generate color palettes.}

	\item{cpal}{Color palette vector when coloring probability sequences. The \command{\link[colorspace]{rainbow_hcl}} command is used to generate  a color palette if none is specified.}

	\item{main}{Title for the graphic. Default is \emph{Probability Transition Matrix}.}

	\item{xlab}{Label for the x axis. Default is \emph{Time}.}

	\item{ylab}{Label for the y axis. Default is \emph{States}.}

	\item{ylim}{Numeric vector of length 2 giving the y coordinates range.}

	\item{xtlab}{Labels for the x axis ticks. Default is time (t, t+1,...).}
	
	\item{ytlab}{Labels for the y axis ticks.}

	\item{pfilter}{Probability filter. Can be specified as one of \code{"tmax", "tmin", "smax", "smin"}. See Details.}

	\item{shade.col}{The color for sequences shaded out using the \code{pfilter} argument. Default is \code{"grey80"}. See Details.}

	\item{num}{Numeric. The number of sequences to be highlighted. Only applicable when using \code{pfilter="tmax"} or \code{pfilter="tmin"}. Default is 1.}

	\item{hide.col}{The color for sequences shaded out using the \code{filter} argument. Default is \code{"grey80"}. See Details.}

	\item{lorder}{Line order. Either \code{"background"} or \code{"foreground"}. When \code{pfilter} is used \code{lorder} is set by default.}

	\item{plot}{Logical. Should the object be plotted. Default is \code{TRUE}.}

	\item{verbose}{Logical. Reports extra information on progress. Default is \code{FALSE}.}

	\item{title}{Deprecated. Use \code{main} instead.}

	\item{...}{Additional arguments, such as graphical parameters, to be passed on. See \command{\link[graphics]{par}} and \command{\link[TraMineR]{seqpcplot}}.}
}
\value{
	\code{trmatplot} returns an object of class \command{\link{trmatplot}}. Some of the arguments are inherent in \command{\link[graphics]{par}} and \command{\link[TraMineR]{seqpcplot}}.
}

\details{

The object \code{d} to be plotted, is a probability matrix which, when of class \command{\link[base]{matrix}} or \command{\link[base]{array}}, can be user defined. Alternatively, if the probability transition matrix is the output packages \pkg{msm}, \pkg{HiddenMarkov} and \pkg{HMM}, it is of class \command{\link[base]{matrix}}. A probability matrix obtained using the package \pkg{depmixS4} is embedded in an object of class \command{\link[depmixS4]{depmix.fitted}}. A probability matrix obtained using the  \pkg{seqHMM} functions \command{\link[seqHMM]{build_hmm}} or \command{\link[seqHMM]{build_mm}} is embedded in an object of class \code{hmm}. Thus, objects of classes \command{\link[depmixS4]{depmix.fitted}} and \code{hmm} are also accepted.

Setting a \code{seed} allows the graphic to be replicated.

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

trmat<-matrix( c (0.1, 0.05, 0.05, 0.80,
                  0.06, 0.02, 0.03, 0.89,
                  0.03, 0.01, 0.01, 0.95,
                  0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

trmatplot(trmat)

#--- Setting a seed so that the graphic can be replicated

trmatplot(trmat, seed = 2)

#--- Adjusting line width

trmatplot(trmat, seed = 2, lwd = 0.8)


###########################################################
# cspal: ready-to-use color palettes using colorspace
###########################################################

#--- Color palette "dynamic"

trmatplot(trmat, seed = 2, cspal = "dynamic")

#--- Color palette "harmonic"

trmatplot(trmat, seed = 2, cspal = "harmonic")

#--- Color palette "cold"

trmatplot(trmat, seed = 2, cspal = "cold")

#--- Color palette "warm"

trmatplot(trmat, seed = 2, cspal = "warm")

#--- Color palette "heat"

trmatplot(trmat, seed = 2, cspal = "heat")

#--- Color palette "terrain"

trmatplot(trmat, seed = 2, cspal = "terrain")

###########################################################
# pfilter: Filtering out most (or least) probable sequences
###########################################################

#--- The most probable sequence given a state

trmatplot(trmat, seed = 2, pfilter = "smax")

#--- The least probable sequence given a state

trmatplot(trmat, seed = 2, pfilter = "smin")

#--- The two most probable sequnces

trmatplot(trmat, seed = 2, pfilter = "tmax", num = 2 )

#--- The ten least probable sequences

trmatplot(trmat, seed = 2, pfilter = "tmin", num = 10 )

######################################################
# filter: Highlighting a specific sequence of interest
######################################################

#--- Highlight the probability that a sequence is initially in state 1 and then in state 4

trmatplot(trmat, seed = 2, filter = list(type = "sequence", value = "(1)-(4)"))
}
\author{
Pauline Poulcheria Adamopoulou
}
\references{
	Buergin, R. and G. Ritschard (2014), "A decorated parallel coordinate plot for categorical longitudinal data", \emph{The American Statistician}. Vol. 68(2), pp. 98-103.

	Zeileis, A.; Hornik, K. and P. Murrell (2009), "Escaping RGBland: Selecting Colors for Statistical Graphics", \emph{Computational Statistics & Data Analysis}. Vol. 53, pp. 3259-3270.
}
\seealso{
	\command{\link{march.Dcmm.trmatplot}}, 
	\command{\link[TraMineR]{seqpcplot}}, 
	\command{\link[graphics]{par}}.
}

