\name{march.Dcmm.trmatplot}
\alias{march.Dcmm.trmatplot}
\title{Transition Matrix Plot for march.Dcmm objects}
\author{
Poulcheria Adamopoulou
}
\description{
A coordinate plot which maps each element in the probability transition matrix as a line, where each line is weighted by probability. Users can apply filters to emphasize the most (or least) probable state sequences overall, or by initial state. Various color palettes using the Hue-Chroma-Luminance color scheme can be easily selected by the user. Input is an object of class \code{march.Dcmm} which is the output of \code{march.dcmm.construct}.
}

\usage{
## For class 'march.Dcmm'

march.Dcmm.trmatplot ( d, seed = NULL, type = "hidden", hstate = 1,
cspal = NULL, cpal = NULL, title = NULL, xlab = NULL,
ylab = NULL, ylim = NULL, xtlab = NULL, ytlab = NULL, 
pfilter = NULL, shade.col = "grey80", num = 1, hide.col = NULL, 
lorder = NULL, plot = TRUE, verbose = FALSE, ...)
}
\arguments{
	\item{d}{Object to be plotted. A \code{march.Dcmm} object.}

	\item{seed}{A single value, interpreted as an integer, or \code{NULL} (default). See Details.}

	\item{type}{Character string. Can be specified as either \code{"hidden"}, if the hidden matrix is to be plotted (default) or as \code{"visible"} if the visible matrix is to be plotted.}

	\item{hstate}{Numeric. Valid when \code{type = "visible"}. Specifies from which hidden state the (visible) probability transition matrix should be plotted. Default \code{hstate = 1}.}

	\item{cspal}{A color palette that can be specified as one of: \code{"dynamic", "harmonic", "cold", "warm", "heat", "terrain"}. The \command{\link[colorspace]{rainbow_hcl}}, \command{\link[colorspace]{heat_hcl}}, and \command{\link[colorspace]{terrain_hcl}} commands are used to generate color palettes. See Examples in \command{\link{trmatplot}}.
}

	\item{cpal}{Color palette vector when coloring probability sequences. The \command{\link[colorspace]{rainbow_hcl}} command is used to generate  a color palette if none is specified.}

	\item{title}{Title for the graphic. Default is \emph{Probability Transition Matrix}.}

	\item{xlab}{Label for the x axis. Default is \emph{Time}.}

	\item{ylab}{Label for the y axis. Default is \emph{States}.}

	\item{ylim}{Numeric vector of length 2 giving the y coordinates range.}

	\item{xtlab}{Label for the x axis ticks. Default is time (t, t+1,...).}

	\item{ytlab}{Labels for the y axis ticks.}

	\item{pfilter}{Probability filter. Can be specified as one of \code{"tmax", "tmin", "smax", "smin"}. See Details.}

	\item{shade.col}{The color for sequences shaded out using the \code{pfilter} argument. Default is \code{"grey80"}. See Details.}

	\item{num}{Numeric. The number of sequences to be highlighted. Only applicable when using \code{pfilter="tmax"} or \code{pfilter="tmin"}. Default is 1.}

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
Setting a \code{seed} allows the graphic to be replicated.

The \code{pfilter} argument serves to highlight probability sequences that are either most probable while shading out those that are less probable in \code{shade.col} and vice-versa.
The four options for \code{pfilter} are described below, and are illustrated in Examples in \command{\link{trmatplot}}.
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
 See Examples in \command{\link{trmatplot}}.
}

\references{
	Buergin, R. and G. Ritschard (2014), "A decorated parallel coordinate plot for categorical longitudinal data", \emph{The American Statistician}. Vol. 68(2), pp. 98-103.

	Zeileis, A.; Hornik, K. and P. Murrell (2009), "Escaping RGBland: Selecting Colors for Statistical Graphics", \emph{Computational Statistics & Data Analysis}. Vol. 53, pp. 3259-3270.
}
\seealso{
	\command{\link{trmatplot}},  
	\command{\link{trmatplot.default}},  
	\command{\link{trmatplot.depmix.fitted}}, 
	\command{\link{trmatplot.array}}, 
	\command{\link[TraMineR]{seqpcplot}}, 
	\command{\link[graphics]{par}}.
}

