\name{trmatplot}
\alias{trmatplot}
\title{Transition Matrix Plot}
\usage{
trmatplot(d, seed = NULL, cpal = NULL, title = NULL, xlab = NULL,
  ylab = NULL, ylim = NULL, xtlab = NULL, pfilter = NULL,
  shade.col = "grey80", num = 1, hide.col = NULL, lorder = NULL,
  verbose = FALSE, ...)
}
\arguments{
\item{d}{Object to be plotted. A transition matrix of probabilities.}

\item{seed}{A single value, interpreted as an integer, or \code{NULL} (default). See Details.}

\item{cpal}{Color palette vector when coloring probability sequences. The \command{\link[colorspace]{rainbow_hcl}} command is used to generate  a color palette if none is specified.}

\item{title}{Title for the graphic. Default is \emph{Probability Transition Matrix}.}

\item{xlab}{Label for the x axis. Default is \emph{Time}.}

\item{ylab}{Label for the y axis. Default is \emph{States}.}

\item{ylim}{Numeric vector of length 2 giving the y coordinates range.}

\item{xtlab}{Label for the x axis ticks. Default is time (t, t+1,...).}

\item{pfilter}{Probability filter. Can be specified as one of \code{"tmax", "tmin", "smax", "smin"}. See Details.}

\item{shade.col}{The color for sequences shaded out using the \code{pfilter} argument. Default is \code{"grey80"}. See Details.}

\item{num}{Numeric. The number of sequences to be highlighted. Only applicable when using \code{pfilter="tmax"} or \code{pfilter="tmin"}. Default is 1.}

\item{hide.col}{The color for sequences shaded out using the \code{filter} argument. Default is \code{"grey80"}. See Details.}

\item{lorder}{Line order. Either \code{"background"} or \code{"foreground"}. When \code{pfilter} is used \code{lorder} is set by default.}

\item{verbose}{Logical. Reports extra information on progress. Default is \code{FALSE}.}

\item{...}{Additional arguments, such as graphical parameters, to be passed on. See \command{\link[graphics]{par}} and \command{\link[TraMineR]{seqpcplot}}.}
}
\value{
\code{trmatplot} returns an object of class \command{\link[TraMineR]{seqpcplot}}
}

\details{
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
\description{
A coordinate plot which maps each row in the probability transition matrix as a line, where each line is weighted by probability. Users can apply filters to empasize the most (or least) probable state sequences overall, or by initial state.
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

###########################################################
# pfilter: Filtering out most (or least) probable sequences
###########################################################

#--- The most probable sequence given a state

trmatplot(trmat, seed = 2, pfilter = "smax")

#--- The least probable sequence given a state

trmatplot(trmat, seed = 2, pfilter = "smin")

#--- The four most probable sequnces

trmatplot(trmat, seed = 2, pfilter = "tmax", num = 4 )

#--- The ten least probable sequences

trmatplot(trmat, seed = 2, pfilter = "tmin", num = 10 )

######################################################
# filter: Highlighting a specific sequence of interest
######################################################

trmatplot(trmat, seed = 2, filter = list(type = "sequence", value = "(1)-(4)"))
}
\author{
Poulcheria Adamopoulou
}
\references{
Buergin, R. and G. Ritschard (2014), "A decorated parallel coordinate plot for categorical longitudinal data", \emph{The American Statistician}. Vol. 68(2), pp. 98-103.
}
\seealso{
\command{\link[TraMineR]{seqpcplot}}
}

