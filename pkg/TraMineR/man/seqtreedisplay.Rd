\name{seqtreedisplay}
\alias{seqtreedisplay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical rendering of a sequence regression tree}
\description{
Generate a graphical representation of a regression tree of state sequence data.
}
\usage{
seqtreedisplay(tree, filename = NULL, seqdata = tree$info$object,
       imgLeafOnly = FALSE, sortv = NULL,
       dist.matrix = NULL, title.cex = 3,
       withlegend = "auto", legend.fontsize = title.cex,
       axes = FALSE, imageformat = "png",
       withquality = TRUE, quality.fontsize = title.cex, 
       legendtext = NULL, showtree=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A \code{seqtree} object (as produced by \code{\link{seqtree}}).}
  \item{filename}{The name of a file where to save the plot (overwriting existing file). If \code{NULL}, a temporary file is created.}
  \item{seqdata}{The sequence object containing the state sequences plotted in the nodes.}
  \item{imgLeafOnly}{Logical. If \code{TRUE} sequences are plotted only in terminal nodes.}
  \item{sortv}{Argument passed to \code{\link{seqplot}}}
  \item{dist.matrix}{Argument passed to \code{\link{seqplot}}}
  \item{title.cex}{The \code{cex} value for the node titles (see \code{\link{par}}).}
  \item{withlegend}{Logical. Should the color legend be displayed on the plot?}
  \item{legend.fontsize}{Font \code{cex} value for the legend.}
  \item{axes}{Argument passed to \code{\link{seqplot}}}
  \item{imageformat}{Image format of the output file (\code{filename})}
  \item{withquality}{Logical: Should quality measures be displayed?}
  \item{quality.fontsize}{Font size for displaying quality measures.}
  \item{legendtext}{Character. Optional additional legend information text.}
  \item{showtree}{Logical. Should the tree be shown on the screen?}
  \item{\dots}{additional arguments passed to seqplot}
}
\details{
This function generates a tree image. For each node, it invokes \code{\link{seqplot}} for the selected lines of \code{seqdata} as argument.
You should at least specify the \code{type} of the plot to use (\code{type="d"} for instance, see \code{\link{seqplot}} for more details).
	
The plot is actually not generated as an R plot, but with GraphViz (www.graphviz.org). Hence, \code{seqtreedisplay} only works when GraphViz is correctly installed.

Conversion to image formats other than "jpeg" or "png" is done using ImageMagick (www.imagemagick.org). To use this feature,  ImageMagick (www.imagemagick.org) should hence also be installed.
}
\value{None}

\seealso{
See \code{\link{seqtree} for examples}
}
