\name{seqtreedisplay}
\alias{seqtreedisplay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical rendering of a sequence regression tree}
\description{
Generate a graphical representation of a regression tree of state sequence data.
}
\usage{
seqtreedisplay(tree, filename = NULL, seqdata = tree$info$object, imgLeafOnly = FALSE,
        sortv = NULL, dist.matrix = NULL, title.cex = 3, withlegend = "auto",
        legend.fontsize = title.cex, axes = FALSE, imageformat = "png",
        withquality = TRUE, legendtext = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree to be plotted.}
  \item{filename}{The name of a file where to save the plot. If \code{NULL}, a temporary file is created.}
  \item{seqdata}{The sequence object containing the state sequences plotted in the nodes.}
  \item{imgLeafOnly}{Logical. If \code{TRUE} sequences are plotted only in terminal nodes.}
  \item{sortv}{Argument passed to \code{\link{seqplot}}}
  \item{dist.matrix}{Argument passed to \code{\link{seqplot}}}
  \item{title.cex}{The \code{cex} value for the node titles (see \code{\link{par}}).}
  \item{withlegend}{If \code{TRUE}, a node with the color legend is added on the plot.}
  \item{legend.fontsize}{Font \code{cex} value for the legend.}
  \item{axes}{Argument passed to \code{\link{seqplot}}}
  \item{imageformat}{Image format of the output file (\code{filename})}
  \item{withquality}{If \code{TRUE}, a node displaying fitting measures of the tree is added to the plot.}
  \item{legendtext}{Character. Optional text information that should be added.}
  \item{\dots}{additional arguments passed to seqplot}
}
\details{
The plot is actually not generated as an R plot but with GraphViz (www.graphviz.org). Hence, \code{seqtreedisplay} will only work when GraphViz is correctly installed.

Conversion to image formats other than "jpeg" or "png" is done using ImageMagick (www.imagemagick.org). To use this feature,  ImageMagick should hence also be installed.
}
\value{None}

\seealso{
See \code{\link{seqtree}}
}
