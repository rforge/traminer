\name{plot.seqdiff}
\alias{plot.seqdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting a seqdiff object.
}
\description{
  Plot method for the sliding values returned by \code{\link{seqdiff}}. Plots a statistic (the Pseudo R2 by default) along the position axis.
}
\usage{
\method{plot}{seqdiff}(x, stat = "Pseudo R2", type = "l", ylab = stat, xlab = "",
			legendposition = "top", ylim = NULL, xaxt = TRUE, col = NULL, xtstep=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object produced by \code{\link{seqdiff}}}
  \item{stat}{Character. Name of the statistic to be plotted. Can be any of the statistics returned by \code{\link{seqdiff}} or \code{"discrepancy"}. See details.}
  \item{type}{the line \code{type}, see \code{\link{lines}} }
  \item{ylab}{Character. y-axis label.}
  \item{xlab}{Character. x-axis label.}
  \item{legendposition}{Character. Position of the line legend, see \code{\link{legend}}}
  \item{ylim}{Numeric. If not \code{NULL}, range of the y-axis.}
  \item{xaxt}{Logical. If \code{TRUE} an x-axis is plotted.}
  \item{col}{list of colors to use for each line.}
  \item{xtstep}{optional interval at which the tick-marks and labels of the x-axis are displayed. If unspecified, the \code{xtstep} attribute of the sequence object used in \code{\link{seqdiff}} is used.}
  \item{\dots}{Additional parameters passed to \code{\link{lines}}}
}
\details{
  The function plots the sliding values of the requested statistic.

  You can plot the evolution of two statistics by providing for instance \code{stat=c("Pseudo R2", "Levene")}.

  Use \code{stat="discrepancy"} to plot the within-discrepancies.

  For \code{"discrepancy"}, a separate line is drawn for the whole set of sequences and for each group. Those two values cannot be paired with another statistic.
}

\seealso{\code{\link{seqdiff}}} 