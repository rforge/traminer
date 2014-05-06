\name{otsplot}

\alias{otsplot}
\alias{otsplot.default}
\alias{otsplot.formula}

\title{Time-series plot for ordinal longitudinal data}

\description{Plots ordinal sequences of subjects in a \eqn{x} (usually
  time) versus \eqn{y} (response variable) scatterplot. The sequences
  are displayed by jittered frequency-weighted parallel lines.}

\usage{

otsplot(x, ...)

\method{otsplot}{default}(x, y, subject, weights = NULL, cond,
        plot_gp = gpar(), grid_gp = gpar(),
        strip_gp = gpar(),
        filter = NULL, maxit = 200,
        main, xlab, ylab, xlim, ylim,
        layout = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
        pop = TRUE, newpage = TRUE, ...)

}

\arguments{
  \item{x}{a \code{numeric} or \code{factor} vector. Represents usually
    time.} 
  \item{y}{a ordered factor vector (see \code{ordered} in the base
    package).} 
  \item{subject}{a factor vector that allocates elements in \code{x} and
    \code{y} to the observation unit.} 
  \item{weights}{a numeric vector of weights of length equal the number
    of subjects.}
  \item{cond}{a vector (\code{numeric} or \code{factor}) of group
    memberships of length equal the number of subjects. When specified,
    one panel is generated for each different membership value.}
  \item{plot_gp}{graphical parameters for the plot.}
  \item{grid_gp}{graphical parameters for the translation zones.}
  \item{strip_gp}{graphical paramaters for the strips appearing when the
    \code{cond} argument is used.}
  \item{filter}{a list of line coloring options.}
  \item{maxit}{maximal number of iteration for the algorithm that
    establishes the translation arrangement.}
  \item{main}{a main title for the plot.}
  \item{xlab}{label for the \code{x} axis.}
  \item{ylab}{label for the \code{y} axis.}
  \item{xlim}{a label for the \code{x} axis.}
  \item{ylim}{a label for the \code{y} axis.}
  \item{layout}{a integer vector of length 2 specifying the number of
    rows and columns of the panel arrangement when the \code{cond}
    argument is used.}
  \item{margins}{a numeric vector interpreted in the same way as
    the \code{mar} argument of function \code{\link{par}} in base
    graphics.} 
  \item{pop}{a logical whether the viewport tree should be popped before
    return.}
  \item{newpage}{a logical indicating whether \code{grid.newpage()}
    should be called.}
  \item{\ldots}{additional arguments passed to the calls.}
}

\details{
  The function is a scaled down version of the \code{seqpcplot}
  function of the \pkg{TraMineR} package, implemented in the \pkg{grid}
  graphics environment.
}

\author{Reto Buergin}

\references{
  Reto Buergin and Gilbert Ritschard, G. (2012), Categorical parallel
  coordinate plot,
  \emph{In LaCOSA Lausanne Conference On Sequence Analysis, University of Lausanne, June 6th-8th 2012 (Poster)},
  URL \url{http://www3.unil.ch/wpmu/sequences2012/files/2012/05/Buergin_LaCOSA_Poster.pdf}
}

\examples{
## ------------------------------------------------------------------- #
## Example 1: Artificial ordinal longitudinal data
##
## Applying the otsplot on a set of 4 artificial ordinal sequences. 
## ------------------------------------------------------------------- #

## load the data
data(otsplot_1)

## plot the data
otsplot(x = otsplot_1$time, y = otsplot_1$y, subject = otsplot_1$id,
        margins = c(4, 6, 1, 1), xlab = "time")
}

\keyword{hplot}