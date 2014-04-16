\name{tvcolmm-plot}

\alias{tvcolmm-plot}
\alias{plot.tvcolmm}

\title{\code{plot} method for \code{tvcolmm} objects}

\description{
  \code{plot} method for \code{tvcolmm} for diagnosis.
}

\usage{
\method{plot}{tvcolmm}(x, type = c("default", "coef", 
        "simple", "terms"),
     main = NULL, drop_terminal = TRUE,
     tnex = 1, newpage = TRUE,
     pop = TRUE, gp = gpar(),
     terminal = FALSE, 
     legend = TRUE, legend_gp = gpar(),
     terms = NULL, variable = NULL, ...)
}

\arguments{
  \item{x}{An object of class \code{\link{tvcolmm}}.}
  \item{type}{the type of the plot. Available types are
    \code{"default"} \code{"simple"} and
    \code{"terms"}.}  
  \item{main}{character. A main title for the plot.}
  \item{drop_terminal}{a logical indicating whether all terminal nodes
    should be plotted at the bottom. See also \code{plot.party}.}
  \item{tnex}{a numeric value giving the terminal node extension in
    relation to the inner nodes.}
  \item{newpage}{a logical indicating whether \code{grid.newpage()}
    should be called.}
  \item{pop}{a logical whether the viewport tree should be popped before
    return.}
  \item{gp}{graphical parameters. See \code{gpar}.}
  \item{terminal}{logical specifying if only plots for node models
    should be shown.}
  \item{legend}{logical specifying if a legend should be
    added. Currently not supported by any \code{type}.}
  \item{legend_gp}{graphical parameters for the legend. Currently not
    supported by any \code{type}.}
  \item{terms}{character string. Which coefficients to be
    visualized in the \code{"terms"}-plot.}
  \item{variable}{character string. Which partitioning covariates to
    be visualized in the \code{"terms"} plot.}
  \item{\ldots}{Additional arguments passed to the calls.}
}

\details{
  The available plot functions allow the inspection of fitted
  \code{\link{tvcolmm}} objects. \code{"default"} and \code{"coef"} show
  the tree structure and coefficients in the node. \code{"fluctest"}
  shows the fluctuation tests in the nodes and terminal
  nodes. \code{"simple"} shows the tree structure only. \code{"terms"}
  plots the partial relationships between partitioning variables and
  model coefficients. The last type may be very slow.
}

\seealso{\code{\link{tvcolmm}}, \code{\link{tvcolmm-methods}}}

\author{Reto Buergin}

\keyword{hplot}
