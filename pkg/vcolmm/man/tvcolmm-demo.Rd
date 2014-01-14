\name{vcolmm-demo}

\alias{otsplot_1}
\alias{tvcolmm_1}
\alias{tvcolmm_2}

\docType{data}

\title{Synthetic data sets}

\description{Synthetic data sets for illustrations on
  \code{\link{tvcolmm}} and \code{\link{otsplot}}.} 

\usage{
data(otsplot_1)
data(tvcolmm_1)
data(tvcolmm_2)
}

\format{
  \describe{
    \item{\code{y}}{ordered factor. The response variable}
    \item{\code{id}}{factor. The subject identification vector.} 
    \item{\code{time, treat}}{predictor variables.}
    \item{\code{z1, z2, z3, z4, z5,
	z6}}{partitioning/ varying-coefficients variables.}
  }
}

\details{\code{otsplot_1} is a data.frame with 20 observations on 5 subjects
  at 5 occasions.

  \code{tvcolmm_1} is a data.frame with 2700 observations on 300
  individuals at 9 occasions.

  \code{tvcolmm_2} is a data.frame with 2500 obervations on 500
  individuals at 5 occasions. 
}

\seealso{\code{\link{olmm}}, \code{\link{otsplot}}, \code{\link{tvcolmm}}}

\keyword{datasets}
