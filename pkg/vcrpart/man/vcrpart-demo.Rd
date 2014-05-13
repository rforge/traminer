\name{vcrpart-demo}

\alias{vcrpart_1}
\alias{vcrpart_2}

\docType{data}

\title{Synthetic data sets}

\description{Synthetic data set for illustrations.} 

\usage{
data(vcrpart_1)
data(vcrpart_2)
}

\format{
  \describe{
    \item{\code{y}}{ordered factor. The response variable}
    \item{\code{id}}{factor. The subject identification vector.}
    \item{\code{wave}}{numeric. The wave identification vector.}
    \item{\code{treat}}{a dummy variable. The treatment effect.}
    \item{\code{x1, x2}}{numeric predictor variables.}
    \item{\code{z1, z2, z3, z2}}{partitioning/ varying-coefficients variables.}
  }
}

\seealso{\code{\link{olmm}}, \code{\link{otsplot}}, \code{\link{tvcm}}}

\examples{
## --------------------------------------------------------- #
## generating vcrpart_1
## --------------------------------------------------------- #

## create skeletton
set.seed(1)
vcrpart_1 <- data.frame(id = factor(rep(1:50, each = 4)),
                        wave = rep(1:4, 50),
                        treat = sample(0:1, 200, TRUE))

## add partitioning variables
vcrpart_1$z1 <- rnorm(50)[vcrpart_1$id]
vcrpart_1$z2 <- rnorm(200)
vcrpart_1$z3 <- factor(sample(1:2, 50, TRUE)[vcrpart_1$id])
vcrpart_1$z4 <- factor(sample(1:2, 200, TRUE))

## simulate response
eta <- 2 * vcrpart_1$treat * (vcrpart_1$z4 == "1")
eta <- eta + rnorm(50)[vcrpart_1$id] + rlogis(200)
vcrpart_1$y <- cut(-eta, c(-Inf, -1, 1, Inf), 1:3,
                   ordered_result = TRUE)


## --------------------------------------------------------- #
## generating vcrpart_2
## --------------------------------------------------------- #

set.seed(1)
vcrpart_2 <- data.frame(x1 = rnorm(100),
                  x2 = rnorm(100),
                  z1 = factor(sample(1:3, 100, TRUE)),
                  z2 = factor(sample(1:3, 100, TRUE)))
vcrpart_2$y <- vcrpart_2$x1 * (vcrpart_2$z1 == "2") +
  2 * vcrpart_2$x1 * (vcrpart_2$z1 == "3")
vcrpart_2$y <- vcrpart_2$y + rnorm(100)
}

\keyword{datasets}