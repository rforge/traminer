\name{fvcm}

\alias{fvcm}
\alias{oobrisk.fvcm}
\alias{plot.fvcm}
\alias{predict.fvcm}
\alias{print.fvcm}

\title{Bagging and random forests based on \code{\link{tvcm}}.}

\description{
  Implements bagging and random forests for \code{\link{tvcm}}. 
}

\usage{

fvcm(..., folds = cvfolds(object, "subsampling", 50L), control = fvcm_control())

\method{oobrisk}{fvcm}(object, fun = NULL, ...)
\method{predict}{fvcm}(object, newdata = NULL,
        type = c("link", "response", "prob", "class", "coef"),
        ranef = FALSE,  verbose = FALSE, ...)
\method{print}{fvcm}(x, ...)
}

\arguments{
  \item{...}{for \code{\link{fvcm}}, this part is used to assign
    arguments for \code{\link{tvcm}}, including \code{formula}, 
    \code{data}, \code{family} and \code{fit}. For the other 
    functions, ... are arguments to be passed to other methods.}
  \item{object,x}{an object of class \code{\link{fvcm}}.}
  \item{fun}{the risk function. The default risk function is defined
    as the sum of the deviance residuals.}
  \item{folds}{a matrix of integers defining the folds for
    cross-validation. See \code{\link{cvrisk.tvcm}}.}
  \item{control}{a \code{list} that specifies the maximum control
    parameters to evaluate. See \code{\link{fvcm_control}}.}
  \item{newdata}{an optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{type}{character string indicating the type of prediction. See
    \code{\link{predict.tvcm}}.}  
  \item{ranef}{logical scalar. This argument is used only if the
    \code{\link{olmm}} function is used. See
    \code{\link{predict.olmm}}.}  
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the validation.}
}

\details{
  Implements "bagging" (Breiman, 1996) and "random forests" (Breiman,
  2001) for \code{\link{tvcm}} objects. The method consist in
  growing multiple trees using \code{\link{tvcm}} and aggregating the
  the predictor equations. For bagging, simply use 
  \code{control = list(mtry = Inf)}. The implementation may be slow.

  \code{oobrisk.fvcolmm} estimates the out-of-bag risk by aggregating 
  the linear predictor for each observations with those trees in which
  the observation didn't appear and using that aggregated predictor 
  to predict the outcome.

  \code{predict.fvcolmm} aggregates the predictors of all fitted trees.  
}

\references{
  Leo Breiman (1996). Bagging Predictors. \emph{Machine Learning},
  123--140 
  Leo Breiman (2001). Random Forests. \emph{Machine Learning},
  \bold{45}(1), 5--32. 
}

\examples{

\dontrun{
## Extend Example 2 in '?tvcm'

## load the data
data(tvcm_2)

## fit the model

folds <- cvfolds(tvcm_2$id, type = "subsampling", K = 2)

f1 <- fvcm(y ~ time + vc(z1,z2,z3,z4,by = treat) + re(1|id),
              data = tvcm_1,
              control = tvcm_control(mtry = 2,
                minbucket = 50, maxevalsplit = 5,
                verbose = TRUE),
              folds = folds)
}
}

\value{An object of class \code{fvcm}.}

\author{Reto Buergin}

\keyword{models}