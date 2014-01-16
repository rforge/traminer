\name{tvcolmm}

\alias{tvcolmm}
\alias{tvcolmm-class}

\title{Ordinal linear mixed models with varying (fixed-) coefficients}

\description{
  \code{\link{tvcolmm}} estimates an ordinal linear mixed model with
  varying fixed effects via model-based recursive partitioning. The
  result is a \code{\link{olmm}} object with piecewise constant varying
  effects, where the regions of the piecewise function result from
  recursive partitioning in a preliminary step.
}

\usage{
tvcolmm(formula, data, control = tvcolmm_control(),
        subset, na.action, offset,
        vi = c("none", "po", "npo"),
        linear = NULL, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to fit. This should
    be of form

    \code{y ~ x1 + \ldots + (1 + w1 + \ldots | id) | xk + \ldots | z1 + \ldots} 

    where the first two parts represents the varying coefficient model
    without effect modifiers model and \code{z1 + \ldots} are the
    effect modifiers for the fixed effects in the first part of the
    equation (the middle part may not be necessary,
    cf. \code{\link{olmm}}). If single coefficients should not be
    treated as varying, or varying intercepts are to be estimated, see
    details.}   
  \item{data}{a data frame containing the variables in the model.}
  \item{control}{a list with control parameters as returned by
    \code{\link{tvcolmm_control}}.}
  \item{subset, na.action, offset}{further model
    specification arguments as in \code{lm}.}
  \item{vi}{character string. Defines whether a main effect for the
    partiton structure should be added. \code{po} adds a
    predictor-invariant (proportional odds) and \code{npo} a
    predictor-variable intercept (non-proportional odds). Note that this
    means not that the algorithm is powerful to detect variables that
    alter those intercepts, see details.}
  \item{linear}{vector of characters. Indicates the terms that are
    assumed to be non-varying i.e. truely linear i.}
  \item{\ldots}{additional arguments passed to the \code{\link{olmm}}
    call for the model} 
}

\details{
  The \code{\link{tvcolmm}} function is an procedure for varying fixed
  effects in ordinal linear mixed models. Fixed effects are approximated
  by a piecewise constant function over partitions that are gained from
  a Breadth-first partitioning algorithm. The underlying model can be
  extracted by \code{extract(object, "model")}.
  
  The algorithm is based on the model-based partitioning (MOB) algorithm
  of Zeileis et al. (2008). In difference, it does however incorporate a
  tree structure in the fixed effect specification instead of entirely
  divide the data. As a result, the programm yields one model for the
  entire data set. 

  The algorithm initializes the value space of variables \code{z1,
  \ldots, zl} as start partition and proceeds as follows:

  \enumerate{
    \item Fit an ordinal linear mixed model to all observations with
    separate fixed effects for each partition.
    \item Assess the constancy of fixed effects with respect to
    each partitioning variable \code{z1, \ldots, zl}. If coefficients vary
    with respect to at least one variable, choose the variable \code{zj} 
    the induces the highest variation, otherwise stop.
    \item Search for among all possible splits in each partition the
    optimal split in \code{zj} by maximum likelihood optimization.
    \item Add the best split to the model specification and repeat the
    procedure.
  }

  The tree size is automatically selected by the stopping criteria in
  step 2. This requires a parameter \code{alpha} that can be specified
  using the \code{control} argument. The default is 0.05 and may be
  specified via cross-validation (see \code{\link{cv.tvcolmm}}). 

  By default, the algorithm is not sensitive to direct influences of z
  variables on the response. If this is desired, set \code{vi = "npo"}
  and something like \code{control = vcolmm_control(terms =
  c("(Intercept)", ...))}. See examples.
}

\value{An object of class \code{\link{tvcolmm}}. The \code{\link{tvcolmm}} 
  class itself is based on the \code{\link{party}} class of the
  \pkg{partykit} package and has the following slots:
  \item{node}{an object of class \code{\link{partynode}}.}
  \item{data}{a (potentially empty) \code{\link{data.frame}}.}
  \item{fitted}{an optional \code{\link{data.frame}} with
    \code{nrow(data)} rows and containing at least the fitted terminal
    node identifiers as element \code{(fitted)}. In addition, weights
    may be contained as element \code{(weights)} and responses as
    \code{(response)}.}
  \item{terms}{an optional \code{\link{terms}} object.}
  \item{info}{additional information comprising \code{title}
    (character), \code{call} (an object of class \code{call}), 
    \code{control} (an object of class \code{\link{tvcolmm_control}}) and
    \code{model} (an object of class \code{\link{olmm}}, the
    fitted varying coefficient model).}
}

\references{ 
   Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based
   Recursive Partitioning. \emph{Journal of Computational and Graphical
     Statistics}, \bold{17}(2), 492--514.

   Torsten Hothorn and Achim Zeileis (2013). partykit: A Toolkit for
   Recursive Partytioning. URL
   \url{http://CRAN.R-project.org/package=partykit}.
}

\seealso{\code{\link{tvcolmm_control}}, \code{\link{tvcolmm-methods}},
  \code{\link{tvcolmm-plot}}, \code{\link{plot.tvcolmm}},
  \code{\link{cv.tvcolmm}}}

\examples{

\dontrun{
## ------------------------------------------------------------------- # 
## Example 1: (Structured) heterogeneous treatment effect
##
## Fitting a varying coefficient olmm on synthetic ordinal longitudinal
## data 'tvcolmm_2' including 500 individuals or 2500 observations. Only
## the effect 'treat' id treated as varying effect regarding  'z1', 'z2'
## and 'z3'. The model is
##
## logit(P(Y <= j)) = a_j + u_i + b_1 time +
##                    b_1 (z_1, z_2) treatment,       u_i ~ N(0, s_u)
##
## The true treatment effects are 1.25 for z1 = 1 / z2 = 1, 0.75 for
## z1 = 1 / z2 = 2 and -0.5 for z1 = 2 and z2 = 2. z3 is a true
## nuisance variable. 
## ------------------------------------------------------------------- #

## load the data
data(tvcolmm_2)

## fit the model

mod1 <- tvcolmm(y ~ time + treat + (1|id) | z1 + z2 + z3, data = tvcolmm_2,
                control = tvcolmm_control(verbose = TRUE,
                  minsplit = 30L, terms = c("treat")), 
                linear = "time", vi = "po")

## print the results
print(mod2)

## ------------------------------------------------------------------- # 
## Example 2: Empirical example
## ------------------------------------------------------------------- #

data(schizo)
mod3 <- tvcolmm(imps79o ~ tx + (1|id) | week,
                data = schizo,
                control = tvcolmm_control(verbose = TRUE),
                vi = "po")

## ------------------------------------------------------------------- #
## Example 3: Piecewise constant growth
##
## Fitting a varying coefficient olmm on synthetic ordinal longitudinal
## data 'tvcolmm_1' including 300 individuals or 2700 observations. The
## generating model is a cumulative random intercept model with 
## 'structual breaks' at t = 3.5 and t = 6.5.
## ------------------------------------------------------------------- #

## load the data
data(tvcolmm_1)

## fit the model
mod3 <- tvcolmm(formula = y ~ (1|id) | time,
                data = tvcolmm_1, vi = "npo",
                control = tvcolmm_control(verbose = TRUE,
                  terms = "(Intercept)")) # 'terms = "(Intercept)"' is necessary!

## print the results
print(mod3)
}
}

\author{Reto Buergin}

\keyword{models}