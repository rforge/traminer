\name{tvcm}

\alias{tvcm}
\alias{tvcolmm}
\alias{tvcglm}

\title{Tree-based varying coefficient regression models}

\description{
  Fits a varying coefficient regression model with a tree-structured,
  piecewise constant function for varying coefficients.
}

\usage{
tvcm(formula, data, fit, family,
     weights, subset, na.action,
     control = tvcm_control(), ...)

tvcolmm(formula, data, family = cumulative(), 
        weights, subset, na.action, 
        control = tvcm_control(), ...)

tvcglm(formula, data, family, 
       weights, subset, na.action, 
       control = tvcm_control(), ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to fit, e.g.,

    \code{y ~ x1 + vc(z1, z2, by = x2 + x3)}
    
    where \code{vc} specifies the varying coefficients.See 
    \code{\link{vc}}.}
  \item{fit}{a character string or a function that specifies the 
    fitting function, e.g. \code{\link{olmm}} or \code{\link{glm}}.} 
  \item{family}{the model family, e.g., an object of class 
    \code{family.olmm} or \code{\link{family}}.} 
  \item{data}{a data frame containing the variables in the model.}
  \item{weights}{an optional numeric vector of weights to be used in the 
    fitting process.}
  \item{subset}{an optional logical or integer vector specifying a subset
    of \code{'data'} to be used in the fitting process.}
  \item{na.action}{a function that indicates what should happen if data 
    contain \code{NA}s. See \code{\link{na.action}}.}
  \item{control}{a list with control parameters as returned by
    \code{\link{tvcm_control}}.}
  \item{\ldots}{additional arguments passed to the fitting function
    \code{fit}.}  
}

\details{
  \code{\link{tvcm}} fits a linear regression model with varying
  coefficients, i.e., functions of coefficients spanned over the value
  space of the covariates defined in the \code{\link{vc}} term. 
  \code{\link{tvcm}} fits varying coefficients by incorporating a 
  tree-structured, piecewise constant function in the fixed effect 
  predictor equation. At the time being, \code{\link{olmm}} 
  (ordinal 2-stage linear mixed models) or \code{\link{glm}} (generalized 
  linear models) can are available as fitting functions, for which are 
  provided the interfaces \code{tvcolmm} and \code{tvcglm}. 

  The algorithm works as follows: Starting with \eqn{M = 1}, the algorithm 
  splits in each step one of the current \eqn{M} nodes, i.e., subsets of the 
  values space of partitioning variables, into two new nodes. For selecting 
  the node, the variable and the cutpoint in each step, there are two options 
  (see the \code{method} argument in \code{\link{tvcm_control}}). The first 
  and default option \code{method = "mob"} selects first the node and the 
  variable based on nodewise M-fluctuation coefficient constancy tests. 
  Afterwards, the cutpoint is selected via grid search risk minimization over 
  the previously selected variable in the previously selected node. By 
  default, the \code{deviance} is used as risk function (see the \code{riskfun} 
  argument in \code{\link{tvcm_control}}). The second option 
  \code{method = "partreg"} selects simultaneously the node, the variable 
  and cutpoint by grid search risk minimization. The option 
  \code{method = "partreg"} may be slower than \code{method = "mob"} algorithm 
  and tends to select variables with many split points (variable selection bias). 
  \code{method = "mob"}, on the other hand, shouldn't be used if predictors and 
  moderators coincide in the desired model.

  The principal tuning parameter for \code{method ="mob"} is \code{alpha},
  i.e., the significance threshold at M-fluctuation tests. The algorithm
  will stop as soon as the minimum \eqn{p}-value of the tests in the current
  step are all larger than \code{alpha}. By default, \code{alpha} is chosen
  as 0.05 plus nodewise Bonferroni correction. Both \code{alpha} and
  \code{Bonferroni} can be specified via the \code{control} argument,
  see \code{\link{tvcm_control}}.  

  For \code{method = "partreg"}, the principal tuning parameter is
  \code{maxwidth}, i.e., the maximal number of terminal nodes in the
  tree structure. By default, this \code{maxwidth} is arbitrarily set to 10. 
  For a first fit with \code{method = "partreg"} should be specified 
  appropriately, e.g., with regard to the number of observations
  in \code{'data'} the number of moderators (as specified by the \code{by} 
  argument in the \code{vc} term). For selecting the final model, 
  \code{maxwidth} may be selected empirically via the in-sample AIC or
  cross-validation, see \code{\link{AIC.tvcm}} or \code{\link{cvrisk.tvcm}}.

  Special attention is given to the \code{intercept} argument in \code{\link{vc}}. 
  Conceptually, the intercepts incorporate direct impacts of moderators, 
  as specified by the \code{by} argument in \code{\link{vc}}, on the response. By 
  default, an intercept is fitted separately for each node and the overall intercept 
  is dropped. To get an overall intercept, use \code{vc(..., intercept = "ge")} or 
  \code{vc(..., intercept = "none")}. The first option, \code{intercept = "ge"}, will 
  fit a separate intercept for each node representing the deviation from the intercept
  of the node from the overall intercept. For \code{fit = "olmm"}, the option 
  \code{intercept = "ge"} implies that the intercept is global, i.e., a proportional 
  odds effect. Using \code{intercept = "none"} drops node specific intercepts, i.e., 
  direct impacts of moderators on the response. This is recommended only if the
  moderators are specified as fixed effects in the equation, as well (and therefore 
  using the \code{method = "partreg"} option).
}

\value{An object of class \code{\link{tvcm}}. The \code{\link{tvcm}} 
  class itself is based on the \code{\link{party}} class of the
  \pkg{partykit} package. The important slots are the following:
  \item{node}{an object of class \code{\link{partynode}}.}
  \item{data}{a (potentially empty) \code{\link{data.frame}}.}
  \item{fitted}{an optional \code{\link{data.frame}} with
    \code{nrow(data)} rows and containing at least the fitted terminal
    node identifiers as element \code{(fitted)}. In addition, weights
    may be contained as element \code{(weights)} and responses as
    \code{(response)}.}
  \item{info}{additional information including\code{control} and
    \code{model}.}
}

\references{ 
   Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based
   Recursive Partitioning. \emph{Journal of Computational and Graphical
     Statistics}, \bold{17}(2), 492--514.

   Torsten Hothorn and Achim Zeileis (2013). partykit: A Toolkit for
   Recursive Partytioning. URL
   \url{http://CRAN.R-project.org/package=partykit}.
}

\seealso{\code{\link{tvcm_control}}, \code{\link{tvcm-methods}},
  \code{\link{tvcm-plot}}, \code{\link{plot.tvcm}},
  \code{\link{cvrisk.tvcm}}}

\examples{

\dontrun{
## ------------------------------------------------------------------- # 
## Example 1: (Structured) heterogeneous treatment effect
##
## Fitting a varying coefficient olmm on synthetic ordinal longitudinal
## data 'tvcm_1' including 500 individuals or 2500 observations. Only
## the effect 'treat' id treated as varying effect regarding  'z1', 'z2'
## and 'z3'. The model is
##
## logit(P(Y <= j)) = a_j + u_i + b_1 time +
##                    b_1 (z_1, z_2) treatment,       u_i ~ N(0, s_u)
##
## The true treatment effects are 1.25 for z1 = 1 / z2 = 1, 0.75 for
## z1 = 1 / z2 = 2 and -0.5 for z1 = 2 and z2 = 2. z3 is nuisance. 
## ------------------------------------------------------------------- #

## load the data
data(tvcm_1)

## fit the model

m1 <- tvcolmm(y ~ fe(time) + vc(z1, z2, z3, by = treat) + re(1|id),
              data = tvcm_1, control = tvcm_control(verbose = TRUE))

## print the results
print(m1)

## ------------------------------------------------------------------- # 
## Example 2: Empirical example
##
## The algorithm is used to find out whether 'week' influences the
## coefficient of 'treat' and if 'week' has a direct impact on the
## responses.
## ------------------------------------------------------------------- #

data(schizo)

schizo.vcm <- tvcolmm(imps79o ~ vc(week, by = tx, intercept = "ce") + re(1|id),
                      data = schizo, 
                      control = tvcm_control(method = "mob", verbose = TRUE))
plot(schizo.vcm)
summary(schizo.vcm)

## in this empirical example the data are unbalanced

table(table(schizo$id))
table(schizo$week)

## for theoretical reasons, the nodewise coefficient constancy implemented
## in method = "mob" apply an auto correction in which some observations 
## are omitted to have balanced data. If we want to control this, we can 
## specify 

Nbal <- 3L # the number of obs. per individual considered as 'balanced'
subs <- !schizo$week %in% c(2, 4, 5) # the observations that are kept for tests

control <- tvcm_control(estfun = list(Nbal = Nbal, subset = subs))

## re-fit the model
schizo.vcm <- tvcolmm(imps79o ~ vc(week, by = tx, intercept = "ce") + re(1|id),
                      data = schizo, 
                      control = control)
plot(schizo.vcm)
summary(schizo.vcm)

## ------------------------------------------------------------------- #
## Example 3: Piecewise constant growth
##
## Fitting a varying coefficient olmm on synthetic ordinal longitudinal
## data 'tvcm_2' including 300 individuals or 2700 observations. The
## generating model is a cumulative random intercept model with 
## 'structual breaks' at t = 3.5 and t = 6.5.
## ------------------------------------------------------------------- #

## load the data
data(tvcm_2)

## fit the model
m3 <- tvcolmm(formula = y ~ vc(time, intercept = "ce") + re(1|id),
             data = tvcm_2, control = tvcm_control(verbose = TRUE))
## print the results
print(m3)
}
}

\author{Reto Buergin}

\keyword{models}