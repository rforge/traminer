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
    
    where \code{vc} specifies the varying coefficients. See 
    \code{\link{vcrpart-formula}}.}
  \item{fit}{a character string or a function that specifies the 
    fitting function, e.g. \code{\link{olmm}} or \code{\link{glm}}.} 
  \item{family}{the model family, e.g., an object of class 
    \code{\link{family.olmm}} or \code{\link{family}}.} 
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
  coefficients, which are functions of coefficients spanned over the value 
  space of the covariates defined in the \code{\link{vc}} term. 
  \code{\link{tvcm}} fits varying coefficients by incorporating a 
  tree-structured, piecewise constant function in the fixed effect 
  predictor equation. At the time being, \code{\link{olmm}} 
  (ordinal 2-stage linear mixed models) or \code{\link{glm}} (generalized 
  linear models) are available as fitting functions, for which are 
  provided the convenience functions \code{tvcolmm} and \code{tvcglm}. 

  The algorithm works as follows: Starting with \eqn{M = 1} node, the
  algorithm  splits in each step one of the current \eqn{M} nodes, i.e.,
  subsets of the values space of partitioning variables, into two new
  nodes. For selecting the node, the variable and the cutpoint in each
  step, there are two options (see the \code{method} argument in
  \code{\link{tvcm_control}}). The default option \code{method
  = "mob"} selects first the node and the variable based on nodewise
  M-fluctuation coefficient constancy tests (cf. Zeileis and Hornik,
  2007) and the cutpoint is selected via grid search risk minimization 
  over the selected variable and node. By default, the \code{deviance} 
  is used as risk function (see the \code{riskfun} argument in
  \code{\link{tvcm_control}}). The second option 
  \code{method = "partreg"} selects simultaneously the node, the 
  variable and cutpoint by grid search risk minimization 
  (cf. Wang and Hastie, 2013). The option \code{method = "partreg"} may 
  be slower than \code{method = "mob"} algorithm and tends to select 
  variables with many split points (variable selection bias). 
  \code{method = "mob"}, on the other hand, shouldn't be used if 
  predictors and moderators coincide in the desired model.  

  The principal tuning parameter for \code{method ="mob"} is
  \code{alpha}, i.e., the significance threshold at M-fluctuation
  tests. The algorithm will stop as soon as the minimum \eqn{p}-value of
  the tests in the current step are all larger than \code{alpha}. By
  default, \code{alpha} is chosen as 0.05 plus nodewise Bonferroni
  correction. Both \code{alpha} and \code{Bonferroni} can be specified
  via the \code{control} argument, see \code{\link{tvcm_control}}.
  Notice that coefficient constancy tests are adjusted for intra-subject 
  correlation for 2-stage models, see \code{\link{estfun.olmm}}. 

  For \code{method = "partreg"}, the principal tuning parameter is
  \code{maxwidth}, i.e., the maximal number of terminal nodes in the
  tree structure. By default, this \code{maxwidth} is arbitrarily set to
  10. \code{maxwidth} should therefore be set to an appropriate value 
  using \code{method = "partreg"}, e.g., with regard to the number 
  of observations in \code{'data'} the number of partitioning variables
  (as specified by the \code{by} argument in the \code{vc} term). For 
  selecting the final model, \code{maxwidth} may be selected empirically
  via the in-sample AIC or cross-validation, see \code{\link{tvcm-cv}}. 

  Special attention is given to the \code{intercept} argument in
  \code{\link{vc}}. Conceptually, the intercepts incorporate direct
  impacts of partitioning variables i.e. moderators on the response.  
  By default, an intercept is fitted separately for each node and the 
  overall intercept is dropped. To get an overall intercept, use 
  \code{vc(..., intercept = "ge")} or \code{vc(..., intercept = "none")}. 
  The first option, \code{intercept = "ge"}, will fit a separate  
  intercept for each node representing the deviation from the intercept 
  of the node from the overall intercept. For \code{fit = "olmm"}, the 
  option \code{intercept = "ge"} implies that the intercept is global,
  i.e., a proportional odds effect for the \code{\link{cumulative}} logit model. 
  Using \code{intercept = "none"} drops node specific intercepts, i.e., 
  direct impacts of moderators on the response. This is recommended only 
  if the moderators are specified as fixed effects in the equation. 
  
  The parameter \code{minbucket}, i.e., the threshold for the minimum
  number of observations in each node, is to be set, for both 
  \code{method = "mob"} and \code{method = "partreg"}, with respect to
  the complexity of the regression model and should be large enough to
  ensure that the regression model can be fitted, cf. Wang and Hastie
  (2013). We recommend, in particular for ordinal regression models,
  values larger than 25.
}

\value{An object of class \code{\link{tvcm}}. The \code{\link{tvcm}} 
  class itself is based on the \code{\link{party}} class of the
  \pkg{partykit} package. The most important slots are:
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
   Zeileis, A., Hothorn, T., and Hornik, K. (2008). Model-Based Recursive 
   Partitioning. \emph{Journal of Computational and Graphical Statistics}, 
   \bold{17}(2), 492--514.

   Zeileis, A., Hornik, K. (2007), Generalized M-Fluctuation Tests for 
   Parameter Instability, \emph{Statistica Neerlandica}, \bold{61}, 488--508.
   doi:10.1111/j.1467-9574.2007.00371.x.

   Hothorn, T. and Zeileis A. (2013). partykit: A Toolkit for
   Recursive Partytioning. \url{http://CRAN.R-project.org/package=partykit}.

   Wang, J. C., Hastie, T. (2013), Boosted Varying-Coefficient Regression
   Models for Product Demand Prediction, \emph{Journal of Computational
   and Graphical Statistics}.
}

\seealso{\code{\link{tvcm_control}}, \code{\link{tvcm-methods}},
  \code{\link{tvcm-plot}}, \code{\link{tvcm-plot}},
  \code{\link{tvcm-cv}}}

\examples{

## ------------------------------------------------------------------- # 
## Example 1: Moderated treatment effect
##
## Fitting a varying coefficient olmm on the synthetic ordinal 
## longitudinal data 'vcrpart_1'. 'z1', 'z2', 'z3' and 'z4' are the 
## moderators i.e. partitioning variables and 'treat' is the 
## variable whose effect is moderated. The coefficient of 'wave' is 
## assumed to be homogeneous. 'z4' is the only true moderator
## ------------------------------------------------------------------- #

data(vcrpart_1)

## fit the model
model.1 <- tvcolmm(y ~ fe(wave) + vc(z1, z2, z3, z4, by = treat) + re(1|id),
                  data = vcrpart_1, control = tvcm_control(verbose = TRUE))

plot(model.1)

## ------------------------------------------------------------------- # 
## Example 2: Schizophrenia
##
## The algorithm is used to find out whether 'week' influences the
## coefficient of 'tx' (treatment) and if 'week' has a direct impact on the
## responses.
## ------------------------------------------------------------------- #

data(schizo)

subs <- schizo$week == 0 | schizo$week == 1 | schizo$week == 3
schizo <- schizo[subs,]
schizo$id <- droplevels(schizo$id)

subs <- (table(schizo$id) == 3)[schizo$id]
schizo <- schizo[subs,]
schizo$id <- droplevels(schizo$id)

model.2 <- tvcolmm(imps79o ~ vc(week, by = tx, intercept = "ce") + re(1|id),
                   data = schizo, control = tvcm_control(verbose = TRUE))

plot(model.2, "coef", terms = list(paste("Eta", 1:3, ":(Intercept)", sep = ""), "tx"), 
     tnex = 2)
}

\author{Reto Buergin}

\keyword{models}