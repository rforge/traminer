\name{olmm}

\alias{family.olmm}
\alias{cumulative}
\alias{adjacent}
\alias{baseline}
\alias{olmm}

\title{Fitting ordinal two-stage linear mixed models}

\description{
  Fits different types of two-stage linear mixed models for longitudinal
  (or clustered) ordinal (or multinomial) responses. one-stage models
  are also allowed. Random effects are assumed to be multivariate normal 
  distributed with expectation 0. At the time being, cumulative link 
  models with the logit, probit or cauchy link, the baseline-category 
  logit and the adjacent-category logit model are
  implemented. Coefficients can be category-specific
  (i.e. non-proportional odds effects) or global (i.e. proportional
  odds, or parallel effects). 

  The function solves the score function for coefficients of the
  marginal likelihood by using Gauss-Hermite quadrature (e.g., Hedeker;
  1994). Random effects are predicted by their expectation (see Hartzl
  et al.; 2001). Standard deviations of parameter estimates are, by
  default, based on the expected Fisher-information matrix. 
}

\usage{
cumulative(link = c("logit", "probit", "cauchy"))
adjacent(link = "logit")
baseline(link = "logit")

olmm(formula, data, family = cumulative(), 
     weights, subset, na.action,
     offset, contrasts, control = olmm_control(), ...)
}

\arguments{
  \item{formula}{a symbolic description of the model. This should be
    something like
    
    \code{y ~ ce(x1) + ge(x2) +re(1 + ge(w2) | id)} 

    where \code{ce(x1)} specificies that the predictor \code{x1} has a 
    category-specific i.e. non-proportional odds effect and \code{ge(x2)} 
    has global i.e. proportional odds fixed effect, see
    \code{\link{ge}}, resp. \code{\link{ce}}. Random effects are
    specified within the \code{\link{re}} term, where the variable
    \code{id} above behind the vertical bar \code{|} defines  the
    subject i.e. cluster factor. Notice that only one subject factor is
    allowed.}   
  \item{data}{an optional data frame with the variables
    in \code{formula}. By default the variables are taken from the 
    environment from which \code{olmm} is called.}
  \item{family}{an \code{family.olmm} object produced by
    \code{cumulative}, \code{adjacent} or \code{baseline}.}
 \item{weights}{a numeric vector of weights with length equal the
    number of observations. The weights should be constant for
    subjects.} 
  \item{offset}{a matrix specifying the offset separately for each
    predictor equation, of which there are the number of categories of
    the response minus one.}
  \item{subset, na.action, contrasts}{further model
    specification arguments as in \code{\link{lm}}.} 
  \item{control}{a list of control parameters produced by \code{\link{olmm_control}}.}
  \item{link}{character string. The name of the link function.}
  \item{...}{arguments to be passed to \code{control}.}
}

\details{The function can be used to fit simple ordinal two-stage
  mixed effect models with up to 3-4 random effects. For
  models with higher dimensions on random effects, the procedure may
  not convergence (cf. Tutz; 1996). Coefficients for the adjacent-category 
  logit model are extracted via coefficient transformation 
  (e.g. Agresti; 2010).

  The three implemented families are defined as follows:
  \code{\link{cumulative}} is defined as the link of the sum of
  probabilities of lower categories, e.g., for \code{link = "logit"},
  the logit of the sum of probabilities of lower
  categories. \code{\link{adjacent}} is defined as the logit of 
  the probability of the lower of two adjacent
  categories. \code{\link{baseline}} is defined as the logit of the
  probability of a category with reference to the highest
  category. Notice that the estimated coefficients of cumulative models
  may have the opposite sign those obtained with alternative software.  
  
  For alternative fitting functions, see fore example \code{clmm} of
  package \code{ordinal}, \code{cumlogitRE} of package \code{glmmAK},
  \code{nplmt} of package \code{mixcat}, \code{DPolmm} of package
  \code{DPpackage}, \code{lcmm} of package \code{lcmm}, \code{MCMCglmm}
  of package \code{MCMCglmm}, \code{sabre} of package \code{sabreR} or
  \code{OrdinalBoost} of package \code{GMMBoost}. 
}

\value{\code{\link{olmm}} returns an object of class
  \code{\link{olmm}}. \code{\link{cumulative}},  \code{\link{adjacent}}
  and \code{\link{baseline}} yield an object of class
  \code{family.olmm}. The \code{\link{olmm}} class is a list containing
  the following components:

  \item{env}{environment in which the object was built.}  
  \item{frame}{the model frame.}
  \item{call}{the matched call to the function that created the object
    (class \code{"call"}).} 
  \item{control}{a list of class \code{olmm_control} produced by
    \code{\link{olmm_control}}.}  
  \item{formula}{the formula of the call.}
  \item{terms}{a list of \code{\link{terms}} of the fitted model.}
  \item{family}{an object of class \code{family.olmm} that specifies
    that family of the fitted model.} 
  \item{y}{(ordered) categorical response vector.} 
  \item{X}{model matrix for the fixed effects.}
  \item{W}{model matrix for the random effects.}
  \item{subject}{a factor vector with grouping levels.}
  \item{subjectName}{variable name of the subject vector.}
  \item{weights}{numeric observations weights vector.}
  \item{weights_sbj}{numeric weights vector of length N.}
  \item{offset}{numeric offset matrix}
  \item{xlevels}{(only where relevant) a list of levels of the
    factors used in fitting.} 
  \item{contrasts}{(only where relevant) a list of contrasts
    used.} 
  \item{dims}{a named integer of dimensions. Some of the
    dimensions are \eqn{n} is the number of observations, \eqn{p} is
    the number of fixed effects per predictor and \eqn{q} is the total
    number of random effects.}
  \item{fixef}{a matrix of fixed effects (one column for each
    predictor).}
  \item{ranefCholFac}{a lower triangular matrix. The cholesky
    decomposition of the covariance matrix of the random effects.}
  \item{coefficients}{a numeric vector of several fitted model
    parameters}
  \item{restricted}{a logical vector indicating which elements
    of the \code{coefficients} slot are restricted to an initial value
    at the estimation.}
  \item{eta}{a matrix of unconditional linear predictors of
    the fixed effects without random effects.}  
  \item{u}{a matrix of orthogonal standardized random
    effects (one row for each subject level).} 
  \item{logLik_obs}{a numeric vector of log likelihood value
    (one value for each observation).}
  \item{logLik_sbj}{a numeric vector of log likelihood values
    (one value for each subject level).}
  \item{logLik}{a numeric value. The log likelihood of the
    model.} 
  \item{score_obs}{a matrix of observation-wise partial
    derivates of the marginal log-likelihood equation.} 
  \item{score_sbj}{a matrix of subject-wise partial derivates
    of the marginal log-likelihood equation.} 
  \item{score}{a numeric vector of (total) partial derivates
    of the log-Likelihood function.}
  \item{info}{the information matrix (default is the expected
    information).}
  \item{ghx}{a matrix of quadrature points for the
    Gauss-Hermite quadrature integration.} 
  \item{ghw}{a matrix of weights for the Gauss-Hermite
    quadrature integration.} 
  \item{ranefElMat}{a transformation matrix}
  \item{optim}{a list of arguments for calling the optimizer
    function.}
  \item{control}{a list of used control arguments produced by
    \code{\link{olmm_control}}.} 
  \item{output}{the output of the optimizer (class
    \code{"list"}).}
}

\author{Reto Buergin}

\references{
  Agresti, A. (2010). \emph{Analysis of Ordinal Categorical Data}, 10
  edn, Wiley.
  
  Hartzel, J., Agresti A. and Caffo, B. (2001). Multinomial Logit Random
  Effect Models, \emph{Statistical Modelling} \bold{1}: 81--102
  
  Hedeker, D. and Gibbons, R. (1994). A random-effects ordinal
  regression model for multilevel analysis, \emph{Biometrics} \bold{20}
  (4): 933--944

  Tutz, G. and Hennevogl W. (1996). Random effects in ordinal regression
  models, \emph{Computational Statistics & Data Analysis} \bold{22} (5):
  537--557

  Tutz, G. (2012). \emph{Regression for Categorical Data}, Cambridge
  Series in Statistical and Probabilistic Mathematics.
}

\seealso{\code{\link{olmm-methods}}, \code{\link{olmm_control}},
  \code{\link{ordered}}}  

\examples{
## ------------------------------------------------------------------- #
## Example 1: Schizophrenia
##
## Estimating the cumulative mixed models of
## Agresti (2010) chapters 10.3.1
## ------------------------------------------------------------------- #

data(schizo)

model.10.3.1 <-
  olmm(imps79o ~ tx + sqrt(week) + tx * sqrt(week) + re(1|id),
       data = schizo, family = cumulative())

summary(model.10.3.1)

## ------------------------------------------------------------------- #
## Example 2: Movie critics
##
## Estimating three of several adjacent-categories
## mixed models of Hartzl et. al. (2001)
## ------------------------------------------------------------------- #

data(movie)

## model with category-specific effects for "review"
model.24.1 <- olmm(critic ~  ce(review) + re(1|movie, intercept = "ce"),
                   data = movie, family = adjacent())

summary(model.24.1)
}

\keyword{models}