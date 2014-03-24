\name{olmm}

\alias{olmm}

\title{Fitting ordinal 2-stage linear mixed models}

\description{
  Fits different types of 2-stage linear mixed models for clustered or
  longitudinal ordinal responses. Random effects are assumed to be
  multivariate normal distributed with expectation 0. At the time being,
  cumulative link models with the logit, probit or cauchy link, the
  baseline-category (baseline category = superior category) and the
  adjacent-category logit model are available. Coefficients can be
  variable (i.e. non-proportional odds effects) or invariant
  (i.e. proportional odds effects) across predictors.

  The fitting function uses direct marginal likelihood maximization
  using Gauss-Hermite quadrature (e.g. Hedeker; 1994). Random effects
  are predicted by their expectation (see Hartzl et al.; 2001). Standard
  deviations of parameter estimates are, by default, based on the expected
  Fisher-information matrix.
}

\usage{
olmm(formula, data, weights, start, subset, na.action,
     offset, contrasts,
     family = c("cumulative", "baseline", "adjacent"),
     link = c("logit", "probit", "cauchy"),
     doFit = TRUE, optim = list(),
     numGrad = FALSE, numHess = numGrad, nGHQ = 7L,
     restricted, verbose = FALSE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model. This should be of
    form 
    
    \code{y ~ x1 + \ldots + (1 + w1 + \ldots | id) | xk + \ldots + (\ldots | id)} 
    
    where \code{x1 + \ldots} describes the predictor-invariant
    i.e. proportional odds fixed coefficient terms, 
    \code{(1 + w1 + \ldots | id)} the predictor-invariant random
    coefficient terms and \code{xk + \ldots + (\ldots | id)} the
    predictor-variable i.e. non-proportional odds fixed and random
    coefficient terms. The response (e.g. \code{y}) must be an
    \code{\link{ordered}} factor and the subject variable
    (e.g. \code{id}) a \code{\link{factor}}. Binary responses are
    allowed.} 
  \item{data}{an optional data frame containing the variables named in
    \code{formula}. By default the variables are taken from the
    environment from which \code{olmm} is called.}
  \item{start}{a named numeric vector of initial values for the
    parameters. The parameter must be named in exactly in the way as
    they appear when the model is fitted.}
  \item{weights}{a numeric vector of weights with length equal the
    number of subjects. Longitudinal weights are not supported.}  
  \item{subset, na.action, offset, contrasts}{further model
    specification arguments as in \code{lm}.}
  \item{family}{character string for ordinal link. Available are
    \code{"cumulative"} (default), \code{"baseline"} and
    \code{"adjacent"}.} 
  \item{link}{the link function for cumulative link models
    (\code{family = "cumulative"}). Available are \code{"logit"} (default),
    \code{"probit"} and \code{"cauchy"}.}
  \item{doFit}{logical scalar. When \code{FALSE} an unfitted \code{olmm} 
    object is returned. See details.}
  \item{optim}{a list for specifying the algorithm that minimizes
    the Likelihood function. By default, \code{olmm} uses the
    \code{ucminf} function of package \code{ucminf}. Additionally, the
    \code{optim} and the \code{nlminb} functions can be used. See
    details.} 
  \item{numGrad}{logical scalar indicating whether the score function
    should be retrieved numerically.}
  \item{numHess}{logical scalar. Indicates whether the Hess matrix for
    the variance-covariance matrix should be estimated numerically (this
    is then an approximation of the observed fisher 
    information). Must be \code{TRUE} if \code{numGrad} is
    \code{TRUE}. See details.} 
  \item{nGHQ}{a positive integer specifying the number of quadrature
    points for the approximation of the marginal Likelihood by numerical
    integration.}
  \item{restricted}{a character vector of names of coefficients to be
    fixed during the optimization. The argument is ignored in case of
    adjacent category models.}
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the optimization of the parameter estimates.}
  \item{...}{further undocumented arguments to be passed.}
}

\details{The function can be used to fit simple ordinal 2-stage
  mixed effect models with up to 3-4 random effects. For
  models with higher dimensions on random effects, the procedure may
  shows difficulties with convergence (cf. Tutz; 1996). 

  Cumulative link models are fit with the direct marginal maximum
  likelihood approach, as suggested by Hedeker (1994). Adjacent- and
  baseline-category models are fitted as well using direct marginal
  likelihood maximization, using for both the likelihood of the
  baseline-category model. Coefficients for the adjacent-category model
  are retrieved via coefficient transformation (e.g. Agresti; 2010).

  Initial values may improve the computation time and avoid
  divergence. The \code{start} argument accepts a vector with named
  elements according to the column labels from the
  \code{\link{model.matrix}}. To retrieve these labels, one may pre-run
  the \code{olmm} function with the argument \code{doFit = FALSE} and
  then look at the names in the \code{coefficients} slot. At the time
  being, initial values for adjacent-categories models must be
  transformed into the baseline-category model form.

  The argument \code{optim} may be used to modify the optimizer. For
  example, the Nelder-Mead algorithm of \code{\link{optim}} function can
  be called by \code{olmm(..., optim = list(method = "Nelder-Mead"))} or
  the \code{\link{nlminb}} function by \code{olmm(..., optim =
  list(method = "nlminb")).}

  \code{weights} should be used with caution. \code{olmm} can handle
  time-invariant weights only (one distinct weight per subject). Weights
  are used for computing the log-likelihood and the expected fisher
  information matrix (as suggested in Hedeker(1994), sec. 3.2). For
  numerical approximations of the (observed) fisher information
  (i.e. \code{numGrad = TRUE}, \code{numHess = TRUE}) the weights are
  currently ignored.

  For alternative fitting functions, see fore example \code{clmm} of
  package \code{ordinal}, \code{cumlogitRE} of package \code{glmmAK},
  \code{nplmt} of package \code{mixcat}, \code{DPolmm} of package
  \code{DPpackage}, \code{lcmm} of package \code{lcmm}, \code{MCMCglmm}
  of package \code{MCMCglmm}, \code{sabre} of package \code{sabreR} or
  \code{OrdinalBoost} of package \code{GMMBoost}. 
}

\value{An object of class \code{\link{olmm}}. See
  \code{\link{olmm-class}} for details.}

\author{Reto Buergin}

\references{
  Agresti, A. (2010). \emph{Analysis of Ordinal Categorical Data}, 10
  edn, Wiley.
  
  Hartzel, J., Agresti A. and Caffo, B. (2001). Multinomial Logit Random
  Effect Models, \emph{Statistical Modelling} \bold{1}: 81--102
  
  Hedeker, D. and Gibbons, R. (1994). A random-effects ordinal
  regression model for multilevel analysis. \emph{Biometrics} \bold{20}
  (4): 933--944

  Tutz, G. and Hennevogl W. (1996). Random effects in ordinal regression
  models, \emph{Computational Statistics & Data Analysis} \bold{22} (5):
  537--557
}

\seealso{\code{\link{ordered}}, \code{\link{olmm-class}}}  

\examples{
## ------------------------------------------------------------------- #
## Example 1: Schizophrenia
##
## Estimating the cumulative mixed models of
## Agresti (2010) chapters 10.3.1 (random intercept
## model) and 10.3.3 (random slope model)
## ------------------------------------------------------------------- #

data(schizo)

model.10.3.1 <-
  olmm(imps79o ~ tx + sqrt(week) + tx * sqrt(week) + (1 | id),
       data = schizo, family = "cumulative",
       link = "logit")

summary(model.10.3.1)


model.10.3.3 <-
  olmm(imps79o ~ tx + sqrt(week) + tx * sqrt(week) + (1 + sqrt(week) | id),
       data = schizo, family = "cumulative",
       link = "logit")

summary(model.10.3.3)


## ------------------------------------------------------------------- #
## Example 2: Movie critics
##
## Estimating three of several adjacent-categories
## mixed models of Hartzl et. al. (2001)
## ------------------------------------------------------------------- #

data(movie)

## model with predictor-variable random intercepts
model.22.1 <- olmm(critic ~ review | (1 | movie), data = movie,
                   family = "adjacent")

summary(model.22.1)


## model with predictor-invariant random intercepts
model.22.2 <- olmm(critic ~ review + (1 | movie), data = movie,
                   family = "adjacent")

summary(model.22.2)


## as model 22.1 but with predictor-variable effects for "review"
model.24.1 <- olmm(critic ~ 1 | review + (1 | movie), data = movie,
                   family = "adjacent")

summary(model.24.1)
}

\keyword{models}