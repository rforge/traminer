\name{olmm}

\alias{cumulative}
\alias{adjacent}
\alias{baseline}
\alias{olmm}

\title{Fitting ordinal 2-stage linear mixed models}

\description{
  Fits different types of 2-stage linear mixed models for clustered or
  longitudinal ordinal (or multinomial) responses. 1-stage data are 
  allowed as well. Random effects are assumed to be multivariate normal 
  distributed with expectation 0. At the time being, cumulative link 
  models with the logit, probit or cauchy link, the baseline-category 
  logit (baseline category = superior category) and the adjacent-category
  logit model are implemented. Coefficients can be category-specific
  (i.e. non-proportional odds effects) or global (i.e. proportional odds, 
  or parallel effects).

  The function solves the score function i.e. the partial derivates for
  coefficients of marginal likelihood using Gauss-Hermite quadrature 
  (e.g., Hedeker; 1994). Random effects are predicted by their expectation 
  (see Hartzl et al.; 2001). Standard deviations of parameter estimates 
  are, by default, based on the expected Fisher-information matrix.
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
  \item{formula}{a symbolic description of the model. This should be of
    form 
    
    \code{y ~ ge(x1) + ce(x2) +re(1 + ce(w2) | id)} 

    where \code{ge(x1)} specificies that the predictor \code{x1} has a 
    global i.e. proportional odds fixed effect and \code{ce(x2)} has a
    category-specific i.e. non-proportional odds effect on the logits of 
    the response, see \code{\link{ge}}, resp. \code{\link{ce}}. Random 
    effects are specified within the \code{\link{re}} term, where the 
    variable \code{id} above behind the vertical bar \code{|} defines 
    the subject i.e.cluster factor.} 
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

\details{The function can be used to fit simple ordinal 2-stage
  mixed effect models with up to 3-4 random effects. For
  models with higher dimensions on random effects, the procedure may
  not convergence (cf. Tutz; 1996). Coefficients for the adjacent-category 
  logit model are extracted via coefficient transformation 
  (e.g. Agresti; 2010).

  Initial values may decrease the computation time and avoid
  divergence. The \code{start} argument accepts a vector with named
  elements according to the column labels from the
  \code{\link{model.matrix}}.At the time being, initial values for 
  adjacent-categories models must be transformed into the 
  baseline-category model form.

  The argument \code{optim} may be used to modify the optimizer. For
  example, the Nelder-Mead algorithm of \code{\link{optim}} function can
  be called by \code{olmm(..., optim = list(method = "Nelder-Mead"))} or
  the \code{\link{nlminb}} function by \code{olmm(..., optim =
  list(method = "nlminb")).}

  \code{weights} should be used with caution. \code{olmm} can handle
  time-invariant weights only (one distinct weight per subject).
  
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

\seealso{\code{\link{olmm_control}}, \code{\link{ordered}}, \code{\link{olmm-class}}}  

\examples{
## ------------------------------------------------------------------- #
## Example 1: Schizophrenia
##
## Estimating the cumulative mixed models of
## Agresti (2010) chapters 10.3.1 (random intercept
## model) and 10.3.3 (random slope model)
## ------------------------------------------------------------------- #

data(schizo)

## model with random intercepts
model.10.3.1 <-
  olmm(imps79o ~ tx + sqrt(week) + tx * sqrt(week) + re(1|id),
       data = schizo, family = cumulative())

summary(model.10.3.1)

## model with random slopes
model.10.3.3 <-
  olmm(imps79o ~ tx + sqrt(week) + tx * sqrt(week) + re(1 + sqrt(week)|id),
       data = schizo, family = cumulative())

summary(model.10.3.3)


## ------------------------------------------------------------------- #
## Example 2: Movie critics
##
## Estimating three of several adjacent-categories
## mixed models of Hartzl et. al. (2001)
## ------------------------------------------------------------------- #

data(movie)

## model with global random intercepts
model.22.1 <- olmm(critic ~ review + re(1|movie),
                   data = movie, family = adjacent())

summary(model.22.1)


## model with category-specific random intercepts
model.22.2 <- olmm(critic ~ review + re(1|movie, intercept = "ce"),
                   data = movie, family = adjacent())

summary(model.22.2)

## as model 22.1 but with category-specific effects for "review"
model.24.1 <- olmm(critic ~  ce(review) + re(1|movie, intercept = "ce"),
                   data = movie, family = adjacent())

summary(model.24.1)
}

\keyword{models}