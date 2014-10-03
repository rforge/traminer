\name{tvcm}

\alias{tvcm}
\alias{tvcolmm}
\alias{tvcglm}

\title{Tree-based varying coefficient regression models}

\description{
  TVCM is a tree-based algorithm that aims to estimate varying
  coefficient regression models. TVCM approximates varying coefficients
  by piecewise constant functions, i.e., it estimates linear models with
  stratum specific coefficients. The \command{\link{tvcm}} function
  implements the partitioning algorithms described in Buergin and
  Ritschard (2014b) (default) and Buergin and Ritschard (2014a).
}

\usage{
tvcm(formula, data, fit, family,
     weights, subset, offset, na.action,
     control = tvcm_control(), ...)

tvcolmm(formula, data, family = cumulative(), 
        weights, subset, offset, na.action, 
        control = tvcm_control(), ...)

tvcglm(formula, data, family, 
       weights, subset, offset, na.action, 
       control = tvcm_control(), ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to fit, e.g.,

    \code{y ~ vc(z1, z2) + vc(z1, z2, by = x)}
    
    where \code{vc} specifies the varying coefficients. See 
    \code{\link{vcrpart-formula}}.}
  \item{fit}{a character string or a function that specifies the fitting
    function, e.g., \command{\link{olmm}} or \command{\link{glm}}.}  
  \item{family}{the model family, e.g., an object of class 
    \command{\link{family.olmm}} or \command{\link{family}}.} 
  \item{data}{a data frame containing the variables in the model.}
  \item{weights}{an optional numeric vector of weights to be used in the 
    fitting process.}
  \item{subset}{an optional logical or integer vector specifying a
    subset of \code{'data'} to be used in the fitting process.}
  \item{offset}{this can be used to specify an a priori known component 
    to be included in the linear predictor during fitting.} 	
  \item{na.action}{a function that indicates what should happen if data 
    contain \code{NA}s. See \command{\link{na.action}}.}
  \item{control}{a list with control parameters as returned by
    \command{\link{tvcm_control}}.}
  \item{\ldots}{additional arguments passed to the fitting function
    \code{fit}.}  
}

\details{ 
  TVCM partitioning works as follows: In each iteration we fit the
  current model and select a binary split for one of the current
  terminal nodes. The selection requires 4 decisions: the \code{vc}
  term, the node, the variable and the cutpoint in the selected
  variable. The algorithm starts with \eqn{M_k = 1} node for each for
  each of the \eqn{K} \command{\link{vc}} terms and iterates until the
  criteria defined by \code{control} are reached, see
  \command{\link{tvcm_control}}. The following two paragraphs explain
  the two available procedures to select the split and choose the final
  model. 

  The first and default procedure (cf. Buergin and Ritschard, 2014b)
  uses two stages. The first stage builds too overly fine partitions and
  the second stage selects the best-sized partitions by pruning. For the
  second stage, which is automatically processed, we refer to
  \command{\link{tvcm-assessment}}. The partitioning stage selects, in
  each iteration, the split that maximizes the training error reduction
  compared to the current model. The training function is defined by the
  loss function specified by the \code{lossfun} argument in
  \command{\link{tvcm_control}}. It is also possible to penalize the
  training error, which may be interesting if \code{vc} terms include
  different number of coefficients. See the  arguments \code{cp},
  \code{dfpar} and \code{dfsplit} in \command{\link{tvcm_control}} to
  define such a penalization. The default stopping parameters are
  \code{minsize = 30} (a minimum node size of 30) and \code{mindev = 2}
  (the training error reduction of the best split must be larger than
  two to continue). 

  The second procedure (cf. Buergin and Ritschard, 2014a) is, for
  technical reasons, restricted to the cases where a single
  \command{\link{vc}} term is used and non of the moderators
  (partitioning variables) intersects with predictors. On the other
  hand, the second procedure avoids the variable selection bias that
  exhibits the exhaustive search above and is computationally much less
  burdensome. The procedure selects the split by first choosing the node
  and the variable with M-fluctuation tests (cf. Zeileis and Hornik,
  2007) and second choosing the cutpoint by the training error reduction
  statistic, as above. To use this option it is necessary to set
  \code{sctest=TRUE} in \command{\link{tvcm_control}}. The algorithm is
  stopped as soon as all nodewise Bonferroni corrected p-values of
  M-fluctuation tests reach a prespecified threshold (e.g., 0.05), see
  argument \code{alpha} in \command{\link{tvcm_control}} and no pruning
  is necessary (or optional). Note that, as explained in Buergin and
  Ritschard (2014a), coefficient constancy tests are adjusted for
  intra-subject correlation for 2-stage models,  see
  \command{\link{estfun.olmm}}. The procedure is illustrated below in
  example 2. 

  An alternative tree-based algorithm to \command{\link{tvcm}} are the
  MOB (Zeileis et al., 2008) and the PartReg (Wang and Hastie, 2014)
  algorithms. The MOB algorithm is implemented by the \command{mob}
  function in the packages \pkg{party} and \pkg{partykit}. For
  alternative, smoothing splines and kernel regression approaches to
  varying coefficients, see the packages \pkg{mgcv},
  \pkg{svcm},\pkg{mboost} or \pkg{np}. 

  The \command{\link{tvcm}} function builds on the software
  infrastructure of the \pkg{partykit} package. The authors are grateful
  for these codes. 
}

\value{An object of class \command{\link{tvcm}}. The
    \command{\link{tvcm}} class itself is based on the
    \command{\link{party}} class of the \pkg{partykit} package. The most
    important slots are: \item{node}{an object of class
    \command{\link{partynode}}.} \item{data}{a (potentially empty)
    \command{\link{data.frame}}.} \item{fitted}{an optional
    \command{\link{data.frame}} with \code{nrow(data)} rows and
    containing at least the fitted terminal node identifiers as element
    \code{(fitted)}. In addition, weights may be contained as element
    \code{(weights)} and responses as \code{(response)}.}
    \item{info}{additional information including\code{control} and
    \code{model}.} 
}

\references{ 
   Zeileis, A., Hothorn, T., and Hornik, K. (2008). Model-Based
   Recursive Partitioning. \emph{Journal of Computational and Graphical
   Statistics}, \bold{17}(2), 492--514. 

   Zeileis, A., Hornik, K. (2007), Generalized M-Fluctuation Tests for
   Parameter Instability, \emph{Statistica Neerlandica}, \bold{61},
   488--508. doi:10.1111/j.1467-9574.2007.00371.x. 

   Torsten Hothorn, Achim Zeileis (2014). partykit: A Modular Toolkit
   for Recursive Partytioning in R. Working Paper 2014-10. Working
   Papers in Economics and Statistics, Research Platform Empirical and
   Experimental Economics, Universitaet Innsbruck. URL
   \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2014-10} 

   Wang, J. C., Hastie, T. (2013), Boosted Varying-Coefficient
   Regression Models for Product Demand Prediction, \emph{Journal of
   Computational and Graphical Statistics}.

   Buergin R. and Ritschard G. (2014a), Tree-based varying coefficient
   regression for longitudinal ordinal responses. Submitted article. 

   Buergin R. and Ritschard G. (2014b), Coefficient-wise tree-based
   varying coefficient regression with vcrpart. Article in progress.

   Akaike H. (1974), A new look at the statistical model identification,
   \emph{IEEE Transactions on Automatic Control}, \bold{19}, 716--723. 
}

\seealso{\command{\link{tvcm_control}}, \command{\link{tvcm-methods}},
  \command{\link{tvcm-plot}}, \command{\link{tvcm-plot}},
  \command{\link{tvcm-assessment}}}

\examples{

## ------------------------------------------------------------------- #  
## Example 1: Moderated effect of education on poverty
##
## The algorithm is used to find out whether the effect of high
## education 'EduHigh' on poverty 'Poor' is moderated by the civil
## status 'CivStat'. We specify two 'vc' terms in the logistic
## regression model for 'Poor': a first that accounts for the direct
## effect of 'CivStat' and a second that accounts for the moderation of
## 'CivStat' on the relation between 'EduHigh' and 'Poor'. We use here
## the 2-stage procedure with a partitioning- and a pruning stage as
## described in Buergin and Ritschard (2014b). 
## ------------------------------------------------------------------- #

data(poverty)
poverty$EduHigh <- 1 * (poverty$Edu == "high")

## fit the model
model.Pov <-
  tvcglm(Poor ~ -1 +  vc(CivStat) + vc(CivStat, by = EduHigh) + NChild, 
         family = binomial(), data = poverty, subset = 1:200,
         control = tvcm_control(verbose = TRUE, papply = lapply,
           folds = folds_control(K = 1, type = "subsampling", seed = 7)))

## diagnosis
plot(model.Pov, "cv")
plot(model.Pov, "coef")
summary(model.Pov)
splitpath(model.Pov, steps = 1:3)
prunepath(model.Pov, steps = 1)


## ------------------------------------------------------------------- # 
## Example 2: Moderated effect effect of unemployment
##
## Here we fit a varying coefficient ordinal linear mixed on the 
## synthetic ordinal longitudinal data 'unemp'. The interest is whether 
## the effect of unemployment 'UNEMP' on happiness 'GHQL' is moderated 
## by 'AGE', 'FISIT', 'GENDER' and 'UEREGION'. 'FISIT' is the only true  
## moderator. For the the partitioning we coefficient constancy tests,
## as described in Buergin and Ritschard (2014a)
## ------------------------------------------------------------------- #

data(unemp)

## fit the model
model.UE <-
  tvcolmm(GHQL ~ -1 + 
          vc(AGE, FISIT, GENDER, UEREGION, by = UNEMP, intercept = TRUE) +
          re(1|PID), data = unemp, control = tvcm_control(sctest = TRUE))

## diagnosis (no cross-validation was performed since 'sctest = TRUE')
plot(model.UE, "coef")
summary(model.UE)
splitpath(model.UE, steps = 1, details = TRUE)


}

\author{Reto Buergin}

\keyword{tree}