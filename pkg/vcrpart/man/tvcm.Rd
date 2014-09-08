\name{tvcm}

\alias{tvcm}
\alias{tvcolmm}
\alias{tvcglm}

\title{Tree-based varying coefficient regression models}

\description{
  TVCM is a tree-based algorithm that aims to estimate varying coefficient 
  regression models. TVCM approximates varying coefficients by piecewise
  constant functions, i.e., it estimates linear models with stratum
  specific coefficients. The \code{\link{tvcm}} function implements the
  partitioning algorithms described in Buergin and Ritschard (2014b)
  (default) and Buergin and Ritschard (2014a).
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
  \item{fit}{a character string or a function that specifies the 
    fitting function, e.g., \code{\link{olmm}} or \code{\link{glm}}.} 
  \item{family}{the model family, e.g., an object of class 
    \code{\link{family.olmm}} or \code{\link{family}}.} 
  \item{data}{a data frame containing the variables in the model.}
  \item{weights}{an optional numeric vector of weights to be used in the 
    fitting process.}
  \item{subset}{an optional logical or integer vector specifying a subset
    of \code{'data'} to be used in the fitting process.}
  \item{offset}{this can be used to specify an a priori known component 
    to be included in the linear predictor during fitting.} 	
  \item{na.action}{a function that indicates what should happen if data 
    contain \code{NA}s. See \code{\link{na.action}}.}
  \item{control}{a list with control parameters as returned by
    \code{\link{tvcm_control}}.}
  \item{\ldots}{additional arguments passed to the fitting function
    \code{fit}.}  
}

\details{ 
  The TVCM partitioning algorithm works as follows: Starting with
  \eqn{M_k = 1} stratum (i.e. node) for all \eqn{K} \code{\link{vc}}
  terms, the algorithm splits in each iteration one of the current 
  \eqn{K \sum_{k=1}^K M_k} nodes into two new nodes. For selecting 
  the \code{vc} term, the node, the variable and the cutpoint in each
  iteration, there are two procedures available.

  The first and default procedure (cf. Buergin and Ritschard, 2014b) 
  is a two-stage procedure which builds overly fine partitions in the
  first stage and selects the best-sized partitions by pruning in the
  second stage. For the second stage, which is automatically processed, we 
  refer to \code{\link{tvcm-assessment}}. The partitioning stage selects, 
  in each iteration, the split that maximizes the penalized loss reduction 
  statistic that compares the current with the one-step ahead model 
  (see also argument \code{lossfun} in \code{\link{tvcm_control}}). By  
  default, the penalized loss reduction statistic is the Akaike Information 
  Criterion (Akaike, 1974), but alternatives could be specified by the 
  \code{control} argument. The algorithm is continued to build until the 
  criteria specified in \code{control} are reached. By default, only the 
  \code{minsize} (minimum node size) criteron are specified (and 
  that the penalized loss reduction is larger than 0). For large data sets 
  with many partitioning variables, this can be slow and, therefore, 
  \code{\link{tvcm_control}} provides further criteria. In particular, you may 
  increase the minimum permitted penalized loss reduction by \code{dfsplit}.

  The second procedure (cf. Buergin and Ritschard, 2014a) is, for technical 
  reasons, restricted to the cases where a single \code{\link{vc}} term is used 
  and non of the moderators (partitioning variables) intersects with predictors. 
  On the other hand, the second procedure avoids the variable selection bias 
  that exhibits the exhaustive search above and is computationally much less 
  burdensome. The procedure selects the split by first choosing the node and 
  the variable with M-fluctuation tests (cf. Zeileis and Hornik, 2007) and second 
  choosing the cutpoint by the deviance reduction statistic, as above. To use this 
  option it is necessary to set \code{sctest=TRUE} in \code{\link{tvcm_control}}. 
  The algorithm is stopped as soon as all nodewise Bonferroni corrected p-values 
  of M-fluctuation tests reach a prespecified threshold (e.g., 0.05), see 
  argument \code{alpha} in \code{\link{tvcm_control}} and no pruning is 
  necessary (or optional). Note that, as explained in Buergin and Ritschard (2014a),
  coefficient constancy tests are adjusted for intra-subject correlation for 
  2-stage models,  see \code{\link{estfun.olmm}}. The procedure is illustrated below 
  in example 2.

  An alternative algorithm to \code{\link{tvcm}} is the MOB algorithm of Zeileis 
  et al. (2008), which is implemented by the \code{mob} function in the 
  packages \pkg{party} (url: \url{http://cran.r-project.org/web/packages/party/}) and 
  \pkg{partykit} (url: \url{http://cran.r-project.org/web/packages/party/}).

  The \code{\link{tvcm}} function builds on the software infrastructure of the
  \pkg{partykit} package. The authors are grateful for these codes
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

   Torsten Hothorn, Achim Zeileis (2014). partykit: A Modular Toolkit for Recursive
   Partytioning in R. Working Paper 2014-10. Working Papers in Economics and Statistics,
   Research Platform Empirical and Experimental Economics, Universitaet Innsbruck. URL 
   \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2014-10}

   Wang, J. C., Hastie, T. (2013), Boosted Varying-Coefficient Regression
   Models for Product Demand Prediction, \emph{Journal of Computational
   and Graphical Statistics}.

   Buergin R. and Ritschard G. (2014a), Tree-based varying coefficient 
   regression for longitudinal ordinal responses. Submitted article.

   Buergin R. and Ritschard G. (2014b), Coefficient-wise tree-based varying
   coefficient regression with vcrpart. Article in progress.

   Akaike H. (1974), A new look at the statistical model identification, 
   \emph{IEEE Transactions on Automatic Control}, \bold{19}, 716--723.
}

\seealso{\code{\link{tvcm_control}}, \code{\link{tvcm-methods}},
  \code{\link{tvcm-plot}}, \code{\link{tvcm-plot}},
  \code{\link{tvcm-assessment}}}

\examples{

## ------------------------------------------------------------------- # 
## Example 1: Moderated effect of education on poverty
##
## The algorithm is used to find out whether the effect of high education 
## 'EduHigh' on poverty 'Poor' is moderated by the civil status 'CivStat'.
## We specify two 'vc' terms in the logistic regression model for 'Poor': 
## a first that accounts for the direct effect of 'CivStat' and a second 
## that accounts for the moderation of 'CivStat' on the relation between 
## 'EduHigh' and 'Poor'. We use here the 2-stage procedure with
## a partitioning- and a pruning stage as described in 
## Buergin and Ritschard (2014b).
## ------------------------------------------------------------------- #

data(poverty)

## fit the model
poverty$EduHigh <- 1 * (poverty$Edu == "high")
model.Pov <- tvcglm(Poor ~ -1 +  vc(CivStat) + vc(CivStat, by = EduHigh) + NChild, 
                    family = binomial(), data = poverty, subset = 1:200,
                    control = tvcm_control(verbose = TRUE,
                        folds = folds_control(K = 1, type = "subsampling", seed = 4)))

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
model.UE <- tvcolmm(GHQL ~ -1 + 
                    vc(AGE, FISIT, GENDER, UEREGION, by = UNEMP, intercept = TRUE) +
                    re(1|PID), data = unemp, control = tvcm_control(sctest = TRUE))

## diagnosis (no cross-validation was performed since 'sctest = TRUE')
plot(model.UE, "coef")
summary(model.UE)
splitpath(model.UE, steps = 1, details = TRUE)

}

\author{Reto Buergin}

\keyword{tree}