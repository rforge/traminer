## --------------------------------------------------------- #
## Author:        Reto Buergin, rbuergin@gmx.ch
## Date:          2013-09-18
##
## Description:
## Routine for fitting different types of linear mixed model
## for ordinal responses with one grouping factor (i.e. 2-level
## models). Random effects are assumed to be multivariate normal
## distributed with expectation 0. At the moment, cumulative
## models with the logit, probit or cauchy link, the
## baseline-category (baseline category = upper category) and
## the adjacent-category logit model are available. The routine
## fits fixed effect and random effect distribution parameters
## based on the marginal likelihood (see Hedeker; 1994), random
## effects are predicted by their expectation (see Hartzl; 2001).
## Predictor-invariable effects can be incorporated by adding a
## second part in the model formula i.e.
## y ~ time + (1 + time | id) | x2.
## Standard deviations of parameter estimates base on the
## expected Fisher-information matrix. It may be reasonable
## use Likelihood profile intervals or the bootstrap method instead.
## Further (partly) required functions can be found in olmm-utils.R,
## olmm-methods.R and import.R
##
## References:
## ordinal:     http://cran.r-project.org/web/packages/ordinal/index.html
## lme4:        http://cran.r-project.org/web/packages/lme4/index.html
## matrixcalc:  http://cran.r-project.org/web/packages/matrixcalc/index.html
## statmod:     http://cran.r-project.org/web/packages/statmod/index.html
##
## Dependencies:
## ucminf:      http://cran.r-project.org/web/packages/ucminf/index.html
##
##
## Modifications:
## 2013-09-15: Free() commands were added in olmm.c
## 2013-09-07: C implementation for updating the marginal Likelihood
## 	       and predicting random-effects was stabilized by 
##	       replacing many alloca's by the R built-in function
##	       Calloc, which may slow the estimation
## 2013-07-27: change 'start' handling and add 'restricted'
##             argument
## 2013-07-19: correct use of numGrad argument (from now the slots
##             score_sbj and score_obs remain empty)
## 2013-07-12: improve use of contrasts. There were irritating
##             warnings under correct use and now the slot
##             'contrasts' also contains contrasts from the
##             model matrix for random effects
##
## To do:
## - implement further family options
## - find better initial parameter values (see polr.R)
## - extract covariance matrix directly from optimizer
## - standardized coefficients
## - unconstrained covariance-matrix for random-effects
## --------------------------------------------------------- #

olmm <- function(formula, data, weights, start, subset, na.action,
                 offset, contrasts, 
                 family = c("cumulative", "baseline", "adjacent"),
                 link = c("logit", "probit", "cauchy"),
                 doFit = TRUE, optim = list(),
                 numGrad = FALSE, numHess = numGrad, nGHQ = 7L,
                 restricted, verbose = FALSE, ...) {

  ## check arguments
  
  if (verbose) cat("* checking arguments ... ")
  mc <- match.call(expand.dots = FALSE)
  if (missing(start)) start <- NULL
  if (missing(restricted)) restricted <- NULL
  
  ## family
  family <- switch(match.arg(family), cumulative = 1L,
                   baseline = 2L, adjacent = 3L)

  ## numerical methods for score and hessian matrix
  if (!numHess & numGrad)
    stop("numHess must be TRUE if numGrad is TRUE")
  
  ## link function
  link <- switch(match.arg(link), logit = 1L, probit = 2L, cauchy = 3L)

  ## check for conflicts between link and family arguments
  if ((family != 1L) & (link != 1L))
    stop("for selected family the logit link is available only")
  
  ## formula
  if (missing(formula)) stop("model requires a formula")

  ## contrasts
  if (missing(contrasts)) contrasts <- NULL

  ## number of Gauss-Hermite quadrature points
  if (nGHQ != as.integer(round(nGHQ)))
    warning(paste("'nGHQ' has been set to ", nGHQ, ".", sep = ""))
  nGHQ <- as.integer(round(nGHQ))
  
  ## optimizer control option
  if (!is.list(optim)) stop("argument optim must be a list")
  optim <- olmm_optim_setup(x = optim, numGrad = numGrad,
                            env = environment())
  
  ## extract model frames
  
  if (verbose) cat("OK\n* extracting model frames ... ")
  
  ## decompose model formula
  form <- olmm_formula(eval.parent(mc$formula), env = parent.frame(n = 1L))
      
  ## set full model frame
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mc), 0)
  mf <- mc[c(1L, m)] # mf = model frame
  mf$formula <- form$full # Form = formula
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  fixefmf <- ranefmf <- mf
  fullmf <- eval.parent(mf)

  ## extract responses
  y <- model.response(fullmf)
  if (!is.ordered(y)) stop("response must be an ordered factor")
  if (nlevels(y) < 2L)
    stop("response variable has less than two response categories")
  
  ## extract fixed effect model matrix
  fixefmf$formula <- terms(form$fixefEtaVar, keep.order = TRUE)
  fixefmfEtaVar <- eval.parent(fixefmf)
  fixefmf$formula <- terms(form$fixefEtaInv, keep.order = TRUE)
  fixefmfEtaInv <- eval.parent(fixefmf)
  
  X <- olmm_mergeMm(x = model.matrix(terms(fixefmfEtaVar), fullmf, contrasts[intersect(names(contrasts), all.vars(form$fixefEtaVar))]), y = model.matrix(terms(fixefmfEtaInv), fullmf, contrasts[intersect(names(contrasts), all.vars(form$fixefEtaInv))]), TRUE)
  rownames(X) <- rownames(fullmf)

  X <- olmm_checkMm(X)
  if (ncol(X) < 1L) stop("design matrix for fixed effects is empty")

  ## check for intercept term
  if (attr(X, "assign")[1L] != 0) {
    intVar <- attr(terms(form$fixefEtaVar), "term.labels")[1L]
    if (!is.factor(fullmf[, intVar])) {
      stop("intercept can be deleted only if the first term in the variable predictor is a factor variable")
    }
    intTerms <- colnames(X)[attr(X, "assign") == 1L & attr(X, "merge") == 1L]
  } else {
    intVar <- colnames(X)[1L]
    intTerms <- colnames(X)[1L]
  }
  
  ## extract random effect grouping factor subject
  hasRanef <- TRUE
  if (is.null(form$subjectName)) { # hack to permit models without random effects
    form$subjectName <- "id"
    form$ranefEtaInv <- formula(~ 1)
    fullmf$id <- factor(1:nrow(fullmf))
    start["ranefCholFac1"] <- 0
    restricted <- unique(c(restricted, "ranefCholFac1"))
    nGHQ <- 1L
    hasRanef <- FALSE
  }
  
  subject <- fullmf[, form$subjectName, drop = TRUE]
  if (!is.factor(subject)) stop("subject variable must be a factor")
  if (nlevels(subject) < 3L)
    stop("subject variable must have 3 levels or more")
  
  ## extract random effect model matrix W
  ranefmf$formula <- terms(form$ranefEtaVar, keep.order = TRUE)
  ranefmfEtaVar <- eval.parent(ranefmf)
  ranefmf$formula <- terms(form$ranefEtaInv, keep.order = TRUE)
  ranefmfEtaInv <- eval.parent(ranefmf)

  contrasts[intersect(names(contrasts), attr(terms(fixefmfEtaVar), "term.labels"))]
  
  W <- olmm_mergeMm(x = model.matrix(terms(ranefmfEtaVar), fullmf, contrasts[intersect(names(contrasts), all.vars(form$ranefEtaVar))]), y = model.matrix(terms(ranefmfEtaInv), fullmf, contrasts[intersect(names(contrasts), all.vars(form$ranefEtaInv))]), FALSE)
  rownames(W) <- rownames(fullmf)
  W <- olmm_checkMm(W)

  ## contrasts
  cons <- append(attr(X, "contrasts"), attr(W, "contrasts"))
  if (!is.null(cons)) cons <- cons[!duplicated(names(cons))]
  if (is.null(cons)) storage.mode(cons) <- "list"
  
  ## vector for dimensions etc.
  dims <- as.integer(c(n = nrow(X), N = nlevels(subject), p = (nlevels(y) - 1L) * sum(attr(X, "merge") == 1L) + sum(attr(X, "merge") == 2L), pEta = ncol(X), pInt = length(intTerms), pEtaVar = sum(attr(X, "merge") == 1L), pEtaInv = sum(attr(X, "merge") == 2L), q = (nlevels(y) - 1L) * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L), qEta = ncol(W), qEtaVar = sum(attr(W, "merge") == 1L), qEtaInv = sum(attr(W, "merge") == 2L), J = nlevels(y), nEta = nlevels(y) - 1L, nPar = (nlevels(y) - 1L) * sum(attr(X, "merge") == 1L) + sum(attr(X, "merge") == 2L) + ((nlevels(y) - 1L) * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L)) * (1L + (nlevels(y) - 1L) * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L)) / 2L, nGHQ = nGHQ, nQP = nGHQ^((nlevels(y) - 1L) * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L)), family = family, link = link, verbose = verbose, numGrad = numGrad, numHess = numHess, doFit = doFit, hasRanef = hasRanef))
  names(dims) <- c("n", "N", "p", "pEta", "pInt", "pEtaVar", "pEtaInv", "q", "qEta", "qEtaVar", "qEtaInv", "J", "nEta", "nPar", "nGHQ", "nQP", "family", "link", "verb", "numGrad", "numHess", "doFit", "hasRanef")

  ## parameter names
  parNames <- list(fixef = c(paste("Eta", rep(seq(1L, dims["nEta"], 1L), each = dims["pEtaVar"]), ":", rep(colnames(X)[attr(X, "merge") == 1L], dims["nEta"]), sep = ""), colnames(X)[attr(X, "merge") == 2L]), ranefCholFac = paste("ranefCholFac", 1L:(dims["q"] * (dims["q"] + 1L) / 2L ), sep = ""))
  
  ## weights
  if (is.null(model.weights(fullmf))) {
    weights <- as.double(rep(1.0, dims["n"]))
    weights_sbj <- as.double(rep(1.0, dims["N"]))
  } else { 
    weights <- model.weights(fullmf)
    if (sum(weights) != dims["n"]) {
      warning("sum of weights must be equal the number of observation. Apply auto correction.")
      weights <- weights / sum(weights) * dims["n"]
    }
    weights_sbj <- tapply(weights, subject, unique)
    if (is.list(weights_sbj)) {
      stop("weights must be constant for subjects")
    } else {
      weights_sbj <- as.double(weights_sbj)
    }
    
    if (length(weights_sbj) != dims["N"]) {
      stop("weights must be constant for subjects")
    }
    if (any(weights < 0.0)) stop("negative weights are not allowed")
  }

  ## offset
  if (is.null(model.offset(fullmf))) {
    offset <- as.double(rep(0, dims["n"]))
  } else {
    offset <- as.double(model.offset(fullmf))
  }

  ## setting the transformed random effect matrix
  u <- matrix(0, dims["N"], dims["q"],
              dimnames = list(levels(subject),
                colnames(parNames$ranefCholFac)))
  
  ## weights and nodes for the Gauss-Hermite quadrature integration
  if (hasRanef) {
    gh <- gauss.quad(nGHQ, "hermite")
    ghx <- olmm_expandQP(gh$nodes, dims["q"]) # with correction
    ghw <- olmm_expandQP(gh$weights * 1 / sqrt(2 * pi) * exp((gh$nodes^2) / 2), dims["q"])
  } else {
    ghx <- matrix(0, 1, 1)
    ghw <- matrix(1, 1, 1)
  }
  
  ## elimination matrix for lower triangular matrices
  ranefElMat <- L.matrix(n = dims["q"])
    
  ## Likelihood function
  logLik_sbj <- rep(0.0, dims["N"])
  names(logLik_sbj) <- levels(subject)
  logLik <- 0.0 # total marginal log Likelihood
  
  ## score function
  if (numGrad == 0L) {
    score_obs <- matrix(0, dims["n"], dims["nPar"])
    rownames(score_obs) <- rownames(X)
    colnames(score_obs) <- unlist(parNames)
    score_sbj <- matrix(0, dims["N"], dims["nPar"],
                        dimnames = list(levels(subject), unlist(parNames)))
  } else {
    score_obs <- matrix(, 0L, 0L)
    score_sbj <- matrix(, 0L, 0L)
  }
  score <- rep(0, dims["nPar"])
  names(score) <- unlist(parNames)

  ## info matrix
  info <- matrix(0, dims["nPar"], dims["nPar"],
                 dimnames = list(unlist(parNames), unlist(parNames)))
  
  ## linear predictor (without contributions of random effects)
  eta <- matrix(0, dims["n"], dims["nEta"],
                dimnames = list(rownames(fullmf),
                  paste("Eta", 1L:(dims["J"] - 1L), sep = "")))

  ## inital values
  
  if (verbose) cat("OK\n* setting inital values ... ")

  start <- olmm_start(start, dims, parNames, X, W, eta, ranefElMat)

  ## restricted
  restr <- rep(FALSE, dims["nPar"])
  names(restr) <- unlist(parNames)
  if (!is.null(restricted)) {

    ## checks
    stopifnot(is.character(restricted))
    if (dims["family"] == 3)
      stop("'restricted' argument is not available for adjacent category model")
    if (!all(restricted %in% names(restr)))
      stop(paste("the coefficient(s) ", paste("'", restricted[!restricted %in% names(restr)], "'", sep = "", collapse = ", "), " in 'restricted' were not found. The coefficient names are ", paste("'", names(restr), "'", sep = "", collapse = ", "), ".", sep = ""))
    
    ## set restricted parameters
    restr[restricted] <- TRUE
  }
  
  ## xlevels
  xlevels <- .getXlevels(attr(fullmf, "terms"), fullmf)
  if (is.null(xlevels)) storage.mode(xlevels) <- "list"

  ## terms
  terms <- list(fixefEtaVar = terms(form$fixefEtaVar, keep.order = TRUE),
                fixefEtaInv = terms(form$fixefEtaInv, keep.order = TRUE),
                ranefEtaVar = terms(form$ranefEtaVar, keep.order = TRUE),
                ranefEtaInv = terms(form$ranefEtaInv, keep.order = TRUE))
  
  ## define fit object
  
  if (verbose) cat("OK\n* building the model object ... ")
  
  object <- new(Class = "olmm",
                call = mc,
                frame = fullmf,
                formula = formula,
                terms = terms,
                y = y,
                X = X,
                W = W,
                subject = subject,
                subjectName = form$subjectName,
                weights = weights,
                weights_sbj = weights_sbj,
                offset = offset,
                xlevels = xlevels,
                contrasts = cons,
                dims = dims,
                fixef = start$fixef,
                ranefCholFac = start$ranefCholFac,
                coefficients = start$coefficients,
                restricted = restr,
                eta = eta,
                u = u,
                logLik_sbj = logLik_sbj,
                logLik = logLik,
                score_obs = score_obs,
                score_sbj = score_sbj,
                score = score,
                info = info,
                ghx = ghx,
                ghw = ghw,
                ranefElMat = ranefElMat,
                optim = optim)

  ## delete big data blocks
  rm(list = ls()[!ls() %in% c("doFit", "object", "verbose")])
  
  if (doFit) {

    ## set fitting evironment
    
    if (verbose) cat("OK\n* setting up the fitting environment ... ")
  
    ## set start parameters
    object@optim[[1L]] <- object@coefficients 
    object@optim[[4L]] <- object@restricted
    ## fit the model
    
    if (verbose) cat("OK\n* fitting the model ... ")

    if (!is.null(object@optim$control$trace) &&
        object@optim$control$trace > 0)
      cat("\n")

    ## extract the function for fitting the model
    if (object@optim$method %in% c("nlminb", "ucminf")) {
      FUN <- object@optim$method
      object@optim <- object@optim[-which(names(object@optim) == "method")] 
    } else {
      FUN <- "optim"
    }
    
    systemTime <- system.time(object@output <- do.call(FUN, object@optim))

    if (FUN %in% c("nlminb", "ucminf")) object@optim$method <- FUN
    
    ## to get sure ...
    .Call("olmm_update_marg", object, object@output$par, PACKAGE = "vcolmm")
      
    ## print messages for opimization
    if (verbose) {
      cat(paste("OK\n\toptimization time:",
                signif(systemTime[3L], 3L),
                "seconds", sep = " "))
      if (is.null(object@output$message)) {
        cat("\n\tno message returned by the optimizer")
      } else {
        cat(paste("\n\tmessage: ", object@output$message, sep = ""))
      }
    }

    ## warnings from optimization
    olmm_optim_warnings(object@output, FUN)
    
    ## numeric estimate of fisher information
    if (object@dims["numHess"] == 1L) {
      
      if (verbose)
        cat("\n* computing the approximative hessian matrix ... ")
      
      object@info[] <- # replace the info slot
        - hessian(func = object@optim[[2L]], x = object@coefficients,
                  method.args = list(r = 2, show.details = TRUE))
      if (verbose) cat("OK")
    }

    if (verbose) {
      eigenHess <- eigen(object@info, only.values = TRUE)$values
      condHess <- abs(max(eigenHess) / min(eigenHess))
      cat("\n\tcondition number of Hessian matrix:",
          format(condHess, digits = 2L, scientific = TRUE))
    }  
    
    ## fit / predict random effects   
    
    ## compute expected standardized random effects
    if (object@dims["hasRanef"] > 0) {
      if (verbose) cat("\n* predicting random effects ... ")
      .Call("olmm_update_u", object, PACKAGE = "vcolmm")
    }
    
    ## reset environment of estimation equations
    environment(object@optim[[2L]]) <- baseenv()
    if (!object@dims["numGrad"]) environment(object@optim[[3]]) <- baseenv()
    
    if (verbose) cat("OK\n* computations finished, return model object\n")
    
  } else {

    ## update the object with the current estimates
    .Call("olmm_update_marg", object, object@coefficients,
          PACKAGE = "vcolmm")
    if (object@dims["hasRanef"] > 0L)
      .Call("olmm_update_u", object, PACKAGE = "vcolmm")

    if (verbose) cat("\n* no computations processed, return model object\n")
  } 
  
  return(object) 
}
