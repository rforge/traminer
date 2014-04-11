## --------------------------------------------------------- #
## Author:          Reto Buergin, rbuergin@gmx.ch
## Date:            2014-03-27
##
## Description:
## Utility functions for the olmm function (see olmm.R). Some
## functions are experimental and not listed in the namespace.
##
## References and dependencies:
## statmod:         http://cran.r-project.org/web/packages/statmod/index.html
## ucminf:          http://cran.r-project.org/web/packages/ucminf/index.html
##
## Contents:
##
## Functions for olmm:
## olmm_expandQP:        expand grid for numerical integration
## olmm_fn
## olmm_gn
## olmm_optim_setup:     set up algorithm function
## olmm_optim_warnings:
## olmm_coefShortLabs:   short labels for coefficient names
## olmm_formula:         extract model matrix formula from input formula
## olmm_checkMm:         check and modify model matrix
## olmm_start:           set initial values
## olmm_mergeMm:         merge the predictor-variable and predictor-invariant
##                       model matrices
##                     
## gcauchy:              derivate of dcauchy
## glogis:               derivate of dlogis
## gnorm:                derivate of dnorm
##
##
## Functions for tvcolmm:
## olmm_update2:           updates the model for a new set of coefficients
##
## Functions for estfun.olmm:
## olmm_scoreVar:          computes the variance of the observation scores.
## olmm_scoreCovWin:       computes the intra-subject covariance of
##                         the observation scores.
## olmm_scoreCovBet:       computes the between-subject covariance
##                         of the observation scores.
## olmm_f_scoreTransfMat:  the equation to be optimized to zero:
##                         the difference between the adjusted
##                         intra-subject covariance and the adjusted
##                         between-subject covariance of likelihood scores.
## olmm_g_scoreTransfMat:  the derivation of olmm_fDecorrelate.
## olmm_scoreTransfMat:    computes the transformation matrix for
##                         removing the intra-subject
##                         correlation of ML scores.
## olmm_simPredictors:     extracts randomly observed predictor
##                         combinations
## olmm_addObsScores:      computes scores to complete an unbalance
##                         data set 
##
## To do:
## - replace ranefChol fac initial value with covariance matrix
## - add multiple family arguments
## - make olmm_refit_MC working (compare version 0.1-5)
## --------------------------------------------------------- #


olmm_expandQP <- function(x, q) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Expand quadrature points to the dimension of random
  ## effects.
  ##
  ## Arguments:
  ## x: nodes (or weights) for one dimension.
  ## q: number of random coefficients.
  ##
  ## Value:
  ## A grid matrix of nodes (or weights) of dimension q.
  ## ------------------------------------------------------- #

  if (length(q)^q * q > 2^31 - 1)
    stop("number of quadrature weights is too large. Decrease 'nAGQ'")
  rval <- matrix(0, length(x)^q, q)
  for (i in 1:q) rval[, i] <- x[rep(1:length(x), each = length(x)^(i - 1L), length.out = length(x)^q)]
  return(rval)
}


olmm_fn <- function(par, restricted) {
  if (!exists("object")) object <- new(Class = "olmm")
  parNew <- object@coefficients
  parNew[!restricted] <- par[!restricted]
  .Call("olmm_update_marg", object, as.numeric(parNew), PACKAGE = "vcolmm")
  return(-object@logLik)
}

olmm_gr <- function(par, restricted) {
  if (!exists("object")) object <- new(Class = "olmm")
  scoreNew <- object@score
  scoreNew[restricted] <- 0.0
  return(-scoreNew)
}

olmm_optim_setup <- function(x, numGrad, env = parent.frame()) {  

  ## ------------------------------------------------------- #
  ## Description:
  ## Processes the algorithm argument.
  ##
  ## Arguments:
  ## x:        argument 'optim' of 'olmm' call
  ## numGrad:  argument 'numGrad' of 'olmm' call
  ## env:      environment of the optimization 
  ##
  ## Value:
  ## A prepared list to be assigned to the eval command
  ## envoking the optimization (using the eval function).
  ## ------------------------------------------------------- #
  
  rval <- list(par = NULL,
               fn = olmm_fn,
               gr = if (numGrad) NULL else olmm_gr,
               restricted = NULL)
  
  rval <- append(rval, x)
  
  if (sum(duplicated(names(rval))) > 0) {
    warning(paste("list elements", paste(paste("'", names(rval)[duplicated(names(rval))], "'", sep = ""), collapse = ", "), "of 'optim' are replaced"))
    rval <- rval[!duplicated(names(rval))]
  }

  if (is.null(rval$method)) {
    rval$method <- "ucminf"
  }
  if (!rval$method %in% c("ucminf", "nlminb", eval(formals(optim)$method)))
    stop("algorithm unknown.")

  if (rval$method == "nlminb") {
    names(rval)[1:3] <- names(formals(nlminb)[1:3])
  }
  
  ## set environment for objective function and gradient
  environment(rval[[2]]) <- env
  if (!numGrad) environment(rval[[3]]) <- env
  
  return(rval)
}

olmm_optim_warnings <- function(output, FUN) {

    if (FUN == "optim") {
        switch(as.character(output$convergence),
               "1" = warning("Stopped by small step (xtol)."),
               "20" = warning("Inadmissible initial parameters."),
               "21" = warning("Inadmissible immediate parameters."),
               "10" = warning("Degeneracy of the Nelder-Mead simplex."),
               "51" = warning("Warning from 'L-BFGS-B'."),
               "52" = warning("Error from 'L-BFGS-B'."),
               NULL)
    }

    if (FUN == "nlminb") {
        if (output$convergence != 0) warning(output$message)
    }
    
    if (FUN == "ucminf") {

        switch(as.character(output$convergence),
               "2" =
               if (output$info["maxgradient"] > 1e-3) warning("Stopped by small step (xtol)."),
               "3" = 
               warning("Stopped by function evaluation limit (maxeval)."),
               "4" =
               if (output$info["maxgradient"] > 1e-3) warning("Stopped by zero step from line search."),
               "-2" =
               warning("Computation did not start: length(par) = 0."),
               "-4" =
               warning("Computation did not start: stepmax is too small."),
               "-5" =
               warning("Computation did not start: grtol or xtol <= 0"),
               "-6" =
               warning("Computation did not start: maxeval <= 0"),NULL)
    }
}

olmm_coefShortLabs <- function(object) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Proposes abbreviations for names of coefficients.
  ##
  ## Arguments:
  ## object: an olmm object
  ##
  ## Value:  vector of character strings with abbreviations
  ##         for names of coefficients.
  ## ------------------------------------------------------- #

  ## abbreviations for predictor-variable fixed effects
  abbEtaVar <- names(fixef(object, "predictor-variable"))
  abbEtaVar <- gsub("Eta", "", abbEtaVar, fixed = TRUE)
  abbEtaVar <- gsub("(Intercept)", "I", abbEtaVar, fixed = TRUE)
  abbEtaVar <- abbreviate(abbEtaVar)

  ## abbreviations for predictor-invariant fixed effects
  abbEtaInv <- names(fixef(object, "predictor-invariant"))
  abbEtaInv <- abbreviate(abbEtaInv)

  ## abbreviations for Choleski factors
  q <- object@dims[["q"]]
  abbReCF <- paste("reCF", 1:(q * (q + 1) / 2), sep = "")
  
  return(c(abbEtaVar, abbEtaInv, abbReCF))
}


olmm_formula <- function(formula, env = parent.frame()) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Evaluates the input formula and returns a list with
  ## the necessary parts for constructing model matrices.
  ##
  ## Arguments:
  ## formula: original formula.
  ## env:     environment for model.frame or model.matrix
  ##          evaluation.
  ##
  ## Value:
  ## A list with multiple formulas. For internal use only.
  ## ------------------------------------------------------- #
  
  mc <- match.call()

  if (missing(env)) env <- parent.frame(n = 1)

  ff <- as.Formula(formula)
    
  ## check whether the response is unique
  if (length(ff)[1] != 1L)
    stop("number of left hand formulas must be 1")
  if (length(all.vars(formula(ff, lhs = 1,
                              rhs = rep(FALSE, length(ff)[2])))) > 1)
    stop("only one response variables is allowed")
  yName <- all.vars(formula)[1]
  
  ## check wether the right hand side has 1 or 2 parts
  if (length(ff)[2] > 2L)
    stop("number of right side formula must be 1 or 2")
  
  ## predictor-invariant formula
  
  fTmp1 <- formula(ff, lhs = -1L, rhs = 1L)
  fixefEtaInv <- nobars(fTmp1)
  if (is.null(fixefEtaInv)) {
    fixefEtaInv <- formula(~ 1)
  } else {
    fixefEtaInv <- formula(fixefEtaInv)
  }
  barList <- expandSlash(findbars(fTmp1))
  if (!is.null(barList)) { # if random effects exist
    if (length(barList) > 1L) {
      stop("maximum one grouping factor (subject) is allowed")
    } else {
      barList <- barList[[1]]
    }
    ranefEtaInv <- formula(paste("~", deparse(barList[[2]])))
    subjectName <- deparse(barList[[3L]])
  } else {
    ranefEtaInv <- formula(~ -1)
    subjectName <- NULL
  }
  
  ## predictor-variable formula
  
  if (length(ff)[2] == 2L) {
    
    fTmp1 <- formula(ff, lhs = -1, rhs = 2)
    fixefEtaVar <- nobars(fTmp1)
    if (is.null(fixefEtaVar)) {
      fixefEtaVar <- formula(~ 1)
    } else {
      fixefEtaVar <- formula(fixefEtaVar)
    }
    if (attr(terms(fixefEtaVar, keep.order = TRUE), "intercept") == 0 &&
        length(attr(terms(fixefEtaVar), "term.labels")) == 0) {
      stop("intercept can be deleted only if the first variable of the variable predictor is categorical")
    }
    
    barList <- expandSlash(findbars(fTmp1))      
    if (!is.null(barList)) {
      if (length(barList) > 1L) {
        stop("maximum one grouping factor (subject) is allowed")
      } else {
        barList <- barList[[1]]
      }
      ranefEtaVar <- formula(paste("~", deparse(barList[[2]])))
      if (is.null(subjectName)) {
        subjectName <- deparse(barList[[3L]])
      } else {
        if (subjectName != deparse(barList[[3L]]))
          stop("maximum one grouping factor (subject) is allowed")
      } 
    } else {
      ranefEtaVar <- formula(~ - 1)
    }
  } else {

    fixefEtaVar <- formula(~ 1)
    ranefEtaVar <- formula(~ -1)
  }
  
  fixefEtaVarTerms <-
    attr(terms(fixefEtaVar, keep.order = TRUE), "term.labels")
  fixefEtaInvTerms <-
    attr(terms(fixefEtaInv, keep.order = TRUE), "term.labels")

  ## check whether effects do not appear in both fixed effect equations
  if (length(tmp <- intersect(fixefEtaVarTerms, fixefEtaInvTerms)) > 0) {
    warning("one or more terms appear in both predictor-variable and predictor-invariant fixed effects. Delete the corresponding predictor-invariant effects.")
    fixefEtaInv <-
      update(fixefEtaInv,
             formula(paste("~ . - ", paste(tmp, collapse = " - "))))
  }

  ranefEtaVarTerms <-
    attr(terms(ranefEtaVar, keep.order = TRUE), "term.labels")
  ranefEtaInvTerms <-
    attr(terms(ranefEtaInv, keep.order = TRUE), "term.labels")
  
  ## check whether effects do not appear in both random effect equations
  if (length(tmp <- intersect(ranefEtaVarTerms, ranefEtaInvTerms)) > 0) {
    warning("one or more terms appear in both predictor-variable and predictor-invariant random effects. Delete the corresponding predictor-invariant effects.")
    ranefEtaInv <-
      update(ranefEtaInv,
             formula(paste("~ . - ", paste(tmp, collapse = " - "))))
  }

  if (attr(terms(ranefEtaVar), "intercept") == 1 &
      attr(terms(ranefEtaInv), "intercept") == 1)
    ranefEtaInv <- formula(paste(deparse(ranefEtaInv), "-1"))
  
  ## full formula for model frame
  if (sum(length(all.vars(fixefEtaVar)) + length(all.vars(fixefEtaInv)) +
          length(all.vars(ranefEtaVar)) + length(all.vars(ranefEtaInv)) +
          length(subjectName)) > 0) {
    full <- formula(paste(yName, " ~ ",
                          paste(c(attr(terms(fixefEtaVar, keep.order = TRUE),
                                       "term.labels"),
                                  attr(terms(fixefEtaInv, keep.order = TRUE),
                                       "term.labels"),
                                  attr(terms(ranefEtaVar, keep.order = TRUE),
                                       "term.labels"),
                                  attr(terms(ranefEtaInv, keep.order = TRUE),
                                       "term.labels"),
                                  subjectName), collapse = " + ")))
  } else {
    full <- formula(paste(yName, " ~ 1"))
  }
  
  ## fixef formula (used in predict)
  fixef <- 
    formula(paste(" ~ ",
                  paste(c("-1", # could be deleted???
                          attr(terms(fixefEtaVar, keep.order = TRUE),
                               "term.labels"),
                          attr(terms(fixefEtaInv, keep.order = TRUE),
                               "term.labels")),
                        collapse = " + ")))

  ## ranef
  ranef <-
    formula(paste(" ~ ",
                  paste(c("1",
                          attr(terms(ranefEtaVar, keep.order = TRUE),
                               "term.labels"),
                          attr(terms(ranefEtaInv, keep.order = TRUE),
                               "term.labels")),
                        collapse = " + ")))
  
  ## modify environment of several formulas
  environment(ff) <- environment(full) <-
    environment(fixef) <- 
      environment(fixefEtaVar) <- environment(fixefEtaInv) <-
        environment(ranef) <- 
          environment(ranefEtaVar) <- environment(ranefEtaInv) <- env
  
  ## return several needed formulas
  return(list(original = ff,
              full = full,
              fixef = fixef,
              fixefEtaVar = fixefEtaVar,
              fixefEtaInv = fixefEtaInv,
              ranef = ranef,
              ranefEtaVar = ranefEtaVar,
              ranefEtaInv = ranefEtaInv,
              subjectName = subjectName))  
}

olmm_mergeMm <- function(x, y, deleteIntY = TRUE) {

  ## ------------------------------------------------------- #
  ## Description: Merge a model matrix of predictor variable
  ##              effects with a model matrix of predictor
  ##              invariant effects
  ##
  ## Arguments:
  ## x: a model matrix for predictor-variable effects.
  ## y: a model matrix for predictor-invariant effects. 
  ##
  ## Value:
  ## A model matrix.
  ## ------------------------------------------------------- #

  if (ncol(y) == 0L) deleteIntY <- FALSE
    
  if (deleteIntY) {
    rval <- cbind(x, y[, -1, drop = FALSE])
    attr(rval, "assign") <- c(attr(x, "assign"), attr(y, "assign")[-1])
    attr(rval, "merge") <- rep(c(1L, 2L), c(ncol(x), ncol(y) - 1))
  } else {
    rval <- cbind(x, y)
    attr(rval, "assign") <- c(attr(x, "assign"), attr(y, "assign"))
    attr(rval, "merge") <- rep(c(1L, 2L), c(ncol(x), ncol(y)))
  }
  attr(rval, "contrasts") <- append(attr(x, "contrasts"),
                                    attr(y, "contrasts"))
  attr(rval, "contrasts") <-
    attr(rval, "contrasts")[!duplicated(names(attr(rval, "contrasts")))]

  rownames(rval) <- rownames(x)
  return(rval)
}

olmm_checkMm <- function(x) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Check and modify model matrix.
  ##
  ## Arguments:
  ## x: a model matrix.
  ##
  ## Value:
  ## a model matrix.
  ## ------------------------------------------------------- #
  
  qr.x <- qr(x, LAPACK = FALSE)

  ## automated drops
  if (qr.x$rank < ncol(x)) {
      
    warning("design matrix appears to be column rank-deficient so dropping some coefficients")
    rval <- x[, qr.x$pivot[1:qr.x$rank], drop = FALSE]
    attr(rval, "assign") <- attr(x, "assign")[qr.x$pivot[1:qr.x$rank]]
    attr(rval, "merge") <- attr(x, "merge")[qr.x$pivot[1:qr.x$rank]]
    attr(rval, "contrasts") <- attr(x, "contrasts")
    if ("orig.colnames" %in% names(attributes(x))) {
      attr(rval, "orig.colnames") <- attr(x, "orig.colnames")
    } else {
        attr(rval, "orig.colnames") <- colnames(x)
      }
    rm(x)
    if (qr.x$rank != qr(rval, LAPACK = FALSE)$rank) {
      stop("determination of full column rank design matrix failed")
    }
    
  } else {
      
    rval <- x
    attr(rval, "orig.colnames") <- colnames(rval)

  }
  
  storage.mode(rval) <- "double"
  return(rval)
}

olmm_start <- function(start, dims, parNames, X, W, eta, ranefElMat) {

    ## ------------------------------------------------------- #
    ## Description:
    ## Evaluate user-defined initial values for the parameters.
    ##
    ## Arguments:
    ## start:      named vector with user-defined initial values
    ## dims:       dimension slot of the olmm object
    ## parNames:   list with elements 'fixef' and 'ranefCholFac'
    ##             including the parameter names
    ## X:          model matrix for fixed effects
    ## W:          model matrix for random effects
    ## eta:        matrix for storing the linear predictor
    ## ranefElMat: matrix for transforming the ranefCholFac
    ##             vectors to matrices
    ##
    ## Value:
    ## a list with elements 'fixef', 'ranefCholFac' and
    ## 'coefficients' that are used as initial values.
    ## ------------------------------------------------------- #

    ## checks
    if (!is.null(start)) {
        if (!is.numeric(start) || is.null(names(start)))
            stop("'start' argument should be a named numeric vector")
        if (!all(names(start) %in% unlist(parNames)))
            start <- start[names(start) %in% unlist(parNames)]
      }

    ## set fixed effects

    ## default values
    intDef <- switch(dims["family"],
                     qlogis(ppoints(dims["nEta"])),
                     rep(0.0, dims["nEta"]),
                     rep(0.0, dims["nEta"]))
    fixef <- c(matrix(c(rep(intDef, dims["pInt"]), rep(0.0, dims["nEta"] * (dims["pEtaVar"] - dims["pInt"]))), dims["pEtaVar"], dims["nEta"], byrow = TRUE), rep(0.0, dims["pEtaInv"]))
    names(fixef) <- parNames$fixef

    ## overwrite with 'start'
    subs <- intersect(names(fixef), names(start))
    fixef[subs] <- start[subs]

    ## make a fixed effects matrix
    fixef <- rbind(matrix(as.numeric(fixef[1:(dims["pEtaVar"] * dims["nEta"])]), dims["pEtaVar"] , dims["nEta"], byrow = FALSE), if (dims["pEtaInv"] > 0) matrix(rep(as.numeric(fixef[(dims["pEtaVar"] * dims["nEta"] + 1):dims["p"]]), each = dims["nEta"]), dims["pEtaInv"], dims["nEta"], byrow = TRUE) else NULL)
    rownames(fixef) <- colnames(X)
    colnames(fixef) <-  colnames(eta)

    if (dims["family"] == 3) {
        
        ## transform adjacent category parameters into baseline category
        ## parameters
        T <- 1 * lower.tri(diag(dims["nEta"]), diag = TRUE)
        fixef[1:dims["pEtaVar"],] <- fixef[1:dims["pEtaVar"],] %*% T
    }

    ## set random effect variance components

    ## default values
    ranefCholFac <- c(ranefElMat %*%  c(diag(dims["q"])))
    names(ranefCholFac) <- parNames$ranefCholFac
            
    ## overwrite with 'start'
    subs <- intersect(names(ranefCholFac), names(start))
    ranefCholFac[subs] <- start[subs]

    ## make a matrix
    ranefCholFac <-
        matrix(t(ranefElMat) %*% ranefCholFac, dims["q"], dims["q"])
    ## set row- and column names
    tmp <- c((if (dims["qEtaVar"] > 0) paste("Eta", rep(seq(1, dims["nEta"], 1), each = dims["qEtaVar"]), ":", rep(colnames(W)[attr(W, "merge") == 1], dims["nEta"]), sep = "") else NULL), (if (dims["qEtaInv"] > 0) colnames(W)[attr(W, "merge") == 2] else NULL))
    rownames(ranefCholFac) <- tmp
    colnames(ranefCholFac) <- tmp

    if (dims["family"] == 3 && dims["qEtaVar"] > 0) {

        ## transform adjacent category parameters into baseline category
        ## parameters
        for (i in 1:dims["qEtaVar"]) {
            subs <- seq(i, dims["qEtaVar"] * dims["nEta"], dims["qEtaVar"])
            for (j in 1:(dims["nEta"] - 1)) {
                ranefCholFac[subs[j], ] <-
                    colSums(ranefCholFac[subs[j:dims["nEta"]], ]) 
            }
        }
        ranefCholFac <- t(chol(ranefCholFac %*% t(ranefCholFac)))
    }
    
    ## coefficients
    coefficients <- fixef[1:dims["pEtaVar"],]
    if (dims["pEtaInv"] > 0) coefficients <- c(coefficients,fixef[(dims["pEtaVar"] + 1):dims["pEta"], 1])
    coefficients <- c(coefficients, c(ranefElMat %*% c(ranefCholFac)))
    names(coefficients) <- unlist(parNames)

    return(list(fixef = fixef,
                ranefCholFac = ranefCholFac,
                coefficients = coefficients))
}

olmm_update2 <- function(object, par) {

  ## ------------------------------------------------------- #
  ## Description:
  ## update an olmm object by a new parameter set.
  ##
  ## Arguments:
  ## object: an olmm object. 
  ## par:    a vector that overwrites the coefficients slot
  ##         of the object.
  ## Value:
  ## An updated olmm object.
  ## ------------------------------------------------------- #
  
  if (!(is.vector(par) &&
        is.numeric(par) &&
        length(par) == object@dims["nPar"])) {
    stop("argument 'par' is missspecified")
  }
  .Call("olmm_update_marg", object, par, PACKAGE = "vcolmm")
  .Call("olmm_update_u", object, PACKAGE = "vcolmm")
  return(object)
}

olmm_scoreVar <- function(object, Nmax = 0L) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Computes variance of scores.
  ##
  ## Arguments:
  ## object: a fitted 'olmm' object
  ## Nmax:   the number of observation per subject
  ##         for which the matrix should be computed.
  ##         Typically the maximal number of observations
  ##         observed. 0 means that the unbalanced case
  ##         is ignored
  ## ------------------------------------------------------- #
  
  Ni <- table(object@subject)
  if (Nmax > 0) {
    subsObs <- object@subject %in% levels(object@subject)[Ni == Nmax]
  } else {
    subsObs <- rep(TRUE, object@dims["n"])
  }
  return(crossprod(-object@score_obs[subsObs,,drop=FALSE]) / sum(subsObs))
}

olmm_scoreCovWin <- function(object, Nmax = 0L) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Computes within-subject covariance of scores.
  ##
  ## Arguments:
  ## object: a fitted 'olmm' object
  ## Nmax:   the number of observation per subject
  ##         for which the matrix should be computed.
  ##         Typically the maximal number of observations
  ##         observed. 0 means that the unbalanced case
  ##         is ignored
  ## ------------------------------------------------------- #
  
  Ni <- table(object@subject)
  subsSbj <- if (Nmax > 0) Ni == Nmax else rep(TRUE, object@dims["N"])
  n <- sum(object@subject %in% levels(object@subject)[subsSbj])
  return((crossprod(-object@score_sbj[subsSbj,,drop=FALSE]) -
          n * olmm_scoreVar(object, Nmax)) / sum(Ni[subsSbj] * (Ni[subsSbj] - 1)))
}

olmm_scoreCovBet <- function(object) {
  
  ## ------------------------------------------------------- #
  ## Description:
  ## Computes between-subject covariance of scores.
  ##
  ## Arguments:
  ## object: a fitted 'olmm' object
  ## ------------------------------------------------------- #
  
  n <- object@dims["n"]
  Ni <- table(object@subject)
  return(-crossprod(object@score_sbj) / (n^2 - n - sum(Ni * (Ni-1))))
}

olmm_f_scoreTransfMat <- function(T, Tindex, sVar, sCovWin, sCovBet, Nmax) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Difference between adjusted within- and adjusted between-
  ## subject covariance.
  ## ------------------------------------------------------- #

  adjScoreCovWin <-
    sVar %*% t(T) + T %*% sVar + (Nmax - 2) * T %*% sVar %*% t(T) +
    sCovWin + (Nmax - 2) * sCovWin %*% t(T) + (Nmax - 2) * T %*% t(sCovWin) +
      ((Nmax - 1)^2 - (Nmax - 2)) * T %*% sCovWin %*% t(T)

  adjScoreCovBet <-
    sCovBet +
        (Nmax - 1) * sCovBet %*% t(T) + (Nmax - 1) * T %*% t(sCovBet) +
        (Nmax - 1)^2 * T %*% sCovBet %*% t(T)
  
  rval <- adjScoreCovWin - adjScoreCovBet
  subs <- which(!duplicated(c(Tindex)) & c(Tindex) != 0)
  return(c(rval[subs]))
}

olmm_g_scoreTransfMat <- function(T, Tindex, sVar, sCovWin, sCovBet, Nmax) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Gradient for the olmm_f_scoreTransfMat() function above
  ## ------------------------------------------------------- #

  k <- ncol(T)
  nPar <- max(Tindex)
  subs <- which(!duplicated(c(Tindex)) & c(Tindex) != 0)
  rval <- matrix(0, nPar, nPar)
  ## each iteration evaluates the gradient for one element in T
  ## regarding all elements in the objective function
  for (i in 1:nPar) {
    ind <- 1L * (Tindex == i)
    val <- T[which(ind == 1L)][1]

    gAdjScoreCovWin <-
      sVar %*% t(ind) + ind %*% sVar + 2 * val * (Nmax - 2) * ind %*% sVar %*% t(ind) +
        (Nmax - 2) * sCovWin %*% t(ind) + (Nmax - 2) * ind %*% sCovWin +
          2 * val * ((Nmax - 1)^2 - (Nmax - 2)) * ind %*% sCovWin %*% t(ind)

    gAdjScoreCovBet <-
      (Nmax - 1) * sCovBet %*% t(ind)  + (Nmax - 1) * ind %*% sCovBet +
           2 * val * (Nmax - 1)^2 * ind %*% sCovBet %*% t(ind)
    
    rval[, i] <- (gAdjScoreCovWin + gAdjScoreCovBet)[subs]
  }
  return(rval)
}

olmm_scoreTransfMat <- function(object, method = c("symmetric", "unconstraint"),
                                Nmax = NULL, terms = NULL, control = list(),
                                start = NULL, verbose = FALSE, omit.terms = TRUE,
                                silent = FALSE) {

  method <- match.arg(method)
  control <- appendDefArgs(control,
                           list(reltol = 1e-6,
                                maxreltol = 1e-3,
                                maxit = 100L,
                                stopreltol = 1e100))
  
  ## get required characteristics of the model
  n <- object@dims["n"]
  Ni <- table(object@subject)
  if (is.null(Nmax)) Nmax <- max(Ni)
  if (!any(Ni == Nmax))
    stop("at least one subject must have ",
         Nmax, " observations")
  if (!silent && length(table(Ni)) > 1)
        warning("computation is based on scores of the ", sum(Ni == Nmax),
                " of ", nlevels(object@subject), " subjects with ", Nmax,
                " observations only. ")
  
  sVar <- olmm_scoreVar(object, Nmax = Nmax)
  sCovWin <- olmm_scoreCovWin(object, Nmax = Nmax)
  sCovBet <- olmm_scoreCovBet(object)
  
  ## reduce to coefficient subset if intended
  if (!is.null(terms)) {
    sVar <- sVar[terms, terms, drop = FALSE]
    sCovWin <- sCovWin[terms, terms, drop = FALSE]
    sCovBet <- sCovBet[terms, terms, drop = FALSE]
  } else {
    terms <- 1:ncol(sVar)
  } 
  k <- length(terms)
  
  ## set initial values
  T <- if (is.null(start)) {
    matrix(0, k, k, dimnames = list(rownames(sVar), colnames(sVar)))
  } else {
    start
  }
  Tindex <- matrix(0, k, k, dimnames = list(rownames(sVar), colnames(sVar)))
  subs <- if (method == "symmetric") lower.tri(T, TRUE) else matrix(TRUE, k, k)
  ## omit off-diagonal terms which are zero
  if (omit.terms) {
    subs[abs(sVar) + abs(sCovWin) + abs(sCovBet) < .Machine$double.eps] <- FALSE
  }  
  nPar <- sum(subs)
  Tindex[subs] <- 1:nPar
  if (method == "symmetric")
    Tindex[upper.tri(Tindex, FALSE)] <- t(Tindex)[upper.tri(Tindex, FALSE)]
  par <- rep(0, nPar)
  fEval <- rep(Inf, nPar)
  nit <- 0; error <- FALSE; eps <- 2 * control$reltol;
  
  ## optimize by Newton's algorithm
  while (!error && nit < control$maxit && eps >= control$reltol) {
    nit <- nit + 1
    fEval <- olmm_f_scoreTransfMat(T, Tindex, sVar, sCovWin, sCovBet, Nmax)
    gEval <- olmm_g_scoreTransfMat(T, Tindex, sVar, sCovWin, sCovBet, Nmax)
    par <- try(solve(gEval, -fEval) + par, silent = TRUE)
    eps <- max(abs(fEval / sCovBet[subs]))
    if (class(par) != "try-error") {
      T[] <- c(0, par)[Tindex + 1]
      if (verbose & nit > 1)
        cat("\nnit =", nit,
            "max|f| =", format(max(abs(fEval)), digits = 3, scientific = TRUE),
            "max|diff/f| =", format(eps, digits = 3, scientific = TRUE))
    }
    if (class(par) == "try-error" || eps > control$stopreltol || is.nan(eps))
      error <- TRUE
  } 
  
  ## check convergence
  if (nit >= control$maxit | error) {
    mess <- paste("optimization not converged: nit =", nit,
                  "reltol =", format(eps, digits = 3, scientific = TRUE), "\n")
    if (verbose) { cat("\n"); cat(mess); }
    if (!silent) warning(mess)
  } else {
    if (verbose)
      cat("\noptimization converged: nit =", nit,
          "max|diff/f| =", format(eps, digits = 3, scientific = TRUE), "\n")
  }
  
  attr(T, "conv") <- as.integer((nit != control$maxit) & !error)
  attr(T, "nit") <- nit
  attr(T, "eps") <- eps
  rownames(T) <- colnames(T) <- colnames(sVar)
  
  ## return transformation matrix
  return(T)
} 


olmm_simPredictors <- function(object, nobs) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Simulate predictors from a fitted 'olmm' object.
  ##
  ## Arguments:
  ## object:    an 'olmm' object
  ## nobs:      the number of observations to be simulated
  ##
  ## Value:
  ## A list with the frame and model matrices of the simulated
  ## predictors
  ## ------------------------------------------------------- #

  if (nobs > 0L) {
    frame <- model.frame(object)
    X <- model.matrix(object, "fixef")
    W <- model.matrix(object, "ranef")
    index <- sample(1:nrow(frame), nobs, replace = TRUE)
    frame <- frame[index,,drop = FALSE]
    frame <- frame[, -c(1, which(colnames(frame) == object@subjectName)), drop = FALSE]
    X <- X[index,,drop = FALSE]
    W <- W[index,,drop = FALSE]
    rownames(frame) <- rownames(X) <- rownames(W) <- 1:nrow(frame)
    rval <- list(frame = frame, X = X, W = W)
  } else {
    rval <- list()
  }
  return(rval)
}

olmm_addScores <- function(object) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Simulates data so that the data become balanced.
  ##
  ## Arguments:
  ## object:    an 'olmm' object
  ##
  ## Value:
  ## A list with slots 'score_obs' and 'subject'.
  ## ------------------------------------------------------- #
  
  ## extract data to be completed
  frame <- model.frame(object)
  subjectName <- object@subjectName

  ## collect information and set subset
  subjectLevs <- levels(object@subject)
  Ni <- table(object@subject)
  Nmax <- max(Ni)
  subset <- object@subject %in% levels(object@subject)[Ni < Nmax]
  frame <- frame[subset,, drop = FALSE]
  frame[, subjectName] <- droplevels(frame[, subjectName])

  ## get the number of observations per subject to be inputed
  Ninpute <- Nmax - table(frame[, subjectName])
  if (max(Ninpute) == 0L) return(list(score_obs= NULL, subject = NULL))
  
  ## simulate predictors (arbitrary, we just need values)
  predictors <- olmm_simPredictors(object, sum(Ninpute))
  newdata <- predictors$frame
  newdata[, subjectName] <- factor(rep(levels(frame[, subjectName]), Ninpute),
                                   levels = levels(frame[, subjectName]))

  ## simulate responses
  ranef <- ranef(object)[levels(newdata[, subjectName]),, drop = FALSE]
  ranef <- ranef[sort(unique(newdata[, subjectName])),, drop = FALSE]
  newdata <- cbind(newdata,
                   simulate.olmm(object, newdata = newdata, ranef = ranef))
  newdata <- newdata[, colnames(model.frame(object)), drop = FALSE]
  
  ## build a new 'olmm' object with additional observations
  newObj <- object
  newObj@frame <- rbind(frame, newdata)
  newObj@y <- newObj@frame[,1]
  newObj@subject <- newObj@frame[, subjectName]
  newObj@X <- rbind(object@X[subset,,drop=FALSE], predictors$X)
  newObj@W <- rbind(object@W[subset,,drop=FALSE], predictors$W)
  newObj@weights <- object@weights_sbj[as.integer(newObj@subject)]
  newObj@offset <- c(object@offset[subset], rep(0, nrow(newdata)))
  newObj@dims["n"] <- nrow(newObj@frame)
  newObj@dims["N"] <- nlevels(newObj@subject)
  newObj@eta <- rbind(object@eta[subset,,drop=FALSE],
                    matrix(0, nrow(newdata), ncol(object@eta)))
  newObj@score_obs <-
    rbind(object@score_obs[subset,,drop=FALSE],
          matrix(0, nrow(newdata), ncol(object@score_obs)))
  newObj@score_sbj <- newObj@score_sbj[levels(newObj@subject),,drop=FALSE]
  newObj@u <- newObj@u[levels(newObj@subject),,drop=FALSE]
  .Call("olmm_update_marg", newObj, newObj@coefficients, PACKAGE = "vcolmm")
  subs <- (nrow(frame) + 1):newObj@dims["n"]
  score_obs <- newObj@score_obs[subs,,drop=FALSE]
  subject <- factor(newObj@subject[subs], levels = subjectLevs)
  return(list(score_obs = score_obs, subject = subject))
}
