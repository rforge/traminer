## --------------------------------------------------------- #
## Author:          Reto Buergin, rbuergin@gmx.ch
## Date:            2013-07-12
##
## Description:
## Utility functions for the olmm function (see olmm.R). Some
## are experimental and not listed in the namespace.
##
## References and dependencies:
## statmod:         http://cran.r-project.org/web/packages/statmod/index.html
## ucminf:          http://cran.r-project.org/web/packages/ucminf/index.html
##
## Contents:
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
  parNew <- object@coefficients
  parNew[!restricted] <- par[!restricted]
  .Call("olmm_update_marg", object, as.numeric(parNew), PACKAGE = "vcolmm")
  return(-object@logLik)
}

olmm_gr <- function(par, restricted) {
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
  yName <- rownames(attr(terms(ff), "factors"))[1L]
  
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

