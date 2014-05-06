## --------------------------------------------------------- #
## Author:          Reto Buergin, rbuergin@gmx.ch
## Date:            2014-05-03
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
## olmm_check_mm:        check and modify model matrix
## olmm_start:           set initial values
## olmm_mergeMm:         merge the predictor-variable and predictor-invariant
##                       model matrices
##                     
## gcauchy:              derivate of dcauchy
## glogis:               derivate of dlogis
## gnorm:                derivate of dnorm
##
## Functions for estfun.olmm:
## olmm_scoreVar:        computes the variance of the observation scores.
## olmm_scoreCovWin:     computes the intra-subject covariance of
##                       the observation scores.
## olmm_scoreCovBet:     computes the between-subject covariance
##                       of the observation scores.
## olmm_f_decormat:      the equation to be optimized to zero:
##                       the difference between the adjusted
##                       intra-subject covariance and the adjusted
##                       between-subject covariance of likelihood scores.
## olmm_g_decormat:      the derivation of olmm_fDecorrelate.
## decormat.olmm:        computes the transformation matrix for
##                       removing the intra-subject
##                       correlation of ML scores.
## olmm_simPredictors:   extracts randomly observed predictor
##                       combinations
## olmm_addObsScores:    computes scores to complete an unbalance
##                       data set 
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
  parNew <- slot(object, "coefficients")
  parNew[!restricted] <- par[!restricted]
  .Call("olmm_update_marg", object, as.numeric(parNew), PACKAGE = "vcrpart")
  return(-slot(object, "logLik"))
}

olmm_gr <- function(par, restricted) {
  if (!exists("object")) object <- new(Class = "olmm")
  scoreNew <- slot(object, "score")
  scoreNew[restricted] <- c(0.0)
  return(-scoreNew)
}

olmm_optim_setup <- function(x, env = parent.frame()) {  

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

  numGrad <- x$numGrad
  x <- x[!names(x) %in% c("start", "restricted")]
  rval <- list(par = NULL, # first 4 arguments must be exactly in that order!
               fn = olmm_fn,
               gr = if (numGrad) NULL else olmm_gr,
               restricted = NULL,
               fit = x$fit)
  rval <- append(rval, x[intersect(names(formals(rval$fit)), names(x))])
  if (rval$fit == "nlminb") names(rval)[1:3] <- names(formals(nlminb)[1:3])
  
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
  abbCe <- names(fixef(object, "predictor-variable"))
  abbCe <- gsub("Eta", "", abbCe, fixed = TRUE)
  abbCe <- gsub("(Intercept)", "I", abbCe, fixed = TRUE)
  abbCe <- abbreviate(abbCe)

  ## abbreviations for predictor-invariant fixed effects
  abbGe <- names(fixef(object, "predictor-invariant"))
  abbGe <- abbreviate(abbGe)

  ## abbreviations for Choleski factors
  q <- slot(object, "dims")["q"]
  abbReCF <- paste("reCF", 1:(q * (q + 1L) / 2L), sep = "")
  
  return(c(abbCe, abbGe, abbReCF))
}

olmm_merge_mm <- function(x, y, deleteIntY = TRUE) {

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
  if (deleteIntY & "(Intercept)" %in% colnames(y)) {
    rval <- cbind(x, y[, -1L, drop = FALSE])
    attr(rval, "assign") <- c(attr(x, "assign"), attr(y, "assign")[-1L])
    attr(rval, "merge") <- rep(c(1L, 2L), c(ncol(x), ncol(y) - 1L))
  } else {
    rval <- cbind(x, y)
    attr(rval, "assign") <- c(attr(x, "assign"), attr(y, "assign"))
    attr(rval, "merge") <- rep(c(1L, 2L), c(ncol(x), ncol(y)))
  }
  attr(rval, "contrasts") <-
    append(attr(x, "contrasts"), attr(y, "contrasts"))
  attr(rval, "contrasts") <-
    attr(rval, "contrasts")[!duplicated(names(attr(rval, "contrasts")))]
  rownames(rval) <- rownames(x)
  return(rval)
}

olmm_check_mm <- function(x) {

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
  rank <- qr.x$rank
  
  ## automated drops
  if (rank < ncol(x)) {
      
    warning("design matrix appears to be column rank-deficient ",
            "so dropping some coefs")
    
    dropterms <- function(x, keep) {
      rval <- x[, keep, drop = FALSE]
      attr(rval, "assign") <- attr(x, "assign")[keep]
      attr(rval, "merge") <- attr(x, "merge")[keep]
      attr(rval, "contrasts") <- attr(x, "contrasts")
      if ("orig.colnames" %in% names(attributes(x))) {
        attr(rval, "orig.colnames") <- attr(x, "orig.colnames")
      } else {
        attr(rval, "orig.colnames") <- colnames(x)
      }
      return(rval)
    }
   
    subs <- qr.x$pivot[1:qr.x$rank]
    x <- dropterms(x, subs)
    
    if (rank != qr(x, LAPACK = FALSE)$rank)
      stop("determination of full column rank design matrix failed")
  } else {

    attr(x, "orig.colnames") <- colnames(x)
  } 
  storage.mode(x) <- "double"
  return(x)
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
    fixef <- c(matrix(c(rep(intDef, dims["pInt"]), rep(0.0, dims["nEta"] * (dims["pCe"] - dims["pInt"]))), dims["pCe"], dims["nEta"], byrow = TRUE), rep(0.0, dims["pGe"]))
    names(fixef) <- parNames$fixef

    ## overwrite with 'start'
    subs <- intersect(names(fixef), names(start))
    fixef[subs] <- start[subs]

    ## make a fixed effects matrix
    fixef <- rbind(matrix(as.numeric(fixef[1:(dims["pCe"] * dims["nEta"])]), dims["pCe"] , dims["nEta"], byrow = FALSE), if (dims["pGe"] > 0) matrix(rep(as.numeric(fixef[(dims["pCe"] * dims["nEta"] + 1):dims["p"]]), each = dims["nEta"]), dims["pGe"], dims["nEta"], byrow = TRUE) else NULL)
    rownames(fixef) <- colnames(X)
    colnames(fixef) <-  colnames(eta)

    if (dims["family"] == 3L && dims["pCe"] > 0L) {
        
        ## transform adjacent category parameters into baseline category
        ## parameters
        T <- 1 * lower.tri(diag(dims["nEta"]), diag = TRUE)
        fixef[1:dims["pCe"],] <- fixef[1:dims["pCe"],] %*% T
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
    tmp <- c((if (dims["qCe"] > 0) paste("Eta", rep(seq(1, dims["nEta"], 1), each = dims["qCe"]), ":", rep(colnames(W)[attr(W, "merge") == 1], dims["nEta"]), sep = "") else NULL), (if (dims["qGe"] > 0) colnames(W)[attr(W, "merge") == 2] else NULL))
    rownames(ranefCholFac) <- tmp
    colnames(ranefCholFac) <- tmp

    if (dims["family"] == 3 && dims["qCe"] > 0) {

        ## transform adjacent category parameters into baseline category
        ## parameters
        for (i in 1:dims["qCe"]) {
            subs <- seq(i, dims["qCe"] * dims["nEta"], dims["qCe"])
            for (j in 1:(dims["nEta"] - 1)) {
                ranefCholFac[subs[j], ] <-
                    colSums(ranefCholFac[subs[j:dims["nEta"]], ]) 
            }
        }
        ranefCholFac <- t(chol(ranefCholFac %*% t(ranefCholFac)))
    }
    
    ## coefficients
    coefficients <- numeric()
    if (dims["pCe"] > 0L) coefficients <- fixef[1:dims["pCe"],]
    if (dims["pGe"] > 0) coefficients <- c(coefficients,fixef[(dims["pCe"] + 1):dims["pEta"], 1])
    coefficients <- c(coefficients, c(ranefElMat %*% c(ranefCholFac)))
    names(coefficients) <- unlist(parNames)

    return(list(fixef = fixef,
                ranefCholFac = ranefCholFac,
                coefficients = coefficients))
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
  
  Ni <- table(slot(object, "subject"))
  if (Nmax > 0) {
    subsObs <- slot(object, "subject") %in% levels(slot(object, "subject"))[Ni == Nmax]
  } else {
    subsObs <- rep(TRUE, slot(object, "dims")["n"])
  }
  return(crossprod(-slot(object, "score_obs")[subsObs,,drop=FALSE]) / sum(subsObs))
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
  
  Ni <- table(slot(object, "subject"))
  subsSbj <- if (Nmax > 0) Ni == Nmax else rep(TRUE, slot(object, "dims")["N"])
  n <- sum(slot(object, "subject") %in% levels(slot(object, "subject"))[subsSbj])
  return((crossprod(-slot(object, "score_sbj")[subsSbj,,drop=FALSE]) -
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
  
  n <- slot(object, "dims")["n"]
  Ni <- table(slot(object, "subject"))
  return(-crossprod(slot(object, "score_sbj")) / (n^2 - n - sum(Ni * (Ni-1))))
}

olmm_f_decormat <- function(T, Tindex, sVar, sCovWin, sCovBet, Nmax) {

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

olmm_g_decormat <- function(T, Tindex, sVar, sCovWin, sCovBet, Nmax) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Gradient for the olmm_f_decormat() function above
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


olmm_simPredictors <- function(object, nobs) {

  ## ------------------------------------------------------- #
  ## Description:
  ## Simulate predictors based on the data of a fitted 'olmm'
  ## object.
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
    X <- model.matrix(object, "fe")
    W <- model.matrix(object, "re")
    index <- sample(1:nrow(frame), nobs, replace = TRUE)
    frame <- frame[index,,drop = FALSE]
    frame <- frame[, -c(1, which(colnames(frame) == slot(object, "subjectName"))), drop = FALSE]
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
  subjectName <- slot(object, "subjectName")

  ## collect information and set subset
  subjectLevs <- levels(slot(object, "subject"))
  Ni <- table(slot(object, "subject"))
  Nmax <- max(Ni)
  subset <- slot(object, "subject") %in% levels(slot(object, "subject"))[Ni < Nmax]
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
  slot(newObj, "frame") <- rbind(frame, newdata)
  slot(newObj, "y") <- slot(newObj, "frame")[,1]
  slot(newObj, "subject") <- slot(newObj, "frame")[, subjectName]
  slot(newObj, "X") <- rbind(slot(object, "X")[subset,,drop=FALSE], predictors$X)
  slot(newObj, "W") <- rbind(slot(object, "W")[subset,,drop=FALSE], predictors$W)
  slot(newObj, "weights") <- slot(object, "weights_sbj")[as.integer(slot(newObj, "subject"))]
  slot(newObj, "offset") <- rbind(slot(object, "offset")[subset,,drop=FALSE],
                         matrix(0.0, nrow(newdata), slot(object, "dims")["nEta"]))
  slot(newObj, "dims")["n"] <- nrow(slot(newObj, "frame"))
  slot(newObj, "dims")["N"] <- nlevels(slot(newObj, "subject"))
  slot(newObj, "eta") <- rbind(slot(object, "eta")[subset,,drop=FALSE],
                                    matrix(0, nrow(newdata), ncol(slot(object, "eta"))))
  slot(newObj, "score_obs") <-
    rbind(slot(object, "score_obs")[subset,,drop=FALSE],
          matrix(0, nrow(newdata), ncol(slot(object, "score_obs"))))
  slot(newObj, "score_sbj") <- slot(newObj, "score_sbj")[levels(slot(newObj, "subject")),,drop=FALSE]
  slot(newObj, "u") <- slot(newObj, "u")[levels(slot(newObj, "subject")),,drop=FALSE]
  .Call("olmm_update_marg", newObj, slot(newObj, "coefficients"), PACKAGE = "vcrpart")
  subs <- (nrow(frame) + 1):slot(newObj, "dims")["n"]
  score_obs <- slot(newObj, "score_obs")[subs,,drop=FALSE]
  subject <- factor(slot(newObj, "subject")[subs], levels = subjectLevs)
  return(list(score_obs = score_obs, subject = subject))
}
