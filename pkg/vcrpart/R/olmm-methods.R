## --------------------------------------------------------- #
## Author:          Reto Buergin, rbuergin@gmx.ch
## Date:            2014-05-10
##
## Description:
## methods for olmm objects.
##
## anova:       Likelihood-ratio tests for the comparison of
##              models
## coef, coefficients: Extract model coefficients
## deviance:    -2*Log-likelihood at the ML estimator
## drop1:       Drop single fixed effects
## estfun:      Negative scores
## extractAIC:  Extract the AIC
## fitted:      Extract fitted values from the model
## fixef:       Extract fixed effect parameters
## formula:     Extracts 'formula'
## gefp:        Extract cumulated decorrelated score process
## getCall:     Extracts 'call'
## logLik:      Log-likelihood at the ML estimator
## model.frame: Model frame (all needed variables)
## model.matrix: Model matrix (for the fixed effects)
## predict:     Predict from the fitted model
## print:       Print summary output (method for olmm and
##              olmm.summary objects)
## ranef:       Extract predicted random effects
## ranefCov:    Covariance-matrix of random effect terms
## resid, residuals Extract different types of residuals from the
##              fitted model
## show:        Print summary output (method for olmm and
##              olmm.summary objects)
## simulate:    Simulate responses based on a fitted model
## summary:     Extract summary information
## terms, terms.olmm: Extracting the terms of the model frame
##              for fixed effects
## update:      Refits a model
## VarCorr, print.VarCorr.olmm: Extract variance and standard
##              deviation of random effects and their correlation
## vcov:        Variance-covariance matrix for fixed effect parameters
## weights:     Weights
##
## Modifications:
## 2013-03-17: changed many methods to S3 methods (as in lme4)
## 2013-09-06: modify formula() method. Now the formula slot
##             is called
## 2013-09-06: add S3 for terms() method
## 2013-09-06: add drop1() method
##
## To do:
## - improve update method
## - plot methods
## - estfun.olmm: handle equal zero random effects
## --------------------------------------------------------- #

anova.olmm <- function(object, ...) {

  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- if (length(dots) > 0) {
    sapply(dots, is, "olmm")
  } else {
    logical(0)
  }

  if (any(modp)) {

    ## multiple models
    opts <- dots[!modp]
    mods <- c(list(object), dots[modp])
    names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
                          as.character)
    mods <- mods[order(sapply(lapply(mods, logLik),
                              attr, "df"))]
    calls <- lapply(mods, getCall)
    data <- lapply(calls, "[[", "data")
    if (any(data != data[[1]]))
      stop("all models must be fit to the same data object")
    header <- paste("Data:", data[[1]])
    subset <- lapply(calls, "[[", "subset")
    if (any(subset != subset[[1]]))
      stop("all models must use the same subset")
    if (!is.null(subset[[1]]))
      header <-
        c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
    llks <- lapply(mods, logLik)
    Df <- sapply(llks, attr, "df")
    llk <- unlist(llks)
    chisq <- 2 * pmax(0, c(NA, diff(llk)))
    dfChisq <- c(NA, diff(Df))
    val <- data.frame(Df = Df,
                      AIC = sapply(llks, AIC),
                      BIC = sapply(llks, BIC),
                      logLik = llk,
                      "Chisq" = chisq,
                      "Chi Df" = dfChisq,
                      "Pr(>Chisq)" = pchisq(chisq,
                        dfChisq,
                        lower.tail = FALSE),
                      row.names = names(mods), check.names = FALSE)
    class(val) <- c("anova", class(val))
    attr(val, "heading") <-
      c(header, "Models:",

        paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
                                 "[[", "formula"), deparse), length))),
              unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
              sep = ": "))
    return(val)

  } else {
    ## single model
    stop("single argument anova for 'olmm' objects not yet implemented")
  }
}

coef.olmm <- function(object, which = c("all", "fe"), ...) {

  which <- match.arg(which)
  if (which == "fe") return(fixef(object))
  
  dims <- slot(object, "dims")
  if (slot(object, "family")$family == "adjacent") {
    T <- diag(dims["nPar"])
    subsRows <- seq(1, dims["pCe"] * (dims["nEta"] - 1L), 1L)
    subsCols <- seq(dims["pCe"] + 1L,
                    dims["pCe"] * dims["nEta"], 1L)
    if (length(subsRows) == 1L) {
      T[subsRows, subsCols] <- c(-1.0)
    } else {
      diag(T[subsRows, subsCols]) <- c(-1.0)
    }
    rval <- c(T %*% slot(object, "coefficients"))
    names(rval) <- names(slot(object, "coefficients"))
  } else {
    rval <- slot(object, "coefficients")
  }
  if (dims["hasRanef"] == 0L)
    rval <- rval[!grepl("ranefCholFac", names(rval))]
  
  return(rval)
}

coefficients.olmm <- coef.olmm

decormat.olmm <- function(object, method = c("symmetric", "unconstraint"),
                          Nbal = NULL, parm = NULL, control = list(),
                          verbose = FALSE, drop = TRUE, silent = FALSE) {

  method <- match.arg(method)

  ## set default control parameters
  control_def <- list(reltol = 1e-6, maxreltol = 1e-3,
                      maxit = 100L, stopreltol = 1e100)
  control <- appendDefArgs(control, control_def)
  
  ## get required characteristics of the model
  n <- nobs(object)
  Ni <- table(slot(object, "subject"))
  if (is.null(Nbal)) Nbal <- max(Ni)
  if (!any(Ni == Nbal))
    stop("at least one subject should have ", Nbal, " observations")
  if (verbose)
    cat("\nT is based on scores of the", sum(Ni == Nbal),
        "of", slot(object, "dims")["N"], "subjects with Ni >=", Nbal, " obs.")
  
  sVar <- olmm_scoreVar(object, Nmax = Nbal)
  sCovWin <- olmm_scoreCovWin(object, Nmax = Nbal)
  sCovBet <- olmm_scoreCovBet(object)
  
  ## reduce to coefficient subset if intended
  if (!is.null(parm)) {
    sVar <- sVar[parm, parm, drop = FALSE]
    sCovWin <- sCovWin[parm, parm, drop = FALSE]
    sCovBet <- sCovBet[parm, parm, drop = FALSE]
  } else {
    parm <- 1:ncol(sVar)
  } 
  k <- length(parm)
  
  ## set initial values (currently omitted)
  start <- NULL
  T <- if (is.null(start)) {
    matrix(0, k, k, dimnames = list(rownames(sVar), colnames(sVar)))
  } else {
    start
  }
  Tindex <- matrix(0, k, k, dimnames = list(rownames(sVar), colnames(sVar)))
  subs <- if (method == "symmetric") lower.tri(T, TRUE) else matrix(TRUE, k, k)
  ## omit off-diagonal terms which are zero
  if (drop) {
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
    fEval <- olmm_f_decormat(T, Tindex, sVar, sCovWin, sCovBet, Nbal)
    gEval <- olmm_g_decormat(T, Tindex, sVar, sCovWin, sCovBet, Nbal)
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

deviance.olmm <- function(object, ...) {
  return(-2 * slot(object, "logLik"))
}

estfun.olmm <- function(x, level = c("observation", "subject"),
                        predecor = FALSE, complete = predecor,
                        Nbal = NULL, subset = NULL,
                        nuisance = NULL, control = list(),
                        verbose = FALSE, silent = FALSE, ...) {

  level <- match.arg(level)
  if (verbose) cat("* extract scores from fitted object ... ")
  xold <- x
  
  ## set subset
  if (is.null(subset)) {
    subset <- rep(TRUE, slot(x, "dims")["n"])
  } else {
    if (is.character(subset)) subset <- rownames(model.frame(x)) %in% subset
    if (is.numeric(subset) && length(subset) == slot(x, "dims")["n"] &&
        all(subset %in% c(1, 0))) subset <- as.logical(subset)
    if (is.numeric(subset)) subset <- (1:slot(x, "dims")["n"]) %in% subset
  }
  
  ## check for balanced data and automatically modify subset
  if (complete) {

    Ni <- table(slot(x, "subject")[subset])

    ## set Nbal
    if (!is.null(Nbal)) {
      if (Nbal > max(Ni))
        stop("at least one subject must have equal or more than ",
             Nbal, " observations.")
    } else {
      Nbal <- sort(unique(Ni))[rev(which(table(Ni) > 50))[1]]
      if (is.na(Nbal)) Nbal <- max(Ni)
    }  
    if (any(Ni > Nbal)) {
      FUN <- function(x) x[1:min(Nbal, length(x))]
      subs <- unlist(tapply(1:length(slot(x, "subject")), slot(x, "subject"), FUN))
      subs <- (1:slot(x, "dims")["n"]) %in% subs
      if (verbose)
        cat("\n\tomit ", sum(!subs), " obs. from subjects with Ni > ", Nbal)
      subset <- subset & subs
    }

    ## recompute the scores for the subset
    if (any(!subset)) {
      slot(x, "frame") <- model.frame(x)[subset,,drop = FALSE]
      slot(x, "y") <- slot(x, "y")[subset]
      slot(x, "subject") <- slot(x, "subject")[subset]
      attrX <- attributes(slot(x, "X"))
      slot(x, "X") <- slot(x, "X")[subset,,drop=FALSE]
      attributes(slot(x, "X")) <- appendDefArgs(attributes(slot(x, "X")), attrX)
      attrW <- attributes(slot(x, "W"))
      slot(x, "W") <- slot(x, "W")[subset,,drop=FALSE]
      attributes(slot(x, "W")) <- appendDefArgs(attributes(slot(x, "W")), attrW)
      slot(x, "weights") <- slot(x, "weights_sbj")[as.integer(slot(x, "subject"))]
      slot(x, "offset") <- slot(x, "offset")[subset,,drop=FALSE]
      slot(x, "dims")["n"] <- nrow(slot(x, "frame"))
      slot(x, "eta") <- slot(x, "eta")[subset,,drop=FALSE]
      slot(x, "score_obs") <- slot(x, "score_obs")[subset,,drop=FALSE]
      .Call("olmm_update_marg", x, slot(x, "coefficients"), PACKAGE = "vcrpart")
    }
  }

  if (verbose && !any(Ni > Nbal)) cat("OK")
  
  scores <-  -slot(x, "score_obs")
  subject <- slot(x, "subject")
  Ni <- table(subject)
  
  ## checks
  if (slot(x, "dims")["hasRanef"] == 0L) {
    predecor <- FALSE
    complete <- FALSE
  }
  
  ## get terms
  parm <- 1:slot(x, "dims")["nPar"] # internal variable
  if (!is.null(nuisance) & is.character(nuisance))
    nuisance <- which(names(coef(x)) %in% nuisance)
  parm <- setdiff(parm, nuisance)
  scores <- scores[, parm, drop = FALSE]
    
  if (level == "observation") {

    n <- nrow(scores)
    k <- ncol(scores)
    
    if (predecor) {
        
      ## compute transformation matrix
      if (verbose) cat("\n* compute transformation matrix ...")
      T <- decormat.olmm(object = x, parm = parm, verbose = verbose,
                         control = control, silent = silent, ...)
      
      ## if transformation failed, return raw scores
      if (attr(T, "conv") > 0L) {
        
        if (verbose) cat("\n* transforming scores ... ")
        
        ## transformation matrix for one subject
        Ti <- kronecker(matrix(1, Nbal, Nbal) - diag(Nbal), T)
        diag(Ti) <- 1
      
        if (all(Ni == Nbal) | !complete) {

          if (!silent && any(Ni < Nbal))
            warning("the transformation method works only with balanced data.")

          sbj <- factor(c(as.character(subject), rep(levels(subject), Nbal - Ni)),
                        levels = levels(subject))
          subsAdd <- c(sapply(Ni, function(n) 1:Nbal > n))
          sT <- rbind(scores, matrix(0, length(sbj) - length(subject), k))
          subsOrd <- order(sbj)
          sTmp <- matrix(c(t(sT[subsOrd,,drop=FALSE])), Nbal * k, nlevels(subject))
          sTmp <- matrix(c(Ti %*% sTmp), nrow(sT), k, byrow = TRUE)
          sTmp <- sTmp[!subsAdd,,drop=FALSE]        
          scores[rownames(scores)[order(subject)],] <- sTmp
          if (verbose) cat("OK")
          
        } else {

          sT <- 0.0 * scores

          ## transform complete cases
          subsSbj <- which(Ni == Nbal)
          subsObs <- which(subject %in% levels(subject)[subsSbj])
          subsOrd <- subsObs[order(subject[subsObs])]
          sTmp <- matrix(c(t(scores[subsOrd,,drop=FALSE])), Nbal * k, length(subsSbj))
          sTmp <- matrix(c(Ti %*% sTmp), length(subsObs), k, byrow = TRUE)
          sT[rownames(sT)[subsOrd],] <- sTmp
          
          ## repeatedly transform incomplete cases
          if (verbose) cat("\n* data are unbalanced.",
                           "Apply imputation algorithm with Nbal =", Nbal, "...")
          subsSbj <- which(Ni < Nbal)
          subsObs <- which(subject %in% levels(subject)[subsSbj])
          subsOrd <- subsObs[order(subject[subsObs])]
          sTmp <- 0.0 * scores[subsObs, ]
          subsAdd <- c(sapply(Ni[subsSbj], function(n) 1:Nbal > n))
          control <-
            appendDefArgs(control,
                          list(Rmax = 100L,
                               abstol = sd(scores[abs(scores) > 0]) / 100))
          eps <- 2 * control$abstol
          r <- 0L
          while (r < control$Rmax & eps > control$abstol) { 
            r <- r + 1L
            sAdd <- olmm_addScores(x)
            sR <- rbind(scores[subsObs,, drop = FALSE],
                        -sAdd$score_obs[, parm, drop = FALSE])
            if (r == 1L) # is always the same for a fixed data set
              sbj <- factor(c(as.character(subject[subsObs]),
                              as.character(sAdd$subject)),
                            levels = levels(subject))
            sR <- matrix(c(t(sR[order(sbj),,drop=FALSE])), Nbal*k,length(subsSbj))
            sR <- matrix(c(Ti %*% sR), Nbal * length(subsSbj), k, byrow = TRUE)
            sR <- sR[!subsAdd,,drop=FALSE]
            sR <- sTmp + sR
            if (r > 1) eps <- max(abs(sR / r - sTmp / (r - 1)))
            if (verbose && r > 1)
              cat("\nnit = ", r,
                  "max|diff| =", format(eps, digits = 3, scientific = TRUE))
            sTmp <- sR
          }
          if (eps <= control$abstol) {
            if (verbose)
              cat("\nalgorithm converged: nit =", r,
                  "max|diff| =", format(eps, digits = 3, scientific = TRUE), "\n")
          } else {
            if (!silent) warning("imputation algorithm did not converge.")
          }
          sT[rownames(sT)[subsOrd],] <- (sTmp / r)
          scores[] <- sT[]
          ## add information about the completing data method
          attr(scores, "conv.complete") <- as.integer(eps <= control$abstol)
          attr(scores, "nit.complete") <- r
          attr(scores, "eps.complete") <- eps
        }
      }
      
      ## add attributes about the transformation matrix
      attr(scores, "T") <- T
      attr(scores, "conv.T") <- attr(T, "conv")
      attr(scores, "nit.T") <- attr(T, "nit")
      attr(scores, "eps.T") <- attr(T, "eps")
      
    }

    if (!complete && any(!subset)) scores <- scores[subset,,drop=FALSE]
    
  } else {

    scores <- apply(scores, 2, tapply, subject, sum)
    
    if (!silent && complete)
      warning("handling unbalanced data for subject level scores has ",
              "yet not been implemented.")
    if (verbose) cat("OK")
  }

  ## append attributes
  defAttr <- list(T = matrix(0,k,k), conv.T = 1L, nit.T = 0L,
                  eps.T = 0.05, conv.complete = 1L, nit.complete = 0L,
                  eps.complete = 0.0, predecor = as.integer(predecor))
  attributes(scores) <- appendDefArgs(attributes(scores), defAttr)
  colnames(attr(scores, "T")) <- rownames(attr(scores, "T")) <- colnames(scores)

  x <- xold
  .Call("olmm_update_marg", x, slot(x, "coefficients"), PACKAGE = "vcrpart")
  
  if (verbose) cat("\n* return negative scores\n")
  
  ## return scores
  return(scores)
}
    

## thanks lme4
extractAIC.olmm <- function(fit, scale = 0, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L,"df")
  c(edf,-2*L + k*edf)
}

fitted.olmm <- function(object, ...) {
  return(predict(object, ...))
}

fixef.olmm <- function(object, which = c("all", "ce", "ge"), ...) {

  which <- match.arg(which)
  dims <- slot(object, "dims")
  coef <- coef(object)
  rval <- c()

  ## category-specific coefficients
  if (which %in% c("all", "ce") && dims["pCe"] > 0) {
    subs <- seq(from  = 1, to = dims["pCe"] * dims["nEta"], by = 1)
    rval <- c(rval, coef[subs])
  }

  ## global coefficients
  if (which %in% c("all", "ge") && dims["pGe"] > 0) { 
    subs <- seq(from = dims["pCe"] * dims["nEta"] + 1,
                to = dims["pCe"] * dims["nEta"] + dims["pGe"],
                by = 1)
    rval <- c(rval, coef[subs])
  }
  
  return(rval)
}

formula.olmm <- function(x, ...) as.formula(slot(x, "formula"), env = parent.frame())

gefp.olmm <- function(object, scores = NULL, predecor = TRUE,
                      order.by = NULL, parm = NULL, subset = NULL,
                      center = TRUE, drop = TRUE,
                      silent = FALSE, ...) {
  
  ## extract scores (if scores is not a matrix)
  if (is.null(scores)) {
    estfunArgs <-
      appendDefArgs(list(...),
                    list(silent = silent, force = force,
                         predecor = predecor, complete = TRUE))
    estfunArgs$x <- object
    scores <- try(do.call("estfun.olmm", estfunArgs))
  } else if (is.function(scores)) {    
    scores <- scores(object)
  } else if (is.matrix(scores)) {
    if (!silent && !is.null(attr(scores, "predecor")) &&
        as.integer(predecor) != attr(scores, "predecor"))
      warning("'scores' is not pre-decorrelated")
  }
    
  if (!is.matrix(scores)) stop("extracting the score function failed.")

  ## set 'order.by'
  if (is.null(order.by)) order.by <- 1:nobs(object)
  if (is.factor(order.by)) order.by <- droplevels(order.by)
  
  ## set subset
  if (!is.null(subset)) {
    if (is.character(subset)) subset <- rownames(scores) %in% subset
    if (is.numeric(subset)) subset <- (1:nobs(object)) %in% subset
  } else {
    subset <- rep(TRUE, nobs(object))
  }
  subsScores <- rownames(model.frame(object)) %in% rownames(scores)

  ## create process
  process <- scores[subset[subsScores],,drop=FALSE]
  cn <- colnames(process) 
  order.by <- order.by[subset & subsScores]
  
  ## get dimensions
  n <- nrow(scores)
  k <- ncol(scores)
  
  ## if necessary, subtract the column means
  if (center & max(abs(cMeans <- colMeans(process))) > 1e-6)
    process <- process - matrix(cMeans, nrow(process), ncol(process), byrow = TRUE)
  
  ## scale scores by the number of observations
  process <- process / sqrt(n)

  ## multiply scores with the inverse of the square root of their crossproduct
  J12Inv <- try(chol2inv(chol(root.matrix(crossprod(process)))), silent = TRUE)
  if (class(J12Inv) == "try-error" && drop) {
      J12 <- crossprod(process)
      subs <- rep(TRUE, k)
      subs[diag(J12) / max(diag(J12)) < 1e-2] <- FALSE
      J12Inv <- matrix(0, k, k)
      J12Inv[subs, subs] <- chol2inv(chol(root.matrix(J12[subs, subs, drop = FALSE])))
      if (!silent) warning("covariance matrix is not positive semidefinite. ",
                           "Omit terms: ",
                           paste(colnames(process)[!subs], collapse = ", "))
  } 
  process <- t(J12Inv %*% t(process))

  ## order and cumulate the process
  index <- order(order.by)
  process <- apply(process[index, , drop = FALSE], 2, cumsum)
  process <- rbind(0, process)
  colnames(process) <- cn

  ## transform process to a multivariate time series
  time <- order.by[index]
  if (is.factor(time)) time <- as.numeric(droplevels(time))
  time <- suppressWarnings(c(time[1] - as.numeric(diff(time[1:2])), time))

  ## extract terms
  if (!is.null(parm)) process <- process[, parm, drop = FALSE]

  ## return a list of class "gefp"
  rval <- list(process = suppressWarnings(zoo(process, time)),
               nreg = k,
               nobs = n,
               call = match.call(),
               fit = NULL,
               scores = NULL, 
               fitted.model = getCall(object),
               par = NULL,
               lim.process = "Brownian bridge",
               type.name = "M-fluctuation test",
               order.name = deparse(substitute(order.by)),
               subset <- rownames(model.frame(object))[subset & subsScores],
               J12 = NULL)
  class(rval) <- "gefp"
  return(rval)
}

getCall.olmm <- function(x, ...) {
  return(slot(x, "call"))
}

logLik.olmm <- function(object, ...) {
  dims <- slot(object, "dims")
  rval <- slot(object, "logLik")
  attr(rval, "nall") <- attr(rval, "nobs") <- dims[["n"]]
  nPar <- dims[["nPar"]] - (1 - dims[["hasRanef"]])
  attr(rval, "df") <- nPar
  class(rval) <- "logLik"
  return(rval)
}

model.frame.olmm <- function(formula, ...) slot(formula, "frame")

model.matrix.olmm <- function(object, which = c("fe", "fe-ce", "fe-ge",
                                        "re", "re-ce", "re-ge"), ...) {
  which <- match.arg(which)
  rval <- switch(substr(which, 1L, 2L),
                 fe = slot(object, "X"),
                 re = slot(object, "W"))
  if (!which %in% c("fe", "re")) {
    attr <- attributes(rval)
    subs <- attr(rval, "merge") == switch(substr(which, 4, 5), ce = 1, ge = 2)
    attr$assign <- attr$assign[subs]
    attr$merge <- attr$merge[subs]
    attr$orig.colnames <- attr$orig.colnames[subs]
    rval <- rval[, subs, drop = FALSE]
    attributes(rval) <- appendDefArgs(attributes(rval), attr)
  }
  return(rval)
}


nobs.olmm <- function(object, ...) slot(object, "dims")[["n"]]

predict.olmm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class", "ranef"),
                         ranef = FALSE, na.action = na.pass, ...) {

  ## extract data
  type <- match.arg(type) # retrieve type

  ## resolve conflicts with the 'ranef' argument
  if (type == "ranef" & !is.null(newdata))
    stop("prediction for random effects for 'newdata' is not implemented.")
  if (type == "ranef") return(ranef(object, ...))
  
  if (type == "prob") type <- "response"
  formList <- vcrpart_formula(formula(object)) # extract formulas
  offset <- list(...)$offset
  subset <- list(...)$subset
  dims <- slot(object, "dims")
  
  ## checks
  if (!is.null(newdata) && !class(newdata) == "data.frame")
    stop("'newdata' must be a 'data.frame'.")
  if (!class(ranef) %in% c("logical", "matrix"))
    stop("'ranef' must be a 'logical' or a 'matrix'.")
  
  if (dims["hasRanef"] < 1L) ranef <- FALSE
  
  if (is.null(newdata)) {
    
    ## extract data from the model fit
    X <- model.matrix(object, "fe")
    W <- model.matrix(object, "re")
    subject <- slot(object, "subject")
    offset <- slot(object, "offset")
    
    ## check and extract random effects
    if (is.logical(ranef)) {
      if (ranef) ranef <- ranef(object)
    } else {
      if (any(dim(ranef) != c(dims["N"], dims["qEta"])))
        stop("'ranef' matrix has wrong dimensions")
    }
    
  } else {
    
    ## data preparation
    if (is.matrix(ranef)) {
      mfForm <- formList$all # whole equation
    } else {
      ## fixed effects only
      getTerms <- function(x) attr(terms(x, keep.order = TRUE), "term.labels")
      terms <- lapply(formList$fe$eta, getTerms)
      mfForm <- formula(paste("~", paste(unlist(terms), collapse = "+")))
    }
    mf <- model.frame(object)
    Terms <- delete.response(terms(mfForm))
    xlevels <- .getXlevels(attr(mf, "terms"), mf)
    xlevels <- xlevels[names(xlevels) %in%  all.vars(Terms)]
    
    newdata <- as.data.frame(model.frame(Terms, newdata,
                                         na.action = na.action,
                                         xlev = xlevels))    
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, mf)   
 
    ## extract fixed effect model matrix from newdata
    X <- olmm_merge_mm(model.matrix(terms(formList$fe$eta$ce, keep.order = TRUE),
                                    newdata, attr(slot(object, "X"), "contrasts")),
                       model.matrix(terms(formList$fe$eta$ge, keep.order = TRUE),
                                    newdata, attr(slot(object, "X"), "contrasts")),
                       TRUE)
    rownames(X) <- rownames(newdata)

    ## delete columns of dropped terms
    X <- X[, colnames(slot(object, "X")), drop = FALSE]

    if (is.null(offset))
      offset <- matrix(0.0, nrow(X), dims["nEta"])
    
    if (is.logical(ranef) && ranef)
      ranef <- matrix(0.0, nrow(X), dims["qEta"])      
      
    ## random effects
    if (is.matrix(ranef)) {
      if (!slot(object, "subjectName") %in% colnames(newdata)) {
        stop(paste("column '", slot(object, "subjectName"),
                   "' not found in the 'newdata'.", sep = ""))
      } else {
        subject <- factor(newdata[, slot(object, "subjectName")])
      }
      
      ## extract model formulas
      W <- olmm_merge_mm(model.matrix(terms(formList$re$eta$ce, keep.order = TRUE),
                                      newdata, attr(slot(object, "W"), "contrasts")),
                         model.matrix(terms(formList$re$eta$ge, keep.order = TRUE),
                                      newdata, attr(slot(object, "W"), "contrasts")),
                                      FALSE)
      rownames(W) <- rownames(newdata)
      
      ## check entered random effects
      if (any(dim(ranef) != c(nlevels(subject), ncol(W))))
        stop("'ranef' matrix has wrong dimensions")
      
      if (any(!levels(subject) %in% rownames(ranef))) {
        stop(paste("random effects missing for subjects",
                   paste(setdiff(levels(subject), rownames(ranef)),
                         collapse = ", ")))
      } else {
        ranef <- ranef[levels(subject),, drop = FALSE]
      }
      
    } else {

      ## set random effects to zero
      W <- matrix(0.0, nrow(X), dims["q"])
      subject <- factor(rep(1, nrow(W)))
    }
  }
  
  if (!is.null(subset)) {
    X <- X[subset, , drop = FALSE]
    W <- W[subset, , drop = FALSE]
    if (!is.null(subject)) subject <- subject[subset]
    offset <- offset[subset,,drop = FALSE]
    if (is.matrix(ranef))
      ranef <- ranef[sort(unique(as.integer(subject))), , drop = FALSE]
    if (!is.null(subject)) subject <- factor(subject)
  }

  yLevs <- levels(slot(object, "y"))
  if (nrow(X) > 0) {

    ## compute linear predictor
    eta <- offset + X %*% slot(object, "fixef")
    
    ## predict marginal ...
    if (is.logical(ranef) && !ranef & type %in% c("response", "class")) {

      probs <- matrix(0, nrow(X), ncol(eta) + 1L)
      colnames(probs) <- levels(slot(object, "y"))
      rownames(probs) <- rownames(X)
      .Call("olmm_pred_marg", object, eta, W, nrow(X), probs, PACKAGE = "vcrpart")
      
    } else {
      ## or conditional probabilities (or linear predictor)
      
      ## extend linear predictor
      if (is.matrix(ranef))
        eta <- eta +
          matrix(rowSums(W * ranef[as.integer(subject), , drop = FALSE]),
                 nrow(X), ncol(eta))
      rownames(eta) <- rownames(X)
      colnames(eta) <- colnames(slot(object, "eta"))
        
      if (type == "link") return(eta)

      probs <- slot(object, "family")$linkinv(eta)
      colnames(probs) <- yLevs
      rownames(probs) <- rownames(X)
    }
    
    if (type == "response") {

      ## return probabilites
      rval <- probs
    } else if (type == "class") {
      
      ## return most probable category
      rval <- apply(probs, 1, which.max)
      rval <- ordered(yLevs[rval], levels = yLevs)
      names(rval) <- rownames(X)
    }
    
  } else { # useful for predict.tvcm
    if (type == "response") {
      rval <- matrix(, 0, dims["J"], dimnames = list(NULL, yLevs))
    } else if (type == "class") {
      rval <- NULL
    }
  }
  return(rval)
}

print.olmm <- function(x, ...) {

  so <- summary.olmm(x, silent = TRUE)
  
  if (length(so$methTitle) > 0) cat(so$methTitle, "\n\n")
  if (length(so$family) > 0) cat(" Family:", so$family, "\n")
  if (length(so$formula) > 0) cat("Formula:", so$formula, "\n")
  if (length(so$data) > 0) cat("   Data:", so$data, "\n")
  if (length(so$subset) > 0) cat(" Subset:", so$subset ,"\n")

  if (length(so$na.action) > 0L) { cat("\n"); cat(so$na.action, "\n"); }
  
  if (length(so$AICtab) > 0) {
    cat("\nGoodness of fit:\n")
    print(so$AICtab, ...)
  }

  if (length(so$REmat) > 0) {
    cat("\nRandom effects:\n")
    print.VarCorr.olmm(so$REmat, ...)
    cat(sprintf("Number of obs: %d, subjects: %d\n", so$dims["n"], so$dims["N"]))
  }

  if (length(so$feMatGe) > 0 && nrow(so$feMatGe) > 0) {
    cat("\nGlobal fixed effects:\n")
    print(fixef(x, "ge"), ...)
  }
  
  if (length(so$feMatCe) > 0 && nrow(so$feMatCe) > 0) {
    cat("\nCategory-specific fixed effects:\n")
    print(fixef(x, "ce"), ...)
  }
}

ranef.olmm <- function(object, norm = FALSE, ...) {
  if (!norm) {
    rval <- slot(object, "u") %*% t(slot(object, "ranefCholFac"))
  } else {
    rval <- slot(object, "u")
  }
  return(rval)
}

ranefCov.olmm <- function(object, ...) {

  dims <- slot(object, "dims")
    
  ## transformation for the adjacent-categories model
  if (slot(object, "family")$family == "adjacent" & dims["qCe"] > 0) {
    
    ranefCholFac <- rval <- slot(object, "ranefCholFac")
    
    ## row-wise subtraction of category-specific effects
    
    for (i in 1L:dims["qCe"]) {
      subs <- seq(i, dims["qCe"] * dims["nEta"], dims["qCe"])
      for (j in 1L:(dims["nEta"] - 1L)) {
        rval[subs[j], ] <- ranefCholFac[subs[j], ] -
          ranefCholFac[subs[j + 1L], ]
      }
    }
    
    for (i in 1L:(dims["nEta"] - 1L)) {
      rval[dims["qCe"]*(i-1)+1:dims["qCe"], ] <-
        ranefCholFac[dims["qCe"]*(i-1L)+1:dims["qCe"], ] -
          ranefCholFac[dims["qCe"]*i+1L:dims["qCe"], ]
    }
    return(rval %*% t(rval))
    
  } else {
    
    ## cumulative-link model or baseline-category model
    return(slot(object, "ranefCholFac") %*% t(slot(object, "ranefCholFac")))
  }
}


resid.olmm <- function(object, norm = FALSE, ...) {
  args <- list(...)
  args$object <- object
  args$type <- "response"
  fitted <- do.call("predict", args)
  y <- as.integer(model.response(model.frame(object)))
  J <- slot(object, "dims")["J"]
  n <- length(y)
  rval <- sapply(1:n, function(i) {
    sum(fitted[i, 1L:J > y[i]]) - sum(fitted[i, 1L:J < y[i]])
  })
  if (norm) {
    var <- (1.0 - apply(fitted^3, 1, sum)) / 3.0
    rval <- rval / sqrt(var)
  }
  return(rval)
}

residuals.olmm <- resid.olmm

setMethod(f = "show", signature = "olmm",
          definition = function(object) print.olmm(object))

simulate.olmm <- function(object, nsim = 1, seed = NULL,
                          newdata = NULL, ranef = TRUE, ...) {
  dotArgs <- list(...)
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1) 
  RNGstate <- .Random.seed
  dotArgs$type <- "response"
  pred <- predict(object, newdata = newdata, type = "prob", ranef = ranef, ...)
  FUN <- function(x) sample(levels(slot(object, "y")), 1, prob = x)
  rval <- as.data.frame(replicate(nsim, apply(pred, 1L, FUN)))
  for (i in 1:nsim)
    rval[,i] <- factor(rval[, i], levels = levels(slot(object, "y")), ordered = TRUE)
  if (nsim == 1) {
    colnames(rval) <- colnames(model.frame(object))[1]
  } else {
    colnames(rval) <- paste(colnames(model.frame(object))[1], 1:nsim, sep = ".")
  }
  attr(rval, "seed") <- RNGstate
  return(rval)
}

summary.olmm <- function(object, silent = FALSE, ...) {
            
  dims <- slot(object, "dims")
            
  ## goodness of fit measures
  lLik <- logLik(object)
  AICtab <- data.frame(AIC = AIC(lLik),
                       BIC = BIC(lLik),
                       logLik = as.vector(lLik),
                       deviance = deviance(object),
                       row.names = "")
  
  ## fixed-effect coefficients
  fixef <- fixef(object)
  
  ## fixed-effect coefficient-covariance matrix
  vcov <- try(vcov(object), silent = TRUE)
  validVcov <- class(vcov) != "try-error" && min(diag(vcov)) > 0
  if (!silent && !validVcov)
    warning("computation of variance-covariance matrix failed")
  
  ## global fixed effects
  if (dims["pGe"] > 0) {
    subs <- seq(dims["pCe"] * dims["nEta"] + 1L,
                dims["pCe"] * dims["nEta"] + dims["pGe"], 1L)
    feMatGe <- 
      cbind("Estimate" = fixef[subs],
            "Std. Error" = rep(NaN, length(subs)),
            "t value" = rep(NaN, length(subs)))
    if (validVcov) {
      feMatGe[, 2L] <- sqrt(diag(vcov)[subs])
      feMatGe[, 3L] <- feMatGe[, 1L] / feMatGe[, 2L]
    }
    
  } else { # empty matrix
    feMatGe <- matrix(, 0L, 3L, dimnames =
                      list(c(), c("Estimate", "Std. Error", "t value")))
  }
  
  ## category-specific fixed effects
  if (dims["pCe"] > 0L) {
    subs <- seq(1L, dims["pCe"] * dims["nEta"], 1L)
    feMatCe <-
      cbind("Estimate" = fixef[subs],
            "Std. Error" = rep(NaN, length(subs)),
            "t value" = rep(NaN, length(subs)))
    if (validVcov) {
      feMatCe[, 2L] <- sqrt(diag(vcov)[subs])
      feMatCe[, 3L] <- feMatCe[, 1L] / feMatCe[, 2L]
    }
  } else { # empty matrix
    feMatCe <- matrix(, 0L, 3L, dimnames =
                      list(c(), c("Estimate", "Std. Error", "t value")))
  }
  
  ## random effects
  if (dims["hasRanef"] > 0L) {
    VarCorr <- VarCorr(object)
  } else {
    VarCorr <-
      matrix(, 0L, 3L, dimnames = list(c(), c("Variance", "StdDev", "")))
  }
  
  ## title
  methTitle <- "Ordinal linear"
  if (dims["hasRanef"] > 0L) methTitle <- paste(methTitle, "mixed")
  methTitle <- paste(methTitle, "model")
  if (dims["hasRanef"] > 0L)
    paste(methTitle, " fit by marginal maximum\n",
          "likelihood with Gauss-Hermite quadrature", sep = "")
  family <- paste(slot(object, "family")$family, slot(object, "family")$link)

  na.action <- naprint(attr(model.frame(object), "na.action"))
  na.action <- if (na.action == "") character() else paste("(", na.action, ")", sep = "")
  call <- getCall(object)
  
  ## return a 'summary.olmm' object
  return(structure(
           list(methTitle = methTitle,
                family = family,
                formula = paste(deparse(formula(object)), collapse = "\n"),
                data = deparseCall(call$data),
                subset = deparseCall(call$subset),
                AICtab = AICtab,
                feMatCe = feMatCe,
                feMatGe = feMatGe,
                REmat = VarCorr,
                na.action = na.action,
                dims = dims,
                dotargs = list(...)), class = "summary.olmm"))
}

print.summary.olmm <- function(x, ...) {

  args <- appendDefArgs(list(...), x$dotargs)
  
  if (length(x$methTitle) > 0L) cat(x$methTitle, "\n\n")
  if (length(x$family) > 0L) cat(" Family:", x$family, "\n")
  if (length(x$formula) > 0L) cat("Formula:", x$formula, "\n")
  if (length(x$data) > 0L) cat("   Data:", x$data, "\n")
  if (length(x$subset) > 0L) cat(" Subset:", x$subset ,"\n")
  
  if (length(x$AICtab) > 0L) {
    cat("\nGoodness of fit:\n")
    args$x <- x$AICtab
    do.call("print", args)
  }
  
  if (length(x$REmat) > 0L) {
    cat("\nRandom effects:\n")
    args$x <- x$REmat
    do.call("print", args)
    cat(sprintf("Number of obs: %d, subjects: %d\n", x$dims["n"], x$dims["N"]))
  }

  if (length(x$na.action) > 0L) cat(x$na.action, "\n")
  
  if (length(x$feMatGe) > 0 && nrow(x$feMatGe) > 0L) {
    cat("\nGlobal fixed effects:\n")
    printCoefmat(x$feMatGe, digits = args$digits)
  }
  
  if (length(x$feMatCe) > 0 && nrow(x$feMatCe) > 0L) {
    cat("\nCategory-specific fixed effects:\n")
    printCoefmat(x$feMatCe, digits = args$digits)
  }
}

terms.olmm <- function(x, which = c("fe-ce", "fe-ge",
                            "re-ce", "re-ge"), ...) {
  which <- match.arg(which)
  which <- switch(which,
                  "fe-ge" = "feGe",
                  "fe-ce" = "feCe",
                  "re-ge" = "reGe",
                  "re-ce" = "reCe")
  return(slot(x, "terms")[[which]])
}


update.olmm <- function(object, formula., evaluate = TRUE, ...) {

  call <- getCall(object)
  if (is.null(call))
    stop("need an object with call slot")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }  
  if (evaluate)
    eval(call, parent.frame())
  else call
}

VarCorr.olmm <- function(x, sigma = 1., rdig = 3) {
            
  ## create formatted output according to VarCorr
  RECovMat <- ranefCov(x)
  REmat <- cbind(Variance = diag(RECovMat),
                 StdDev = sqrt(diag(RECovMat)))
  rval <- cbind(REmat, cov2cor(RECovMat))
  attr(rval, "title") <- paste("Subject:", slot(x, "subjectName"))
  class(rval) <- "VarCorr.olmm"
  return(rval)
}


print.VarCorr.olmm <- function(x, ...) { # S3 method

  rval <- format(x, ...)

  if (nrow(rval) > 1L) {

    ## prettify Corr matrix
    Corr <- rval[, 3L:ncol(rval)]
    Corr[upper.tri(Corr)] <- ""
    diag(Corr) <- colnames(Corr)
    Corr <- Corr[, -nrow(Corr), drop = FALSE]
    dimnames(Corr) <- NULL
    colnames(Corr) <- c("Corr", rep("", nrow(Corr) - 2L))
    rval <- cbind(rval[, 1L:2L, drop = FALSE], Corr)
  } else {

    ## random intercept models need no correlation terms
    rval <- rval[, 1L:2L, drop = FALSE]
  }
  
  ## print the output
  if (!is.null(attr(x, "title"))) {
    cat(attr( x, "title" ), "\n")
    attr(x, "title") <- NULL
  }
  print(rval, quote = FALSE)
  invisible(x)
}

vcov.olmm <- function(object, ...) {

  dims <- slot(object, "dims")
  info <- slot(object, "info")
  if (dims["hasRanef"] == 0L)
    info <- info[1L:dims["p"], 1:dims["p"]]
  
  ## extract inverse of negative info-matrix
  rval <- chol2inv(chol(-info)) 
  dimnames(rval) <- dimnames(info)
  
  ## parameter transformation for adjacent-category models
  if (slot(object, "family")$family == "adjacent") {
    
    ## matrix T with partial derivates of transformation 
    T <- diag(nrow(rval))
    subsRows <- seq(1, dims["pCe"] * (dims["nEta"] - 1L), 1L)
    subsCols <- seq(dims["pCe"] + 1L,
                    dims["pCe"] * dims["nEta"], 1L)
    if (length(subsRows) == 1L) {
      T[subsRows, subsCols] <- c(-1.0)
    } else {
      diag(T[subsRows, subsCols]) <- c(-1.0)
    }
    subsRows <- seq(dims["pCe"] + 1L,
                    dims["pCe"] * dims["nEta"], 1L)
    subsCols <- seq(1,dims["pCe"] * dims["nEta"], 1L)
    if (length(subsRows) == 1L) {
      T[subsRows, subsCols] <- c(-1.0)
    } else {
      diag(T[subsRows, subsCols]) <- c(-1.0)
    }
    
    ## transform covariance-matrix
    rval <- T %*% rval %*% t(T)
    dimnames(rval) <- dimnames(info)
  }
  
  return(rval)
}

weights.olmm <- function(object, level = c("observation", "subject"), ...) {
  return(switch(match.arg(level),
                observation = slot(object, "weights"),
                subject = slot(object, "weights_sbj")))
}
