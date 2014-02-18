## --------------------------------------------------------- #
## Author:          Reto Buergin, rbuergin@gmx.ch
## Date:            2013-09-18
##
## Description:
## methods for olmm objects.
##
## anova:       Likelihood-ratio tests for the comparison of
##              models
## coef, coefficients: Extract model coefficients
## deviance:    -2*Log-likelihood at the ML estimator
## estfun:      Negative scores
## fitted:      Extract fitted values from the model
## fixef:       Extract fixed effect parameters
## formula:     Extracts object@formula
## getCall:     Extracts object@call
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
## reweight:    Refits a model with new weights
## show:        Print summary output (method for olmm and
##              olmm.summary objects)
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
## 2013-09-06: modify formula() method. Now the formula slot
##             is called
## 2013-09-06: add S3 for terms() method
## 2013-09-06: add drop1() method
##
## To do:
## - implement new type of residuals (article ...)
## - improve update method, implement residuals for the
##   baseline-category- and the adjacent category model and
##   implement plot methods
## - predict marginal distribution
## --------------------------------------------------------- #


olmm_anova <- function(object, ...) {

  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- if (length(dots) > 0) {
    sapply(dots, is, "olmm") |
    sapply(dots, is, "polr") |
    sapply(dots, is, "vcolmm")
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

    stop("single argument anova for olmm objects' not yet implemented")
  }
}

setMethod(f = "anova",
          signature = "olmm",
          definition = olmm_anova)


olmm_coef <- function(object, ...) {
  dims <- object@dims
  if (dims["family"] == 3L) {
    T <- diag(dims["nPar"])
    subsRows <- seq(1, dims["pEtaVar"] * (dims["nEta"] - 1), 1)
    subsCols <- seq(dims["pEtaVar"] + 1,
                    dims["pEtaVar"] * dims["nEta"], 1)
    if (length(subsRows) == 1) {
      T[subsRows, subsCols] <- -1
    } else {
      diag(T[subsRows, subsCols]) <- -1
    }
    rval <- c(T %*% object@coefficients)
    names(rval) <- names(object@coefficients)
  } else {
    rval <- object@coefficients
  }
  if (dims["hasRanef"] == 0)
    rval <- rval[!grepl("ranefCholFac", names(rval))]
  
  return(rval)
}


setMethod(f = "coef",
          signature = "olmm",
          definition = olmm_coef)


setMethod(f = "coefficients",
          signature = "olmm",
          definition = olmm_coef)

setMethod(f = "deviance",
          signature = "olmm",
          definition = function(object, ...) {
            return(-2 * object@logLik)
          })


estfun.olmm <- function(x, level = c("observation", "subject"), ...) {
  level <- match.arg(level)
  if (level == "subject") return(-x@score_sbj) else return(-x@score_obs)
}


setMethod(f = "fitted", signature = "olmm",
          definition = function(object, ...) {
            return(predict(object, ...))
          })


olmm_fixef <- function(object, which = c("all", "npo", "po"), ...) {

  which <- match.arg(which)
  dims <- object@dims
  coef <- coef(object)
  rval <- c()

  ## predictor-variable coefficients
  if (which %in% c("all", "npo") && dims["pEtaVar"] > 0) {
    subs <- seq(from  = 1, to = dims["pEtaVar"] * dims["nEta"], by = 1)
    rval <- c(rval, coef[subs])
  }

  ## predictor-invariant coefficients
  if (which %in% c("all", "po") && dims["pEtaInv"] > 0) { 
    subs <- seq(from = dims["pEtaVar"] * dims["nEta"] + 1,
                to = dims["pEtaVar"] * dims["nEta"] + dims["pEtaInv"],
                by = 1)
    rval <- c(rval, coef[subs])
  }
  
  return(rval)
}


setMethod(f = "fixef",
          signature = "olmm",
          definition = olmm_fixef)


setMethod("formula",
          signature = "olmm",
	  function(x, ...) {
            formula(x@formula)
	  })


rhoTerm <- function(object, order.by) {
  subject <- as.integer(object@subject[order(order.by)])
  n <- length(subject)
  N <- nlevels(subject)
  Ti <- table(subject)
  Tit <- rep(0, length(Ti))
  FUN <- function(t) {
    Tit[subject[t]] <<- Tit[subject[t]] + 1
    term1 <- t^2 / n^2 * sum(Ti * (Ti - 1))
    term2 <- 2 * t^2 / n - 2 * t / n * sum(Ti[subject[1:t]])
    term3 <- sum(Tit^2) - t
    return(term1 + term2 + term3)
  }    
  return(sapply(1:length(subject), FUN))
}

gefp.olmm <- function(object, order.by, terms = NULL, subset = NULL) {
  psi <- estfun(object)
  if (is.null(terms)) terms <- 1:ncol(psi)
  psi <- psi[, terms, drop = FALSE]
  if (is.null(subset)) subset <- 1:nrow(psi)
  psi <- psi[subset, , drop = FALSE]
  order.by <- order.by[subset]
  subject <- object@subject[subset]
  if (is.factor(order.by)) order.by <- droplevels(order.by)
  J_obs <- as.matrix(crossprod(psi))
  J_sbj <- as.matrix(crossprod(apply(psi, 2, tapply, droplevels(subject), sum)))
  rho <- (J_sbj - J_obs) / sum(table(subject) * (table(subject) - 1))
  n <- nrow(psi)
  k <- ncol(psi)
  index <- order(order.by)
  process <- apply(psi[index, , drop = FALSE], 2, cumsum)
  if (object@dims["hasRanef"] > 0L) {
    term <- rhoTerm(object, order.by)
    J <- sapply(1:n, function(t) J_obs / n + n / (t * (n - t)) * rho * term[t])
    J <- array(J, dim = c(ncol(psi), ncol(psi), n))
    J[, , n] <- J[, , 1]
    for (i in 1:n) {
      J12i <- try(chol2inv(chol(root.matrix(as.matrix(J[,,i])))), silent = TRUE)
      if (class(J12i) == "try-error") {
        process[i, ] <- NA
      } else {
        process[i, ] <- J12i %*% process[i, ] / sqrt(n)
      }
    }
  } else {
    J12 <- root.matrix(J_obs / n)
    process <- t(chol2inv(chol(J12)) %*% t(process)) / sqrt(n)
  }
  process <- rbind(0, process)
  time <- order.by[index]
  if (is.factor(time)) time <- as.numeric(droplevels(time))
  time <- suppressWarnings(c(time[1] - as.numeric(diff(time[1:2])), time))
  rval <- list(process = suppressWarnings(zoo(process, time)),
               nreg = k,
               nobs = n,
               call = match.call(),
               fit = NULL,
               scores = NULL,
               fitted.model = NULL,
               par = NULL,
               lim.process = "Brownian bridge",
               type.name = "M-fluctuation test",
               order.name = deparse(substitute(order.by)),
               J12 = NULL)
  class(rval) <- "gefp"
  return(rval)
}

setMethod(f = "getCall",
          signature = "olmm",
          definition = function(x, ...) {
            return(x@call)
            })

setMethod(f = "logLik",
          signature = "olmm",
          definition = function(object, ...) {
            dims <- object@dims
            rval <- object@logLik
            attr(rval, "nall") <- attr(rval, "nobs") <- dims[["n"]]
            nPar <- dims[["nPar"]] - (1 - dims[["hasRanef"]])
            attr(rval, "df") <- nPar
            class(rval) <- "logLik"
            return(rval)
          })

setMethod(f = "model.frame", signature(formula = "olmm"),
	  function(formula, ...) formula@frame)


setMethod(f = "model.matrix", signature(object = "olmm"),
	  function(object, which = c("fixef", "ranef"), ...) {
            return(switch(match.arg(which),
                          fixef = object@X,
                          ranef = object@W))
          })

neglogLik.olmm <- function(object, ...)
    return(-as.numeric(logLik(object)))


olmm_predict <- function(object, newdata = NULL,
                         type = c("prob", "class", "link"),
                         na.action = na.pass, offset = NULL,
                         ranef = FALSE, subset = NULL, ...) {
  
  type <- match.arg(type) # retrieve type
  form <- vcolmm:::olmm_formula(object@formula) # extract formulas

  if (object@dims["hasRanef"] < 1L) ranef <- FALSE
  
  if (is.null(newdata)) {
    
    ## extract data from the model fit
    X <- model.matrix(object, "fixef")
    W <- model.matrix(object, "ranef")
    subject <- object@subject
    offset <- object@offset
    
    ## check and extract random effects
    if (is.logical(ranef)) {
      if (ranef) ranef <- ranef(object)
    } else {
      stop("ranef must be a 'logical'.")
    }
    
  } else {
        
    ## data preparation
    if (is.matrix(ranef)) {
      mfForm <- form$full # whole equation
    } else {
      mfForm <- form$fixef # fixed effects only
    } 
    mf <- model.frame(mfForm, model.frame(object))
    Terms <- delete.response(attr(mf, "terms"))
    xlevels <- .getXlevels(attr(mf, "terms"), mf)
    if (is.matrix(ranef)) { # delete terms of subject vector
      xlevels <- xlevels[names(xlevels) != form$subjectName]
    }    
    newdata <- as.data.frame(model.frame(Terms, newdata,
                                         na.action = na.action,
                                         xlev = xlevels))    
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, mf)   
    
    ## extract fixed effect model matrix from newdata
    X <- vcolmm:::olmm_mergeMm(model.matrix(terms(form$fixefEtaVar),
                                   newdata, attr(object@X, "contrasts")),
                      model.matrix(terms(form$fixefEtaInv), newdata,
                                   attr(object@X, "contrasts")), TRUE)
    rownames(X) <- rownames(newdata)

    ## delete columns of dropped terms
    X <- X[, colnames(object@X), drop = FALSE]

    if (is.logical(ranef) && ranef)
      stop("enter a matrix with random effects") 
    
    ## random effects
    if (is.matrix(ranef)) {
      if (!form$subjectName %in% colnames(newdata)) {
        stop(paste("column '", form$subjectName, "' not found in the 'newdata'.", sep = ""))
      } else {
        subject <- factor(newdata[, form$subjectName])
      }
      
      ## extract model formulas
      W <- vcolmm:::olmm_mergeMm(model.matrix(terms(form$ranefEtaVar), newdata,
                                     attr(object@W, "contrasts")),
                        model.matrix(terms(form$ranefEtaInv), newdata,
                                     attr(object@W, "contrasts")), FALSE)
      rownames(W) <- rownames(newdata)
      
      ## check entered random effects
      if (any(dim(ranef) != c(nlevels(subject), ncol(W))))
        stop("ranef matrix has wrong dimensions")

      if (any(!levels(subject) %in% rownames(ranef))) {
        stop(paste("random effects missing for subjects",
                   paste(setdiff(levels(subject), rownames(ranef)),
                         collapse = ", ")))
      } else {
        ranef <- ranef[levels(subject),, drop = FALSE]
      }
            
    } else {

      ## set random effects to zero
      W <- matrix(0.0, nrow(X), object@dims["q"])
      subject <- factor(rep(1, nrow(W)))
    }
    if (is.null(offset)) offset <- rep(0.0, nrow(newdata))
  }
  
  if (!is.null(subset)) {
    X <- X[subset, , drop = FALSE]
    W <- W[subset, , drop = FALSE]
    if (!is.null(subject)) subject <- subject[subset]
    offset <- offset[subset]
    if (is.matrix(ranef))
      ranef <- ranef[sort(unique(as.integer(subject))), , drop = FALSE]
    if (!is.null(subject)) subject <- factor(subject)
  }
    
  n <- nrow(X)
  J <- object@dims["J"]
  
  if (n > 0) {

    ## compute linear predictor
    eta <- X %*% object@fixef
    if (length(offset) > 0)
      eta <- eta + matrix(offset, n, object@dims["nEta"])
    
    ## predict marginal ...
    if (is.logical(ranef) && !ranef & type %in% c("prob", "class")) {

      probs <- matrix(0, nrow(X), object@dims["J"])
      colnames(probs) <- levels(object@y)
      rownames(probs) <- rownames(X)
      .Call("olmm_pred_marg", object, eta, W, n, probs, PACKAGE = "vcolmm")
      
      ## or conditional probabilities (or linear predictor)
    } else {
            
      ## compute linear predictor
      if (is.matrix(ranef))
        eta <- eta + matrix(rowSums(W * ranef[as.integer(subject), , drop = FALSE]), n, object@dims["nEta"])
      rownames(eta) <- rownames(X)
      colnames(eta) <- colnames(object@eta)
        
      if (type == "link") return(eta)
        
      if (object@dims["family"] == 1L) {    
        linkFUN <- switch(object@dims["link"], plogis, pnorm, pcauchy)
        
        cumProbs <- linkFUN(eta)
        cumProbs <- cbind(cumProbs, rep(1, n))
        tmp <- apply(cumProbs, 1, diff)
        if (ncol(cumProbs) > 2) tmp <- t(tmp)
        probs <- cbind(cumProbs[, 1], tmp)
        
      } else if (object@dims["family"] %in% c(2L, 3L)) {
            
        probs <- exp(eta) / (1 + matrix(rowSums(exp(eta)),
                                        n, object@dims["nEta"]))
        probs <- cbind(probs, 1 - rowSums(probs))
      }
      colnames(probs) <- levels(object@y)
      rownames(probs) <- rownames(X)
    }
        
    if (type == "prob") {
      ## return probabilites
      rval <- probs
    } else if (type == "class") {
      ## return most probable category
      rval <- apply(probs, 1, which.max)
      rval <- factor(levels(object@y)[rval], levels = levels(object@y),
                     ordered = TRUE)
      names(rval) <- rownames(X)
    }
    
  } else { # useful for predict.vcolmm
    if (type == "prob") {
      rval <- matrix(, 0, J, dimnames = list(NULL, levels(object@y)))
    } else if (type == "class") {
      rval <- NULL
    }
  }
  return(rval)
}


setMethod(f = "predict",
          signature = "olmm",
          definition = olmm_predict)


olmm_print <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  if (length(x@methTitle) > 0)
    cat(x@methTitle, "\n\n")
  
  if (length(x@formula) > 0)
    cat("Formula:", x@formula, "\n")

  if (length(x@data) > 0)
    cat("   Data:", x@data, "\n")

  if (length(x@subset) > 0)
    cat(" Subset:", x@subset ,"\n")

  if (length(x@AICtab) > 0) {
    cat("\nGoodness of fit:\n")
    print(x@AICtab, digits = digits)
  }

  if (length(x@REmat) > 0) {
    cat("\nRandom effects:\n")
    print.VarCorr.olmm(x@REmat, digits = digits, ...)
    cat(sprintf("Number of obs: %d, subjects: %d\n", x@dims["n"],
                x@dims["N"]))
  }

  if (length(x@FEmatEtaInv) > 0 && nrow(x@FEmatEtaInv) > 0) {
    cat("\nPredictor-invariant fixed effects:\n")
    printCoefmat(x@FEmatEtaInv, digits = digits)
  }
  
  if (length(x@FEmatEtaVar) > 0 && nrow(x@FEmatEtaVar) > 0) {
    cat("\nPredictor-variable fixed effects:\n")
    printCoefmat(x@FEmatEtaVar, digits = digits)
  }
}


setMethod(f = "print", signature = "olmm",
          definition = function(x, ...) olmm_print(summary(x), ...))


setMethod(f = "print", signature = "summary.olmm", definition = olmm_print)


setMethod(f = "ranef", signature = "olmm",
          definition = function(object, norm = FALSE, ...) {
            if (!norm) {
              rval <- object@u %*% t(object@ranefCholFac)
            } else {
              rval <- object@u
            }
            return(rval)
          })


setMethod(f = "ranefCov",
          signature = "olmm",
          definition = function(object, ...) {

            ## transformation for the adjacent-categories model
            if (object@dims["family"] == 3L & object@dims["qEtaVar"] > 0) {

              dims <- object@dims
              ranefCholFac <- rval <- object@ranefCholFac

              ## row-wise subtraction of predictor-variable effects
              
              for (i in 1:dims["qEtaVar"]) {
                  subs <- seq(i, dims["qEtaVar"] * dims["nEta"], dims["qEtaVar"])
                  for (j in 1:(dims["nEta"] - 1)) {
                      rval[subs[j], ] <- ranefCholFac[subs[j], ] -
                          ranefCholFac[subs[j + 1], ]
                  }
              }
              
              for (i in 1:(dims["nEta"] - 1)) {
                rval[dims["qEtaVar"]*(i-1)+1:dims["qEtaVar"], ] <-
                  ranefCholFac[dims["qEtaVar"]*(i-1)+1:dims["qEtaVar"], ] -
                    ranefCholFac[dims["qEtaVar"]*i+1:dims["qEtaVar"], ]
              }
              return(rval %*% t(rval))
              
            } else {

              ## cumulative-link model or baseline-category model
              return(object@ranefCholFac %*% t(object@ranefCholFac))
            }
          })


olmm_resid <- function(object, norm = FALSE, ...) {
  fitted <- predict(object, ...)
  y <- as.integer(model.response(model.frame(object)))
  J <- object@dims["J"]
  n <- length(y)
  rval <- sapply(1:n, function(i) {
    sum(fitted[i, 1:J > y[i]]) - sum(fitted[i, 1:J < y[i]])
  })
  if (norm) {
    var <- (1 - apply(fitted^3, 1, sum)) / 3
    rval <- rval / sqrt(var)
  }
  return(rval)
}


setMethod(f = "resid",
          signature = "olmm",
          definition = olmm_resid)


setMethod(f = "residuals",
          signature = "olmm",
          definition = olmm_resid)


olmm_reweight <- function(object, weights, verbose, silent = FALSE, ...) {

  rval <- object
  
  ## check verbose argument  
  if (missing(verbose)) {
    verbose <- object@dims["verb"] == 1
  } else if (!is.logical(verbose)) {
    stop("verbose argument must be logical")
  }
  
  ## checking new weights
  
  if (verbose) cat("* checking new weights ... ")
  
  ## return if the new weights concord with the old weights
  if (identical(weights, weights(object))) {
    if (!silent) warning("new weights confound with the current weights. Return the old object")
    return(object)
  }

  ## ... to get sure
  weights <- as.double(weights)
  
  ## stop if negative weights
  if (any(weights < 0)) stop("negative weights")
  
  ## weights for subjects
  FUN <- function(x) {
    rval <- unique(x)
    if (length(rval) > 1)
      rval <- rval[rval != 0] # allows the omitting of an observation
    return(rval)
  }
  weights_sbj <- tapply(weights, object@subject, FUN)
  if (is.list(weights_sbj)) {
    stop("weights must be constant within subjects")
  } else {
    weights_sbj <- c(weights_sbj)
  }
  
  ## update model frames  
  object@frame <- object@frame[weights > 0, , drop = FALSE]
  object@y <- object@y[weights > 0]
  object@X <- object@X[weights > 0, , drop = FALSE]
  object@W <- object@W[weights > 0, , drop = FALSE]
  object@subject <- factor(object@subject[weights > 0])
  object@weights <- weights[weights > 0]
  object@weights_sbj <- weights_sbj[weights_sbj> 0]
  object@offset <- object@offset[weights > 0]
  object@dims["n"] <- nrow(object@X)
  object@dims["N"] <- nlevels(object@subject)
  object@eta <- object@eta[weights > 0, , drop = FALSE]
  object@u <- object@u[weights_sbj > 0, , drop = FALSE]
  object@logLik_sbj <- object@logLik_obs[weights > 0]
  object@logLik_sbj <- object@logLik_sbj[weights_sbj > 0]
  object@score_obs <- object@score_obs[weights > 0, , drop = FALSE]
  object@score_sbj <- object@score_sbj[weights_sbj > 0, , drop = FALSE]
  
  ## modify optim slot

  if (verbose) cat("OK\n* setting up the fitting environment ... ")

  optim <- object@optim
  optim[[1]] <- object@coefficients
  environment(optim[[2]]) <- environment()
  if (!object@dims["numGrad"]) environment(optim[[3]]) <- environment()
  
  ## fitting the model

  if (verbose) cat("OK\n* fitting the model ... ")

  if (!is.null(optim$control$trace) &&
      optim$control$trace > 0)
      cat("\n")
  
  ## extract the function for fitting the model
  if (optim$method %in% c("nlminb", "ucminf")) {
      FUN <- optim$method
      optim <- optim[-which(names(optim) == "method")] 
  } else {
      FUN <- "optim"
  }
    
  systemTime <- system.time(object@output <- do.call(FUN, optim))
  
  ## print messages for opimization
  if (verbose) {
    cat(paste("OK\n\toptimization time:",
              signif(systemTime[3], 3),
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
    if (verbose) cat("\n* computing the approximative hessian matrix ... ")
    object@info[] <- # replace the info slot
      - hessian(func = optim[[2]], x = object@coefficients,
                method.args = list(r = 2, show.details = TRUE))
  }

  if (verbose) {
      eigenHess <- eigen(object@info, only.values = TRUE)$values
      condHess <- abs(max(eigenHess) / min(eigenHess))
      cat("\n\tcondition number of Hessian matrix:",
          format(condHess, digits = 2, scientific = TRUE))
    }

  ## prepare output model
  rval@weights <- weights
  rval@weights_sbj <- weights_sbj
  .Call("olmm_update_marg", rval, object@output$par, PACKAGE = "vcolmm")
  if (object@dims["numHess"] == 1L) rval@info[] <- object@info
  
  ## compute expected standardized random effects
  if (verbose) cat("\n* predicting random effects ... ")
  .Call("olmm_update_u", rval, PACKAGE = "vcolmm")

  ## reset environment of estimation equations
  environment(object@optim[[2]]) <- baseenv()
  if (!object@dims["numGrad"]) environment(object@optim[[3]]) <- baseenv()
  
  if (verbose) cat("OK\n* computations finished, return model object\n")
  
  return(rval)
}


setMethod(f = "reweight",
          signature = "olmm",
          definition = olmm_reweight)


setMethod(f = "show", signature = "olmm",
          definition = function(object) olmm_print(summary(object)))


setMethod(f = "show", signature = "summary.olmm",
          definition = function(object) olmm_print(object))


setMethod(f = "summary",
          signature = "olmm",
          definition = function(object, silent = FALSE, ...) {
            
            dims <- object@dims

            ## goodness of fit measures
            lLik <- logLik(object)
            AICframe <- data.frame(AIC = AIC(lLik),
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

            ## predictor-invariant fixed effects
            if (dims["pEtaInv"] > 0) {
              subs <- seq(dims["pEtaVar"] * dims["nEta"] + 1,
                          dims["pEtaVar"] * dims["nEta"] + dims["pEtaInv"],
                          1)
              FEmatEtaInv <- 
                cbind("Estimate" = fixef[subs],
                      "Std. Error" = rep(NaN, length(subs)),
                      "t value" = rep(NaN, length(subs)))
              if (validVcov) {
                FEmatEtaInv[, 2] <- sqrt(diag(vcov)[subs])
                FEmatEtaInv[, 3] <- FEmatEtaInv[, 1] / FEmatEtaInv[, 2]
              }
              
            } else { # empty matrix
              FEmatEtaInv <- matrix(, 0, 3, dimnames = list(c(), c("Estimate", "Std. Error", "t value")))
            }
            
            ## predictor-variable fixed effects
            if (dims["pEtaVar"] > 0) {
              subs <- seq(1, dims["pEtaVar"] * dims["nEta"], 1)
              FEmatEtaVar <-
                cbind("Estimate" = fixef[subs],
                      "Std. Error" = rep(NaN, length(subs)),
                      "t value" = rep(NaN, length(subs)))
              if (validVcov) {
                FEmatEtaVar[, 2] <- sqrt(diag(vcov)[subs])
                FEmatEtaVar[, 3] <- FEmatEtaVar[, 1] / FEmatEtaVar[, 2]
              }
              rownames(FEmatEtaVar) <- sub(pattern = "(Intercept)",
                                           replacement = "(Int)",
                                           x = rownames(FEmatEtaVar),
                                           fixed = TRUE)
            } else { # empty matrix
              FEmatEtaVar <- matrix(, 0, 3, dimnames = list(c(), c("Estimate", "Std. Error", "t value")))
            }
            
            ## random effects
            if (dims["hasRanef"] > 0) {
              VarCorr <- unclass(VarCorr(object))
            } else {
              VarCorr <- matrix(, 0, 3, dimnames = list(c(), c("Variance", "StdDev", "")))
            }
            ## title
            methTitle <- c("Cumulative", "Baseline",
                           "Adjacent-Categories")[dims["family"]]
            if (dims["family"] == 1L)
              methTitle <- paste(methTitle,
                                 c("Logit", "Probit",
                                   "Cauchy")[dims["link"]])
            if (dims["hasRanef"] > 0) {
              methTitle <- paste(methTitle,
                                 " Mixed-Effect Model fit by the\n",
                                 "Gauss-Hermite approximation of the\n",
                                 "Marginal Likelihood", sep = "")
            } else {
              methTitle <- paste(methTitle, "Model")
            }
            
            ## build a summary.olmm object
            new(Class = "summary.olmm",
                methTitle = methTitle,
                formula = if (!is.null(object@call$formula)) deparse(object@call$formula, nlines = 1, width.cutoff = 45)[1] else character(),
                data = if (!is.null(object@call$data)) deparse(object@call$data, nlines = 1, width.cutoff = 45)[1] else character(),
                subset = if (!is.null(object@call$subset)) deparse(object@call$subset, nlines = 1, width.cutoff = 45)[1] else character(),
                AICtab = AICframe,
                FEmatEtaVar = FEmatEtaVar,
                FEmatEtaInv = FEmatEtaInv,
                REmat = VarCorr,
                dims = dims)
          })


setMethod(f = "terms", signature(x = "olmm"),
	  function(x, which = c("fixefEtaInv", "fixefEtaVar",
                        "ranefEtaInv", "ranefEtaVar"), ...) {
            which <- match.arg(which)
            return(x@terms[[which]])
          })


olmm_update <- function(object, formula., evaluate = TRUE, ...) {

  call <- object@call
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

setMethod(f = "update",
          signature = "olmm",
          definition = olmm_update)


setMethod(f = "VarCorr",
          signature = "olmm",
          definition = function(x, sigma = 1., rdig = 3) {
            
            ## create formatted output according to VarCorr
            RECovMat <- ranefCov(x)
            REmat <- cbind(Variance = diag(RECovMat),
                          StdDev = sqrt(diag(RECovMat)))
            rval <- cbind(REmat, cov2cor(RECovMat))
            attr(rval, "title") <- paste("Subject:", x@subjectName)
            class(rval) <- "VarCorr.olmm"
            return(rval)
          })


print.VarCorr.olmm <- function(x, ...) { # S3 method

  rval <- format(x, ...)

  if (nrow(rval) > 1) {

    ## prettify Corr matrix
    Corr <- rval[, 3:ncol(rval)]
    Corr[upper.tri(Corr)] <- ""
    diag(Corr) <- colnames(Corr)
    Corr <- Corr[, -nrow(Corr), drop = FALSE]
    dimnames(Corr) <- NULL
    colnames(Corr) <- c("Corr", rep("", nrow(Corr) - 2))
    rval <- cbind(rval[, 1:2, drop = FALSE], Corr)
  } else {

    ## random intercept models need no correlation terms
    rval <- rval[, 1:2, drop = FALSE]
  }
  
  ## print the output
  if (!is.null(attr(x, "title"))) {
    cat(attr( x, "title" ), "\n")
    attr(x, "title") <- NULL
  }
  print(rval, quote = FALSE)
  invisible(x)
}


setMethod(f = "vcov",
          signature = "olmm",
          definition = function(object, ...) {

            dims <- object@dims
            info <- object@info
            if (dims["hasRanef"] == 0L)
              info <- object@info[1:dims["p"], 1:dims["p"]]
              
            ## extract inverse of negative info-matrix
            rval <- chol2inv(chol(-info)) 
            dimnames(rval) <- dimnames(info)

            ## parameter transformation for adjacent-category models
            if (dims["family"] == 3L) {
              
              ## matrix T with partial derivates of transformation
              T <- diag(dims["nPar"])
              subsRows <- seq(1, dims["pEtaVar"] * (dims["nEta"] - 1), 1)
              subsCols <- seq(dims["pEtaVar"] + 1,
                              dims["pEtaVar"] * dims["nEta"], 1)
              if (length(subsRows) == 1) {
                T[subsRows, subsCols] <- -1
              } else {
                diag(T[subsRows, subsCols]) <- -1
              }
              subsRows <- seq(dims["pEtaVar"] + 1,
                              dims["pEtaVar"] * dims["nEta"], 1)
              subsCols <- seq(1,dims["pEtaVar"] * dims["nEta"], 1)
              if (length(subsRows) == 1) {
                T[subsRows, subsCols] <- -1
              } else {
                diag(T[subsRows, subsCols]) <- -1
              }
              
              ## transform covariance-matrix
              rval <- T %*% rval %*% t(T)
              dimnames(rval) <- dimnames(object@info)
            }
            
            return(rval)
          })


setMethod(f = "weights",
          signature = signature(object = "olmm"),
          definition = function(object,
            level = c("observation", "subject"), ...) {
            return(switch(match.arg(level),
                          observation = object@weights,
                          subject = object@weights_sbj))
          })

terms.olmm <- function(x, ...) terms(x, ...)

## thanks lme4
extractAIC.olmm <- function(fit, scale = 0, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L,"df")
  c(edf,-2*L + k*edf)
}

## thanks lme4 (is modified)
drop1.olmm <- function(object, scope, scale = 0, test = c("none", "Chisq"),
                       k = 2, trace = FALSE, ...) {

  termsFixefEtaInv <- terms(object, "fixefEtaInv")
  tl <- factor.scope(attr(termsFixefEtaInv, "factor"),
                     list(drop = numeric()))$drop
  termsFixefEtaVar <- terms(object, "fixefEtaVar")
  tl <- c(tl, factor.scope(attr(termsFixefEtaVar, "factor"),
                           list(drop = numeric()))$drop)

  ff <- Formula(object@formula)
  
  if (missing(scope)) {
    scope <- tl
  } else {
    if (!is.character(scope)) {
      newff <- update(ff, scope)
      scope <- setdiff(tl, all.vars(newff))
    }
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L,
                dimnames =  list(c("<none>", scope), c("df", "AIC")))
  ans[1, ] <- extractAIC.olmm(object, scale, k = k, ...)
  ## BMB: avoid nobs, to avoid dependence on 2.13
  ## n0 <- nobs(object, use.fallback = TRUE)
  
  n0 <- nrow(object@frame)
  env <- environment(object@formula)
  ff <- Formula(as.formula(formula(object)))
  for (i in seq(ns)) {
    tt <- scope[i]
    if(trace > 1) {
      cat("trying -", tt, "\n", sep='')
      utils::flush.console()
    }
    if (length(ff)[2] == 2L &&
        length(all.vars(formula(ff, lhs = 0L, rhs = 2L))) > 0) {
      newff <- update(ff, as.formula(paste(". ~ . -", tt, "| . - ", tt)))
    } else {
      newff <- update(ff, as.formula(paste(". ~ . -", tt)))
    }
    nfit <- update(object, evaluate = FALSE)
    nfit$formula <- formula(newff)
    nfit <- eval(nfit, envir = env) # was  eval.parent(nfit)
    ans[i+1, ] <- extractAIC.olmm(nfit, scale, k = k, ...)
    ## BMB: avoid nobs, to avoid dependence on 2.13
    ## nnew <- nobs(nfit, use.fallback = TRUE)
    nnew <- nrow(nfit@frame)
    if(all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  
  dfs <- ans[1L , 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[,2])
  test <- match.arg(test)

  if (test == "Chisq") {
    dev <- ans[, 2L] - k*ans[, 1L]
    dev <- dev - dev[1L] ; dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    ## BMB: hack to extract safe_pchisq
    P[nas] <- stats:::safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
  } else if (test == "F") {
    stop("F test STUB -- unfinished maybe forever")
    dev <- ans[, 2L] - k*ans[, 1L]
    dev <- dev - dev[1L] ; dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    ## BMB: hack to extract safe_pchisq
    P[nas] <- stats:::safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(F)")] <- list(dev, P)
  }
  
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
            if(scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  return(aod)
}
