## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-02
##
## Description:
## Cross-validation, stability paths and random forests for
## 'tvcm' objects.
##
## Contents:
## AIC.tvcm:
## print.AIC.tvcm:
## BIC.tvcm
## oobrisk.tvcm
## cvfolds:             create cross-validation folds
## cvrisk.tvcmm:        cross-validation for 'tvcm' objects
## print.cvrisk.tvcm:   print for 'cv.tvcm' objects
## plot.cvrisk.tvcm:    plot fot 'cv.tvcm' objects
## stabpath.tvcm:       stability selection for 'tvcm' objects
## print.stabpath.tvcm: print for 'stabpath.tvcm' objects
## plot.stabpath.tvcm:  plot for 'stabpath.tvcm' objects
##
## Last modifications:
##
## To do:
## --------------------------------------------------------- #

AIC.tvcm <- function(object, sub = FALSE, maxalpha = NULL, maxwidth = NULL,
                        k = 2,  verbose = FALSE, ...) { 

  criterion <- list(...)$criterion
  if (is.null(criterion)) criterion <- "AIC"
  
  if (!sub) return(AIC(extract(object, "model"), k = k))
  
  control <- object$info$control
  control$verbose <- FALSE
  data <- model.frame(object)
  n <- nrow(data)
  partvar <- colnames(object$data)
  curDf <- attr(logLik(object), "df")
  if (!is.null(maxalpha)) control$alpha <- maxalpha
  if (!is.null(maxwidth)) control$maxwidth <- maxwidth
  
  ## refit the model if necessary
  if (any(c(control$alpha != object$info$control$alpha,
            control$maxwidth != object$info$control$maxwidth))) {
    
    ## prepare formulas and calls
    args <- list()
    args$formula <- formula(object, "original")
    args$control <- control
    args$fit <- object$info$fit
    args$family <- object$info$family
    args <- append(args, object$info$dotargs)
    args$data <- data
    args$weights <- weights(object)
    if (verbose) cat("* refitting the model ...")
    object <- try(do.call("tvcm", args), silent = TRUE)
    if (class(object)[1] == "try-error") stop("refitting failed.")
    if (verbose) cat("OK\n")
  }
  
  ## get steps to evaluate
  nstep <- object$info$nstep
  alpha <- width <- NULL

  splitpath <- object$info$splitpath
  steps <- rev(unlist(lapply(splitpath, function(x) x$step)))
  
  if (control$method == "mob") {
    alpha <- c(extract(object, "p.value", step = steps), 0.0)[-1L]
    if (length(alpha) > 1)
      for (i in 1:(length(alpha) - 1L))
        if (any(alpha[(i + 1):length(alpha)] > alpha[i])) steps <- steps[-i]
    alpha <- rev(alpha)[steps]
  } else {
    width <- steps + 1L
  }
  
  logLik <- unlist(lapply(splitpath[steps], function(x) x$logLik))
  df <- unlist(lapply(splitpath[steps], function(x) attr(x$logLik, "df")))
  
  ## output
  AICtab <- data.frame(Df = df, AIC = - 2 * logLik + k * df)
  colnames(AICtab)[2L] <- criterion
  if (control$method == "mob") {
    maxpar <- control$alpha
    parname <- "alpha"
    rownames(AICtab) <- paste("alpha", format(alpha, digits = 2))
  } else {
    maxpar <- control$maxwidth
    parname <- "maxwidth"
    rownames(AICtab) <- paste("width", width)
  }
  if (any(subs <- AICtab$Df == curDf))
    rownames(AICtab)[subs] <- paste(rownames(AICtab)[subs], "<current>")
  if (verbose) {
    cat("\n\n"); print(AICtab[order(AICtab[,criterion]),,drop=FALSE]);
  }
  par <- if (control$method == "mob") alpha else width
  bestpar <- unlist(par[which.min(AICtab[,criterion])])
  rval <- list(AICtab = AICtab,
               call = deparseCall(getCall(object)),
               criterion = criterion,
               par = par,
               bestpar = bestpar,
               maxpar = maxpar,
               parname = parname)
  class(rval) <- "AIC.tvcm"
  return(rval)
}

print.AIC.tvcm <- function(x, ...) {
    cat("In-Sample", x$criterion, "\n\n")
    cat("Call:", x$call, "\n\n")
    print(x$AICtab[order(x$AICtab[,x$criterion]),,drop=FALSE], ...)
    cat(paste("\nOptimal ",x$parname, ": ",
              format(x$bestpar, digits = 2), "\n", sep = ""))
    return(invisible(x))
}

BIC.tvcm <- function(object, ...) {
    args <- list(...)
    args$object <- object
    args$k <- log(nobs(object))
    args$criterion <- "BIC"
    return(do.call("AIC.tvcm", args))
}

oobrisk.tvcm <- function(object, newdata = NULL, weights = NULL, 
                            fun = NULL, ...) {
  
  if (is.null(fun)) {
    fun <- function(y, mu, wt)
      sum(object$info$family$dev.resids(y, mu, wt))
  }
  
  if (missing(newdata)) stop("require 'newdata'.")
  if (is.null(weights)) weights <- rep(1.0, nrow(newdata))
  yName <- all.vars(object$info$formula$original)[1]
  yMat <- model.matrix(formula(paste("~ -1 + ", yName)), data = newdata)
  if (object$info$family$family == "binomial" && nrow(yMat) > 1L)
    yMat <- yMat[,2L,drop = FALSE]
  mu <- predict(object, newdata, type = "response", ...)
  rval <- fun(yMat, mu, weights)
  
  return(rval)
}


cvfolds <- function(x, type = c("kfold", "subsampling", "bootstrap"),
                    K = ifelse(type == "kfold", 10L, 50L),
                    prob = 0.5) {
  
  ## check input
  if ("bootstrapping" %in% type) type <- "bootstrap"
  type <- match.arg(type)
  
  if (K < 0) stop("'K' must be a positive number.")
  if (K != as.integer(round(K)))
    warning(paste("'K' has been set to ", K, ".", sep = ""))
  K <- as.integer(round(K))
  if (prob <= 0 | prob >= 1)
    stop("'prob' must be within the interval (0, 1).")

  if (inherits(x, "tvcm")) {
      subject <- vcrpart_slot(extract(x, "model"), "subject")
      if (is.null(subject)) subject <- factor(1:nobs(x))
  } else if (is.factor(x)) {
      subject <- x
  } else {
      subject <- factor(1:nobs(x))
  }
  
  subject <- droplevels(subject)
  N <- nlevels(subject)
  
  getSample <- function(x, prob, replace) {
      x <- factor(x)
      xLevs <- levels(x)
      nLevs <- nlevels(x)
      rSample <- sample(x = xLevs,
                        size = ceiling(prob * nLevs),
                        replace = replace)
    return(table(factor(rSample, levels = xLevs)))
  }
  
  
  if (type == "subsampling") {
    
    selected <- replicate(K, getSample(subject, prob, FALSE))
    folds <- 1L * selected[subject, , drop = FALSE]
    dimnames(folds) <- NULL
    
  } else if (type == "kfold") {
    
    if (type == "kfold" & (K > length(subject)) || (K <= 1)) 
      stop("'K' outside allowable range")
    
    selected <- sample(levels(subject), N)
    folds <- matrix(, length(subject), K)
    for (i in 1:K) {
      subs <- seq(ceiling(N / K) * (i - 1) + 1,
                  min(ceiling(N / K) * i, N), 1)
      folds[, i] <- 1.0 * (!subject %in% selected[subs])
    }
    
  } else if (type == "bootstrap") {
    
    selected <- replicate(K, getSample(levels(subject), 1, TRUE))
    folds <- 1L * selected[subject, , drop = FALSE]
    dimnames(folds) <- NULL
  }
  
  return(folds)
}

cvrisk.tvcm <- function(object, folds = cvfolds(object, "kfold", 10),
                           fun = NULL, maxalpha = NULL, maxwidth = NULL,
                           verbose = FALSE, ...) {
  
  type <- list(...)$type
  if (is.null(type)) type <- "risk"
  
  fixed <- list(...)$fixed
  if (is.null(fixed)) fixed <- FALSE
  
  control <- object$info$control
  control$verbose <- FALSE

  if (is.null(maxalpha)) maxalpha <- control$alpha
  if (is.null(maxwidth)) maxwidth <- control$maxwidth
  if (type %in% c("risk", "stabpath")) control$alpha <- maxalpha
  
  ## get data of 'object'
  data <- model.frame(object)
  weights <- weights(object)
  partvar <- colnames(object$data)
  modelvar <- colnames(model.frame(extract(object, "model")))
  subjectName <- vcrpart_slot(extract(object, "model"), "subjectName")
  
  ## check folds
  if (!is.matrix(folds) | (is.matrix(folds) && nrow(folds) != nrow(data)))
    stop(sQuote("folds"), " must be a ", sQuote("matrix"), " with the number of rows equal the number of observations in ", sQuote("object"), ".")
  K <- ncol(folds)
  if (nrow(folds) != nrow(data)) {
    dataCall <- eval.parent(object$info$call$data)
    if (nrow(folds) != nrow(dataCall))
      stop(sQuote("folds"), " has wrong number of rows.")
    folds <- folds[rownames(dataCall) %in% rownames(data), ]
    rm(dataCall)
  }
  
  ## prepare formulas and calls
  ibArgs <- list()
  ibArgs$formula <- formula(object, "original")
  ibArgs$fit <- object$info$fit
  ibArgs$family <- object$info$family
  ibArgs$control <- control
  ibArgs <- append(ibArgs, object$info$dotargs)

  ## number of rows of output matrices
  nr <- switch(type,
               risk = 1,
               stabpath = length(partvar),
               forest = 1)
  
  cv <- vector(mode = "list", length = K)
  
  ## process cross-validation
  nfails <- 0L
  for (i in 1:K) {

    if (verbose) {
      if (i > 1L) cat("\n")
      cat("* evaluating fold", i, "")
    }
    
    ## re-fit the tree with training sample (in-bag)
    ibSubs <- folds[, i]
    ibSubs <- rep(1:length(ibSubs), ibSubs) # rep. bootstrap replications
    ibArgs$data <- data[ibSubs, , drop = FALSE] # learning data
    if (!is.null(subjectName))
      ibArgs$data[, subjectName] <- factor(paste(ibArgs$data[, subjectName], unlist(lapply(folds[, i], function(x) if (x > 0) 1:x)), sep = ".")) # treat replicated subjects as distinct
    ibArgs$weights <- weights[ibSubs]        
    ibTree <- try(do.call("tvcm", ibArgs), silent = TRUE)
    
    ## prepare call for model with test sample (out-of-bag)
    oobSubs <- folds[, i] <= 0
    
    ## get sequence of alpha values
    if (!inherits(ibTree, "try-error")) {
      
      if (type %in% c("risk", "stabpath")) {

        if (fixed) {
          steps <- 1
          alpha <- control$alpha
        } else {
          alpha <- width <- NULL
          steps <- rev(unlist(lapply(splitpath(object), function(x) x$step)))          
          if (control$method == "mob") {
            alpha <- c(extract(object, "p.value", step = steps), 0.0)[-1L]
            if (length(alpha) > 1)
              for (j in 1:(length(alpha) - 1L))
                if (any(alpha[(j + 1):length(alpha)] > alpha[j])) steps <- steps[-j]
            alpha <- rev(alpha)[steps]
          } else {
            width <- steps + 1L
          }
        }
        cv[[i]] <- vector(mode = "list", length = 2)   
        cv[[i]][[1]] <- if (control$method == "mob") alpha else width
        cv[[i]][[2]] <- matrix(, nr, length(cv[[i]][[1]]))

        if (length(steps) > 0L) {
          for (j in 1:length(steps)) {
            
            ibTreePr <- if (!fixed) try(prune(ibTree, maxstep = steps[j])) else ibTree
            
            if (inherits(ibTreePr, "try-error")) {
              
              nfails <- nfails + 1L
              cv[[i]] <- vector(mode = "list", length = 2L)
              warning("computations for fold ", i, " failed. Omit.")
              break;
              
            } else if (type %in% c("risk")) {        
              
              cv[[i]][[2]][, j] <-
                oobrisk(ibTreePr, newdata = data[oobSubs, ,drop = FALSE],
                        weights = weights[oobSubs], fun = fun)        
              
            } else if (type == "stabpath") {
              
              ## get selected variables in current tree
              cv[[i]][[2]][, j] <-
                as.integer(partvar %in% extract(ibTreePr, "selected"))
            }
            
            if (verbose) cat(".")
            
          }
        }
      } else if (type == "forest") {
        
        if (verbose) cat("...")
        ibTreePr <- ibTree
        cv[[i]] <- vector(mode = "list", length = 2)
        cv[[i]][[1]] <- ibTree$node
        cv[[i]][[2]] <- coef(extract(ibTree, "model"))        
      }
    }
    
    if (verbose)
      if (inherits(ibTreePr, "try-error")) cat("failed") else cat(" OK")
  }
  if (type %in% c("risk", "stabpath")) {
    
    ## delete fails
    cv <- cv[!sapply(cv, function(x) sum(sapply(x, length)) == 0)]
    if (length(cv) < 1L) stop("no valid results.")
    
    ## function that evaluates the risk at each alpha of 'grid'
    getVals <- function(x, grid) {
      rval <- matrix(0, nrow(x[[2]]), length(grid))
      for (i in 1:length(grid)) {
        subs <- which(x[[1L]] <= grid[i])
        if (length(subs) > 1) subs <- which(x[[1L]] == max(x[[1L]][subs]))
        if (length(subs) > 0) rval[, i] <- x[[2L]][,subs]
      }
      return(rval)
    }
    
    ## compute results
    grid <- sort(unique(c(unlist(lapply(cv, function(x) x[[1]])))))      
    rval <- list(par = grid)
    
    if (type == "risk") {
      
      ## make a matrix with the 'risk' for each fold
      ## at each alpha in 'grid'
      
      cv <- lapply(cv, function(x) {
        x[[2]][is.nan(x[[2]]) | is.infinite(x[[2]])] <- NA;
        return(x)
      })
      oobRisk <- t(sapply(cv, getVals, grid = grid))
      if (length(grid) == 1) oobRisk <- t(oobRisk)
      
      oobWeights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
      oobWeights[folds > 0] <- 0
      oobRisk <- oobRisk / colSums(oobWeights)
      rownames(oobRisk) <- paste("fold", 1L:length(cv), sep = "")
      colnames(oobRisk) <- 1:length(grid)
      rval <- append(rval, list(risk = oobRisk))
      meanRisk <- colMeans(rval$risk, na.rm = TRUE)
      subs <- which(meanRisk == min(meanRisk))
      if (length(subs) > 1L) subs <- max(subs)
      rval$bestpar <- rval$par[subs]    
      class(rval) <- "cvrisk.tvcm"
      
    } else if (type == "stabpath") {
      
      ## make a matrix with the probability of each variable
      ## to be selected at each alpha in 'grid'
      
      tmp <- lapply(cv, getVals, grid = grid)
      phat <- nValid <- matrix(0, ncol(object$data), length(grid))
      for (i in 1:length(cv)) {
        nValid <- nValid + 1 * !is.na(tmp[[i]])
        tmp[[i]][is.na(tmp[[1]])] <- 0
        phat <- phat + tmp[[i]]
      }
      phat <- phat / nValid
      rownames(phat) <- colnames(object$data)
      colnames(phat) <- 1:length(grid)
      rval <- append(rval, list(phat = phat))
      
      nselected <- rowMeans(sapply(tmp, function(x) colSums(x)))
      rval <- append(rval, list(nselected = nselected))
    }
    
    if (control$method == "mob") {
      rval$maxpar <- maxalpha
      rval$parname <- "alpha"
    } else {
      rval$maxpar <- maxwidth
      rval$parname <- "maxwidth"
    }   
    rval$call <- deparseCall(getCall(object))
    
  } else if (type == "forest") {
    
    rval <- vector(mode = "list", length = 3)
    names(rval) <- c("node", "coefficients", "error")
    rval$error$which <-
      which(sapply(cv, function(x) inherits(x, "try-error")))
    rval$error$message <- unlist(cv[rval$error$which])
    cv[rval$error$which] <- NULL
    rval$node <- lapply(cv, function(x) x[[1]])
    rval$coefficients <- lapply(cv, function(x) x[[2]])
  }
  
  if (verbose) cat("\n")
  
  return(rval)
}

print.cvrisk.tvcm <- function(x, ...) {
  cat("Cross-validated observation-mean risk", "\n\n")
  cat("Call: ", x$call, "\n\n")
  rval <- colMeans(x$risk)
  names(rval) <- format(x$par, digits = 2)
  print(rval)
  cat(paste("\nOptimal ",x$parname, ": ", format(x$bestpar, digits = 2),
            "\n", sep = ""))
  return(invisible(x))
}

plot.cvrisk.tvcm <- function(x, ...) {
  
  dotList <- list(...)
  if (x$parname == "alpha") {
      defArgs <- list(type = "S", xlab = expression(alpha),
                      col = "grey", lty = 1,
                      ylab = expression(paste("observation-mean risk(",
                          alpha, ")", sep = "")))
  } else {
      defArgs <- list(type = "S", xlab = "maxwidth",
                      col = "grey", lty = 1,
                      ylab = "observation-mean risk(maxwidth)")
  }

  yy <- t(x$risk)
  xx <- matrix(x$par, nrow(yy), ncol(yy))

  ## delete zeros if 'log = "x"' is called
  if ("log" %in% names(dotList) && dotList$log %in% c("x", "xy")) {
    xx <- xx[-1, ]
    yy <- yy[-1, ]
  }   

  ## set plot arguments
  plotArgs <- appendDefArgs(list(x = xx, y = yy), list(...))
  plotArgs <- appendDefArgs(plotArgs, defArgs)
    
  ## plot 
  do.call("matplot", plotArgs)
      
  ## mean curve for 'type = "risk"'
  meanRisk <- colMeans(x$risk, na.rm = TRUE)
  points(x = x$par, meanRisk, type = "S")
    
  ## best alpha parameter
  bestPar <- x$bestpar
  bestRisk <- min(meanRisk)

  segments(bestPar, par()$usr[3], bestPar, bestRisk, lty = 2)
  axis(1, bestPar, format(bestPar, digits = 2),
       line = 1, tick = FALSE)   
}

stabpath.tvcm <- function(object, q, maxalpha = NULL, maxwidth = NULL, ...) {

  object$info$control$nselect <- q
  cv <- cvrisk(object, type = "stabpath", ...)
  if (max(cv$nselected) < q)
    warning(sQuote("alpha"), " too small, the average",
            "number of selected partitioning variables is ",
            max(cv$nselected))  
  mm <- apply(cv$phat, 1L, max)
  rval <- append(cv, list(max = mm, q = q, maxalpha = maxalpha,
                          maxwidth = maxwidth))
  class(rval) <- "stabpath.tvcm"
  return(rval)
}

print.stabpath.tvcm <- function(x, ...) {

  cat("Stability path\n\n")
  cat("call:", x$call, "\n")
  cat("q:", x$q, "\n")
  cat("Selection probabilities:\n\n")
  print(x$max)
  return(invisible(x))
}

plot.stabpath.tvcm <- function(x, ...) {

  dotList <- list(...)
  
  col <- rep(rainbowPalette, length.out = length(x$max))
  lty <- rep(1:4, length.out = length(x$max))

  if (x$parname == "alpha") {
      defArgs <- list(type = "S", xlab = expression(alpha),
                      col = "grey", lty = 1,
                      ylab = expression(paste("observation-mean risk(",
                          alpha, ")", sep = "")))
  } else {
      defArgs <- list(type = "S", xlab = "maxwidth",
                      col = "grey", lty = 1,
                      ylab = "risk(maxwidth)")
  }
  
  yy <- t(x[["phat"]])
  xx <- matrix(x[["par"]], nrow(yy), ncol(yy))
  
  ## delete zeros if 'log = "x"' is called
  if ("log" %in% names(dotList) && dotList$log %in% c("x", "xy")) {
    xx <- xx[-1, ]
    yy <- yy[-1, ]
  } 
  
  ## set plot arguments
  plotArgs <- appendDefArgs(list(x = xx, y = yy), list(...))
  plotArgs <- appendDefArgs(plotArgs, defArgs)
  
  ## plot 
  do.call("matplot", plotArgs)
  max <- apply(x[["phat"]], 1L, max)
  axis(4, max, rownames(x[["phat"]]), las = 1)
}
