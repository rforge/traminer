## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2013-12-06
##
## Description:
## Cross-validation, stability paths and random forests for
## 'tvcolmm' objects.
##
## Contents:
## tvcolmm_folds:     create cross-validation folds
## cv.vcolmm:        cross-validation for 'tvcolmm' objects
## print.cv.tvcolmm: print for 'cv.tvcolmm' objects
## plot.cv.tvcolmm:  plot fot 'cv.tvcolmm' objects
## stabpath.tvcolmm: stability selection for 'tvcolmm' objects
## print.stabpath.tvcolmm: print for 'stabpath.tvcolmm' objects
## plot.stabpath.tvcolmm: plot for 'stabpath.tvcolmm' objects
##
## Last modifications:
##
## To do:
## - create utility function of contents of 'cv.tvcolmm' (by now
##   'stabpath.tvcolmm' and 'forest.tvcolmm' call internally
##   'cv.tvcolmm' which may lead to some confusion
## --------------------------------------------------------- #

tvcolmm_folds <- function(x, type = c("subsampling", "kfold", "bootstrap"),
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

  if (class(x)[1] == "tvcolmm") {
      subject <- extract(x, "model")@subject
  } else if (is.factor(x)) {
      subject <- x
  } else {
      stop("'x' must be a 'tvcolmm' object or a 'factor'.")
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

cv.tvcolmm <- function(object, folds = tvcolmm_folds(object),
                       alpha.max = 0.5, verbose = FALSE, ...) {

  type <- list(...)$type
  if (is.null(type)) type <- "loss"

  fixed <- list(...)$fixed
  if (is.null(fixed)) fixed <- FALSE
  
  control <- object$info$control
  control$verbose <- FALSE
  if (type %in% c("loss", "stabpath")) control$alpha <- alpha.max
  
  ## get data of 'object'
  data <- model.frame(object)
  weights <- weights(extract(object, "model"))
  partvar <- colnames(object$data)
  modelvar <- colnames(model.frame(extract(object, "model")))
  subjectName <- extract(object, "model")@subjectName
  yname <- all.vars(object$info$formula$root)[1]
  yMat <- model.matrix(formula(paste("~ -1 + ", yname)), data = data)
  
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
  ibArgs$formula <- formula(paste(deparse(object$info$formula$root, width.cutoff = 500L), "|", paste(partvar, collapse = " + ")))    
  ibArgs$control <- control
  ibArgs$vi <- object$info$vi
  ibArgs$linear <- object$info$linear
  ibArgs <- append(ibArgs, object$info$dotargs)

  ## number of rows of output matrices
  nr <- switch(type,
               loss = 1,
               stabpath = length(partvar),
               forest = 1)
  
  cv <- vector(mode = "list", length = K)
  
  ## process cross-validation
  nfails <- 0L
  for (i in 1:K) {
    
    if (verbose) cat("\n* evaluating fold", i, "")
    
    ## re-fit the tree with training sample (in-bag)
    ibSubs <- folds[, i]
    ibSubs <- rep(1:length(ibSubs), ibSubs) # rep. bootstrap replications
    ibArgs$data <- data[ibSubs, , drop = FALSE] # learning data
    ibArgs$data[, subjectName] <- factor(paste(ibArgs$data[, subjectName], unlist(lapply(folds[, i], function(x) if (x > 0) 1:x)), sep = ".")) # treat replicated subjects as distinct
    ibArgs$weights <- weights[ibSubs]        
    ibTree <- try(do.call("tvcolmm", ibArgs), silent = TRUE)
    
    ## prepare call for model with test sample (out-of-bag)
    oobSubs <- folds[, i] <= 0
    
    ## get sequence of alpha values
    if (!inherits(ibTree, "try-error")) {
      
      nodes <- setdiff(nodeids(ibTree, terminal = FALSE),
                       nodeids(ibTree, terminal = TRUE))
      
      if (type %in% c("loss", "stabpath")) {

          if (fixed) {
              alpha <- control$alpha
          } else {
              alpha <-
                  c(0, if (depth(ibTree) > 0)
                    sort(extract(ibTree, "pval", nodes)))
          }
          cv[[i]] <- vector(mode = "list", length = 2)   
          cv[[i]][[1]] <- alpha
          cv[[i]][[2]] <- matrix(, nr, length(alpha))
          
      } else if (type == "forest") {

        if (verbose) cat("...")
        cv[[i]] <- vector(mode = "list", length = 2)
        cv[[i]][[1]] <- ibTree$node
        cv[[i]][[2]] <- coef(extract(ibTree, "model"))
        alpha <- c()
      }

    } else {
      
      alpha <- c() # avoids the following loop to start
      
      if (type == "forest")
        cv[[i]] <- ibTree
    }

    if (length(alpha) > 0) {
      for (j in 1:length(alpha)) {
        
        ibTree <- try(prune(ibTree, alpha = rev(alpha)[j]))
        
        if (inherits(ibTree, "try-error")) {
          
          nfails <- nfails + 1L
          cv[[i]] <- vector(mode = "list", length = 2)
          warning("computations for fold ", i, " failed. Omit.")
          break;
          
        } else if (type == "loss") {        
          
          pred <- predict(ibTree, newdata = data[oobSubs, ,drop = FALSE])
          cv[[i]][[2]][, length(alpha) - j + 1L] <-
            sum(-log(pred[yMat[oobSubs,] > 0]))
          
        } else if (type == "stabpath") {
          
          ## get selected variables in current tree
          cv[[i]][[2]][, length(alpha) - j + 1L] <- as.integer(partvar %in% extract(ibTree, "selected"))
          
        }
        
        if (verbose) cat(".")
      }
    }
    if (verbose)
      if (inherits(ibTree, "try-error")) cat("failed") else cat(" OK")
  }
  
  if (type %in% c("loss", "stabpath")) {

    ## delete fails
    cv <- cv[!sapply(cv, function(x) sum(sapply(x, length)) == 0)]
    if (length(cv) < 1L) stop("no valid results.")
    
    ## function that evaluates the outcome at each alpha in grid
    getVals <- function(x, grid) {
      rval <- matrix(0, nrow(x[[2L]]), length(grid))
      for (i in 1:length(grid)) {
        subs <- max(which(x[[1]] <= grid[i]))
        rval[, i] <- x[[2L]][,subs]
      }
      return(rval)
    }
  
    ## compute results
    
    grid <-
        sort(unique(c(unlist(lapply(cv, function(x) x[[1]])), alpha.max)))
    rval <- list(alpha = grid)
    
    if (type == "loss") {
      
      ## make a matrix with the 'loss' for each fold
      ## at each alpha in 'grid'
      
      cv <- lapply(cv, function(x) {
        x[[2]][is.nan(x[[2]]) | is.infinite(x[[2]])] <- NA;
        return(x)
      })
      oobLoss <- t(sapply(cv, getVals, grid = grid))
      if (length(grid) == 1) oobLoss <- t(oobLoss)
          
      oobWeights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
      oobWeights[folds > 0] <- 0
      oobLoss <- oobLoss / colSums(oobWeights)
      rownames(oobLoss) <- paste("fold", 1L:length(cv), sep = "")
      colnames(oobLoss) <- 1:length(grid)
      rval <- append(rval, list(loss = oobLoss))   
      
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
  
  ## compute optimal alpha
  if (type == "loss") {
    
    meanLoss <- colMeans(rval[["loss"]], na.rm = TRUE)
    subs <- which(meanLoss == min(meanLoss))
    if (length(subs) > 1L) subs <- max(subs)
    rval$alpha.best <- rval[["alpha"]][subs]
    rval$call <- deparse(getCall(object), nlines = 1, width.cutoff = 45)[1]
    class(rval) <- "cv.tvcolmm"
  }
  
  return(rval)
}

print.cv.tvcolmm <- function(x, ...) {

    cat("\n\t Cross-validated loss", "\n\t", x$call, "\n\n")
    rval <- colMeans(x$loss)
    names(rval) <- format(x$alpha, digits = 2)
    print(rval)
    cat("\n\t Optimal alpha:", format(x$alpha.best, digits = 2), "\n")
    return(invisible(x))
}

plot.cv.tvcolmm <- function(x, ...) {
  
  dotList <- list(...)
  defArgs <- list(type = "S", xlab = expression(alpha),
                  col = "grey", lty = 1,
                  ylab = expression(paste("loss(", alpha, ")")))

  yy <- t(x[["loss"]])
  xx <- matrix(x[["alpha"]], nrow(yy), ncol(yy))

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
      
  ## mean curve for 'type = "loss"'
  meanLoss <- colMeans(x[["loss"]], na.rm = TRUE)
  points(x = x[["alpha"]], meanLoss, type = "S")
    
  ## best alpha parameter
  bestAlpha <- x$alpha.best
  bestLoss <- min(meanLoss)

  segments(bestAlpha, par()$usr[3], bestAlpha, bestLoss, lty = 2)
  axis(1, bestAlpha, format(bestAlpha, digits = 2),
       line = 1, tick = FALSE)   
}

stabpath.tvcolmm <- function(object, q = 2L, ...) {

  object$info$control$nselect <- q
  
  cv <- cv.tvcolmm(object, type = "stabpath", ...)
  
  if (max(cv$nselected) < q)
    warning(sQuote("alpha"), " too small, the average number of selected partitioning variables is ", max(cv$nselected)) 
  
  mm <- apply(cv$phat, 1L, max)
  rval <- append(cv, list(max = mm, q = q))
  class(rval) <- "stabpath.tvcolmm"
  return(rval)
}

print.stabpath.tvcolmm <- function(x, ...) {

  cat("\tStability Path\n")
  if (length(x$selected) > 0) {
    cat("\nSelected partitioning variables:\n")
    print(x$selected)
  } else {
    cat("\nNo variables selected\n")
  }
  cat("\nSelection probabilities:\n")
  print(x$max[x$max > 0])
  cat("q: ", x$q, "\n\n")
  return(invisible(x))
}

plot.stabpath.tvcolmm <- function(x, ...) {

  dotList <- list(...)
  
  col <- rep(rainbowPalette, length.out = length(x$max))
  lty <- rep(1:4, length.out = length(x$max))
  
  defArgs <- list(type = "S", xlab = expression(alpha),
                  col = col, lty = lty, ylim = c(0, 1),
                  ylab = expression(Pi(alpha)))
  
  yy <- t(x[["phat"]])
  xx <- matrix(x[["alpha"]], nrow(yy), ncol(yy))
  
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
