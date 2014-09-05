## --------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
##' Date:            2014-09-02
##'
##' Description:
##' Function for model selection and assessment for 'tvcm' objects.
##'
##' Contents:
##' oobloss.tvcm:        computes the out.of-bag loss
##' folds:               parameters for cross-validation folds
##' tvcm_folds:          create cross-validation folds
##' cvloss.tvcm:         cross-validation for 'tvcm' objects
##' print.cvloss.tvcm:   print for 'cv.tvcm' objects
##' plot.cvloss.tvcm:    plot fot 'cv.tvcm' objects
##'
##' Last modifications:
##' 2014-09-02: - modifications on 'tvcm_get_node'. The former
##'               implementation was a time-killer an therefore
##'               there is a new argument 'formList' and a
##'               auxiliary function 'tvcm_get_fitted' was added
##' 2014-08-30: - deleted 'nsplit' extension of in 'cvloss'. Reason:
##'               cross-validation for number of split would
##'               require a different pruning procedure
##'             - small justifications for the plot
##' 2014-08-06: - substituted 'cvfolds' by 'folds', which defines
##'               a list of parameters for the new function
##'               'tvcm_folds' that creates the cross-validation
##'               matrix
##'               matrix
##' 2014-07-29: - defined 'cvfolds' as method for 'tvcm' objects
##'             - added new argument 'weights' to 'cvfolds' to
##'               allow for models where the weights represent
##'               counts
##' 2014-07-22: - removed AIC and BIC methods since they do
##'               not apply to the 'tvcm' framework
##' 2014-07-17: - change the column names of the cross-validation
##'               output matrix
##'             - improve the desciption
##' 2014-07-08: - remove the 'sub' argument of AIC.tvcm, BIC.tvcm
##'             - remove additional methods for AIC.tvcm and BIC.tvcm
##'             - remove stabsel.tvcm and and corresponding methods
##'             - cvloss: set 'cp' as the only tuning parameter
##'             - cvloss: add 'direction' as new parameter 
## --------------------------------------------------------- #


oobloss.tvcm <- function(object, newdata = NULL, weights = NULL, 
                         fun = NULL, ...) {
  
  if (is.null(fun)) {
    fun <- function(y, mu, wt)
      sum(object$info$family$dev.resids(y, mu, wt))
  }
  
  if (missing(newdata)) stop("require 'newdata'.")
  if (is.null(weights)) weights <- rep(1.0, nrow(newdata))
  yName <- all.vars(object$info$formula$original)[1]
  yMat <- model.matrix(formula(paste("~ -1 + ", yName)), data = newdata)
  if (object$info$family$family == "binomial" && ncol(yMat) > 1L)
    yMat <- yMat[,2L,drop = FALSE]
  mu <- suppressWarnings(predict(object, newdata, type = "response", ...))
  rval <- fun(yMat, mu, weights)
  
  return(rval)
}


folds_control<- function(type = c("kfold", "subsampling", "bootstrap"),
                         K = ifelse(type == "kfold", 5, 30),
                         prob = 0.5, weights = c("case", "freq"),
                         seed = NULL) {
  if ("bootstrapping" %in% type) type <- "bootstrap"
  type <- match.arg(type)
  stopifnot(is.numeric(K) && length(K) == 1L)
  if (round(K) < 1L) stop("'K' must be a positive number.")
  if (type == "kfold" && K < 2L)
    stop("'K' must be larger than 1 for 'type = 'kfold''")
  if (K != as.integer(round(K)))
    warning(paste("'K' has been set to ", K, ".", sep = ""))
  K <- as.integer(round(K))
  stopifnot(is.numeric(prob))
  if (prob <= 0 | prob >= 1)
    stop("'prob' must be within the interval (0, 1).")
  K <- round(K)
  stopifnot(is.numeric(prob) && length(prob) == 1L)
  weights <- match.arg(weights)
  return(structure(list(type = type,
                        K = K,
                        prob = prob,
                        weights = weights,
                        seed = seed),
                   class = "folds"))
}


## --------------------------------------------------------- #
##' Creates a cross-validation matrix
##'
##' @param object  an object of class \code{tvcm}
##' @param args    a list of arguments as produced by
##'    \code{\link{folds}}.
##'
##' @return A matrix.
## --------------------------------------------------------- #

tvcm_folds <- function(object, control) {

  stopifnot(inherits(control, "folds"))
  
  ## detach input
  type <- control$type
  K <- control$K
  prob <- control$prob
  weights <- control$weights
  seed <- control$seed

  subject <- extract(object, "model")$subject
  if (!is.null(subject) && weights == "freq")
    stop("option 'weights = 'freq'' is not available for 2-stage data")

  freq <- switch(weights,
                 case = rep(1, nobs(extract(object, "model"))),
                 freq = weights(extract(object, "model")))
  if (weights == "freq" && any(!freq == round(freq)))
    stop("some of the weights are not integers.")

  if (is.null(subject)) subject <- factor(rep(1:length(freq), freq))
  N <- nlevels(subject)
    
  getSample <- function(subject, freq, prob, replace, weights) {
    levs <- if (weights == "case") 1:nlevels(subject) else freq = 1:length(subject)
    rSample <- sample(x = levs,
                      size = ceiling(prob * length(levs)),
                      replace = replace)
    rSample <- table(factor(rSample, levels = levs))
    if (weights == "case") rSample <- rSample[subject]
    if (weights == "freq") rSample <- tapply(rSample, subject, sum)
    return(rSample)
  }

  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
  if (!is.null(seed)) set.seed(seed)
  RNGstate <- .Random.seed
  
  if (type == "subsampling") {
    
    folds <- replicate(K, getSample(subject, freq, prob, FALSE, weights))
    
  } else if (type == "kfold") {
    
    if (type == "kfold" & (K > length(subject)) || (K <= 1)) 
      stop("'K' outside allowable range")

    levs <- if (weights == "case") 1:nlevels(subject) else 1:length(subject)
    selected <- sample(levs, length(levs))
    split <- c(0, quantile(levs, (1:(K - 1)) / K, type = 1L), max(levs))
    folds <- matrix(, length(levs), K)
    for (i in 1:K) {
      subs <- levs > split[i] & levs <= split[i + 1L]
      folds[, i] <- 1.0 * (!levs %in% selected[subs])
    }
    if (weights == "case") folds <- folds[subject, ]
    if (weights == "freq") folds <- apply(folds, 2, tapply, subject, sum)
    
  } else if (type == "bootstrap") {
    
    folds <- replicate(K, getSample(subject, freq, prob = 1, TRUE, weights))
  }

  colnames(folds) <- 1:K
  rownames(folds) <- rownames(model.frame(extract(object, "model")))

  assign(".Random.seed", oldSeed, envir=globalenv())
  attr(folds, "type") <- type
  attr(folds, "value") <- ifelse(weights == "freq", "weights", "freq")
  attr(folds, "seed") <- RNGstate
  return(folds)
}


cvloss.tvcm <- function(object, folds = folds_control(),
                        fun = NULL, dfpar = 2, direction = c("backward", "forward"),
                        papply = mclapply, verbose = FALSE, ...) {

  
  stopifnot(inherits(folds, "folds"))
  stopifnot(is.numeric(dfpar) && length(dfpar) == 1L)
  type <- list(...)$type
  direction <- match.arg(direction)
  if (is.null(type)) type <- "loss"
  folds <- tvcm_folds(object, folds)
  
  papplyArgs <- list(...)[names(list(...)) %in% names(formals(papply))]
  keeploss <- list(...)$keeploss
  if (is.null(keeploss)) keeploss <- formals(prune.tvcm)$keeploss
  
  control <- object$info$control
  control$verbose <- FALSE
  control$center <- FALSE
      
  ## get data of 'object'
  data <- model.frame(object)
  
  weights <- weights(object)
  partvar <- colnames(object$data)
  modelvar <- colnames(model.frame(extract(object, "model")))
  subjectName <- extract(object, "model")$subjectName
  
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
  ibCall <- call(name = "tvcm",
                 formula = quote(object$info$formula$original),
                 data = quote(ibData),
                 fit = quote(object$info$fit),
                 family = quote(object$info$family),
                 control = quote(control))
  for (arg in names(object$info$dotargs)) ibCall[[arg]] <- object$info$dotargs[[arg]]
    
  ## process cross-validation

  cvFun <- function(i) {

    cv <- vector(mode = "list", length = 2)
    
    if (verbose) {
      if (!identical(papply, mclapply) && i > 1L) cat("\n")
      if (identical(papply, mclapply)) cat("[", i, "]") else cat("* fold", i, "...")
    }
    
    ## extract subset
    if (attr(folds, "value") == "weights") {
      ibSubs <- folds[, i] > 0.0
      ibData <- data[ibSubs,, drop = FALSE]
      ibCall$weights <- folds[ibSubs, i]
      oobSubs <- rep(TRUE, nrow(folds))
      oobWeights <- weights - folds[, i]
    } else {
      ibSubs <- folds[, i]
      ibSubs <- rep(1:length(ibSubs), ibSubs) # rep. bootstrap replications
      ibData <- data[ibSubs, , drop = FALSE] # learning data
      if (!is.null(subjectName))
        ibData[, subjectName] <- factor(paste(ibData[, subjectName], unlist(lapply(folds[, i], function(x) if (x > 0) 1:x)), sep = ".")) # treat replicated subjects as distinct
      ibCall$weights <- weights[ibSubs]
      oobSubs <- folds[, i] <= 0
      oobWeights <- weights[oobSubs]
    }

    ## re-fit the tree with training sample (in-bag)
    ibTree <- try(eval(ibCall), silent = TRUE)
    dfsplit <- switch(direction, backward = 0.0, forward = 0.0)
    
    if (!inherits(ibTree, "try-error")) {
      
      if (type == "loss") {
        run <- 1L
        
        while (run > 0L) {
          ibTree <- try(prune(ibTree, dfsplit, dfpar, direction,
                              papply = lapply, keeploss = keeploss), TRUE)
          if (!inherits(ibTree, "try-error")) {
            ## save the out-of-bag loss and the current tuning parameter
            oobLoss <- oobloss(ibTree, newdata = data[oobSubs,,drop = FALSE],
                               weights = oobWeights, fun = fun)
            cv[[1L]] <-
              cbind(cv[[1L]], c(dfsplit))
            cv[[2L]] <-
              cbind(cv[[2L]], c(control$lossfun(ibTree), oobLoss))

            ## set a new and stronger tuning parameter
            tab <- ibTree$info$prunepath[[length(ibTree$info$prunepath)]]$tab
            if (nrow(tab) > 1L) dfsplit <- min(tab$dfsplit[-1L])
            if (nrow(tab) == 1L) run <- 0L
            
          } else {
            cv <- NULL
            run <- 0L

          }
        }
  
      } else if (type == "forest") {        
        if (verbose && !identical(papply, mclapply)) cat("...")
        ibTreePr <- ibTree
        cv[[1]] <- ibTree$info$node
        cv[[2]] <- coef(extract(ibTree, "model"))
        
      }

    } else {
      cv <- NULL

    }
    if (verbose) {
      if (identical(papply, mclapply)) {
        if (is.null(cv)) cat("failed")
      } else {
        if (is.null(cv)) cat("failed") else cat(" OK")
      }
    }
    ## return output
    return(cv)
  }

  cv <- do.call(papply, append(list(X = 1:ncol(folds), FUN = cvFun), papplyArgs))
  
  if (type %in% c("loss")) {
    
    ## delete fails
    fails <- sapply(cv, function(x) sum(sapply(x, length)) == 0)
    if (all(fails)) stop("no valid results.")
    cv <- cv[!fails]
    
    if (type == "loss") {

      ## function that evaluates the loss at each alpha of 'grid'
      getVals <- function(x, grid, rowsG = 1L, rowsL = 1L) {
        rval <- matrix(NA, length(row), length(grid))
        for (i in 1:length(grid)) {
          cols <- which(x[[1L]][rowsG,] <= grid[i])
          if (length(cols) > 1)
            cols <- which(x[[1L]][rowsG, ] == max(x[[1L]][rowsG, cols]))
          if (length(cols) > 0)
            rval[, i] <- x[[2L]][rowsL, cols]
        }
        return(rval)
      }
    
      ## compute results
      grid <- sort(unique(c(unlist(lapply(cv, function(x) x[[1]][1, ])))))
      rval <- list(grid = grid)
      
      ## column names
      if (length(grid) > 1L) {
        cn <- c(paste("<=", round(grid[2L], 2)), paste(">", round(grid[-1L], 2)))
      } else {
        cn <- "0"
      }
     
      ## make a matrix with the 'loss' for each fold
      ## at each dfsplit in 'grid'
      
      cv <- lapply(cv, function(x) {
        x[[2]][is.nan(x[[2]]) | is.infinite(x[[2]])] <- NA;
        return(x)
      })
      
      ## ib-loss
      ibLoss <- t(sapply(cv, getVals, grid = grid, rowsG = 1L, rowsL = 1L))
      if (length(grid) == 1L) ibLoss <- t(ibLoss)
     
      if (attr(folds, "value") == "weights") {
        ibWeights <- folds
      } else {
        ibWeights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
        ibWeights[folds <= 0] <- 0
      }
      ibLoss <- ibLoss / colSums(ibWeights)[!fails]
      rownames(ibLoss) <- paste("fold", 1L:length(cv))
      colnames(ibLoss) <- cn
      
      ## oob-loss
      oobLoss <- t(sapply(cv, getVals, grid = grid, rowsG = 1L, rowsL = 2L))
      if (length(grid) == 1L) oobLoss <- t(oobLoss)
       
      oobWeights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
      if (attr(folds, "value") == "weights") {
        oobWeights <- oobWeights - folds
      } else {
        oobWeights[folds > 0] <- 0
      }
      oobLoss <- oobLoss / colSums(oobWeights)[!fails]
      rownames(oobLoss) <- paste("fold", 1L:length(cv))
      colnames(oobLoss) <- cn
      
      rval <- append(rval, list(ibloss = ibLoss, oobloss = oobLoss))
      meanLoss <- colMeans(rval$oobloss, na.rm = TRUE)

      ## dfsplit with minimal loss
      minSubs <- which(meanLoss == min(meanLoss))
      if (length(minSubs) > 1L) minSubs <- max(minSubs)
      rval$dfsplit.hat <- max(0, mean(c(rval$grid, Inf)[minSubs:(minSubs + 1)]))

      rval$folds <- folds
      class(rval) <- "cvloss.tvcm"     
    }

    rval$direction <- direction
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
    rval$folds <- folds
  }
  
  if (verbose) cat("\n")
  
  return(rval)
}


print.cvloss.tvcm <- function(x, ...) {
  cat(ifelse(x$direction == "backward", "Backwards", "Forward"),
      "cross-validated average loss", "\n\n")
  cat("Call: ", x$call, "\n\n")
  rval <- colMeans(x$oobloss, na.rm = TRUE)
  if (length(rval) > 10L)
    rval <- rval[seq(1L, length(rval), length.out = 10)]
  print(rval)
  cat("\n")
  cat("dfsplit with minimal loss:", format(x$dfsplit.hat, digits = 3), "\n")
  return(invisible(x))
}


plot.cvloss.tvcm <- function(x, legend = TRUE, details = TRUE, ...) {
  
  xlab <- "dfsplit"
  ylab <- "average loss(dfsplit)"
  type <- "s"
  lpos <- "topleft"
  lsubs <- if (details) 1L:3L else 1L
  if (x$dfsplit.hat < Inf) lsubs <- c(lsubs, 4L)
  lcol <- c("black", "grey80", "black", "black")
  llty <- c(1, 1, 3, 2)

  dotList <- list(...)
  defArgs <- list(type = type, xlab = xlab, ylab = ylab,
                  col = lcol[1], lty = llty[1])
  
  xx <- matrix(x$grid, length(x$grid), nrow(x$oobloss))  
  xx <- rbind(xx, 1.1 * xx[nrow(xx),])

  yy2 <- t(cbind(x$oobloss, x$oobloss[,ncol(x$oobloss)]))
  yy1 <- rowMeans(yy2, na.rm = TRUE)
  yy3 <- rowMeans(t(cbind(x$ibloss, x$ibloss[,ncol(x$ibloss)])), na.rm = TRUE)
  
  defArgs$ylim <- range(if (details) c(yy2, yy3) else yy1, na.rm = TRUE)
  
  ## set plot arguments
  plotArgs <- appendDefArgs(list(x = xx[, 1], y = yy1), list(...))
  plotArgs <- appendDefArgs(plotArgs, defArgs)
  llty[1L] <- plotArgs$lty
  lcol[1L] <- plotArgs$col
  
  ## plot mean validated prediction error
  do.call("plot", plotArgs)
  
  ## plot details
  if (details) {
    matplot(x = xx, yy2, type = type, lty = llty[2L], col = lcol[2L], add = TRUE)
    points(x = xx[,1], yy3, type = type, lty = llty[3L], col = lcol[3L])
  }
  
  if (legend) {
    ltext <- c("average oob estimate",
               "foldwise oob estimates",
               "in-sample estimate",
               "dfsplit with minimal oob-loss")
    legend(lpos, ltext[lsubs], col = lcol[lsubs], lty = llty[lsubs])
  }
  
  if (x$dfsplit.hat < Inf) {
    subsMin <- max(which(x$grid <= x$dfsplit.hat))
    minLoss <- colMeans(x$oobloss, na.rm = TRUE)[subsMin]
    segments(x$dfsplit.hat, par()$usr[3], x$dfsplit.hat, minLoss, lty = 2)
    axis(1, x$dfsplit.hat, format(x$dfsplit.hat, digits = 3),
         line = 1, tick = FALSE)
  } 
}
