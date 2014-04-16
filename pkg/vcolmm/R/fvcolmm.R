## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-04-13
##
## Description:
## Random forests based on the tvcolmm algorithm.
##
## Contents:
## fvcolmm:         random forest estimator for 'tvcolmm' objects
## predict.fvcolmm: prediction for 'forest.tvcolmm' objects
## print.fvcolmm:   print method for 'forest.tvcolmm' objects
##
## To do:
## --------------------------------------------------------- #


fvcolmm <- function(..., folds = NULL, control = fvcolmm_control(), add = TRUE) {
  
  mc <- match.call()
  
  if (class(eval(mc[[2]]))[1] == "fvcolmm") {
    
    ## take the old structure as prototyp tree
    object <- list(...)[[1]]
    
    ## update control argument
    if (!is.null(mc$control))
      for (i in names(mc$control))
        object$info$control[[i]] <- mc$control[[i]]
    verbose <- object$info$control$verbose
    object$info$control$verbose <- FALSE
    
  } else {
    
    ## modify control parameters temporarily for a first tree
    maxdepth <- control$maxdepth
    control$maxdepth <- 0L
    verbose <- control$verbose
    control$verbose <- FALSE
    
    ## fit a prototyp tree
    if (verbose) cat("\n* fitting an initial tree ... ")
    initargs <- append(list(...), list(control = control))
    object <- do.call("tvcolmm", args = initargs)
    if (verbose) cat("OK\n")
    
    ## reset the depth parameter and set the verbose parameter
    object$info$control$maxdepth <- maxdepth
  }
  
  if (is.null(folds)) folds <- tvcolmm_folds(object)
  
  ## compute trees for subsamples
  args <- list(object = object,
               folds = folds,
               type = "forest",
               verbose = verbose)
  cv <- do.call("cv.tvcolmm", args = args)
  fails <- cv$error$which  
  
  ## add new information to info slot
  if (length(fails) > 0)
    folds <- folds[, -fails, drop = FALSE]

  if (!add) object$info$node <- object$info$coefficients <- NULL
  
  object$info$node <- append(object$info$node, cv$node)
  object$info$coefficients <-
    append(cv$coefficients, object$info$coefficient)
  object$info$folds <- cbind(object$info$folds, folds)
  object$info$error <- cv$error
  object$info$control$verbose <- verbose
  
  ## modifiy class attribute to allow methods
  class(object) <- append("fvcolmm", class(object))
  return(object)
}

fvcolmm_control <- function(alpha = 1.0, bonferroni = TRUE,
                            minsplit = 50L, minbucket = 25L,
                            maxdepth = Inf, maxwidth = 10L,
                            mtry = 5L, nselect = Inf,
                            estfun = list(), sctest = TRUE, 
                            maxevalsplit = 20, lossfun = neglogLik,
                            fast = 0L, verbose = FALSE,...) {

  return(tvcolmm_control(
           alpha = alpha, bonferroni = bonferroni,
           minsplit = minsplit, minbucket = minbucket,
           maxdepth = maxdepth, maxwidth = maxwidth,
           mtry = mtry, nselect = nselect,
           estfun = estfun, sctest = sctest, 
           maxevalsplit = maxevalsplit, lossfun = lossfun,
           fast = fast, verbose = verbose, ...))
}

logLik.fvcolmm <- function(object, ...) {
  yMat <- model.matrix(~ -1 + object$fitted[,"(response)"])
  pred <- predict(object, oob = TRUE, ...)
  weights <- weights(extract(object, "model"))
  return(sum(log(rowSums(yMat * pred) * weights)))
}

neglogLik.fvcolmm <- function(object, ...) {
  return(-logLik(object, ...))
}

predict.fvcolmm <- function(object, newdata = NULL,
                            type = c("prob", "class", "link",
                              "terms"), ranef = FALSE, 
                            verbose = FALSE, ...) {
  
  type <- match.arg(type)
  tp <- if (type == "terms") type else "link"

  if (!"oob" %in% list(...)) oob <- FALSE
  if (!"which" %in% list(...)) which <- 1:length(object$info$node)
  
  ## modifiy class to apply 'tvcolmm' methods
  class(object)[1] <- "tvcolmm"
  if (is.null(newdata))  {
    newdata <- model.frame(object)
    folds <- object$info$folds
    weights <- object$fitted[,"(weights)"]
  } else {
    folds <- matrix(1, nrow(newdata), ncol(object$info$folds))
    if ("weights" %in% list(...)) {
      weights <- list(...)$weights
    } else {
      weights <- rep(1, nrow(newdata))
    }
  }

  if (!all(which %in% 1:length(object$info$node)))
    stop("'which' is missspecified.")
  folds <- folds[, which, drop = FALSE]
  object$info$node <- object$info$node[which]
  object$info$coefficients <- object$info$coefficients[which]
  
  K <- length(object$info$node)
  
  ## names of effect-modifiers and predictors
  partvar <- colnames(object$data)
  modelvar <- colnames(model.frame(extract(object, "model")))
  
  ## prepare arguments to build the model
  args <- list()  
  args$data <- model.frame(object$info$model)
  args$doFit <- FALSE
  args <- appendDefArgs(args, object$info$dotargs)
  
  if (verbose) cat("\n* predicting ")
  
  ## extract each tree and aggregate predicted coefficients
  for (i in 1:K) {

    if (verbose) cat(".")
    
    ## build model
    object$node <- object$info$node[[i]]
    args$start <- object$info$coefficients[[i]]
    if (depth(object$node) < 1L) {
      args$formula <- object$info$formula$root
    } else {
      args$data$Part <-
        tvcolmm_get_part(object, object$data, object$fitted[, "(weights)"])
      args$contrasts <- appendDefArgs(list(Part = contrasts(args$data$Part)),
                                      args$contrasts)
      args$formula <- object$info$formula$tree
    }
    object$info$model <- suppressWarnings(do.call("olmm", args))
    
    predi <- predict(object, newdata = newdata, type = tp, ...)
    if (oob) predi[folds[,i] > 0] <- NA
    
    if (i == 1) {
      pred <- count <- predi
      pred[] <- count[] <- 0
    }
    
    subs <- !is.na(predi)
    pred[subs] <- pred[subs] + predi[subs]
    count <- count + 1 * subs
    
  }
  
  if (verbose) cat(" OK\n")
  
  pred <- pred / count
  pred[count == 0] <- NA

  n <- nrow(pred)

  if (type %in% c("terms", "link")) {
    rval <- pred

  } else {
    
    if (!is.logical(ranef)) stop("'ranef' must be logical.")

    yname <- all.vars(object$info$formula$root)[1]
    ylevs <- levels(object$fitted[,"(response)"])
    
    if (is.logical(ranef) && !ranef) {

      ## marginal prediction
      
      probs <- matrix(0, nrow(pred), ncol(pred) + 1)
      colnames(probs) <- ylevs
      rownames(probs) <- rownames(pred)
      
      ## prepare formula for model without fixed effects
      form <- as.Formula(object$info$formula$root)
      reEtaVar <- findbars(formula(form, rhs = 2))[[1]]
      reEtaVar <- ifelse(is.null(reEtaVar), "1",
                         paste("1 + (", deparse(reEtaVar), ")", sep = ""))
      reEtaInv <- findbars(formula(form, rhs = 1))[[1]]
      reEtaInv <- ifelse(is.null(reEtaInv), "",
                         paste("(", deparse(reEtaInv), ")", sep = ""))
      form <-
        formula(paste(yname, " ~ ", reEtaInv, " | ", reEtaVar, sep = ""))
      
      ## create an 'olmm' object
      newdata[, yname] <- # add an artificial response
        factor(sample(ylevs, nrow(newdata), replace = TRUE),
               levels = ylevs, ordered = TRUE)
      args <- appendDefArgs(list(formula = form, data = newdata,
                                 doFit = FALSE), object$info$dotargs)
      model <- do.call("olmm", args)
      
      ## run the marginal prediction with 'pred' as linear predictor
      .Call("olmm_pred_marg", model, pred, model@W, nrow(pred), probs,
            PACKAGE = "vcolmm")
      
    } else {

      ## conditional prediction
      
      ## evaluate probabilities at ranef = 0
      dims <- extract(object, "model")@dims
      
      if (dims["family"] == 1L) {    
        
        linkFUN <- switch(dims["link"], plogis, pnorm, pcauchy)
        
        cumProbs <- linkFUN(pred)
        cumProbs <- cbind(cumProbs, rep(1, n))
        probs <- cbind(cumProbs[, 1], t(apply(cumProbs, 1, diff)))
        
      } else if (dims["family"] %in% c(2L, 3L)) {
        
        probs <- exp(pred) / (1 + matrix(rowSums(exp(pred)),
                                         n, dims["nEta"]))
        probs <- cbind(probs, 1 - rowSums(probs))
      }
      colnames(probs) <- levels(object$fitted[, "(response)"])
      rownames(probs) <- rownames(pred)
    }
   
    if (type == "class") {
      
      rval <- apply(probs, 1, which.max)
      rval <- factor(ylevs[rval], levels = ylevs, ordered = TRUE)
      names(rval) <- rownames(newdata)
    } else {
      rval <- probs
    }
  }
  
  return(rval)
}

print.fvcolmm <- function(x, ...) {

  cat("\tVarying coefficient olmm forest\n")
  cat("\nNumber of trees:", length(x$info$node), "\n")
  return(invisible(x))
}
