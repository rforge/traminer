## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-09
##
## Description:
## Random forests implementation based on the 'tvcm' algorithm.
##
## Contents:
## fvcm:         random forest algorithm
## fvcm_control: control function for 'fvcm'
## fitted.fvcm:  extracts fitted values
## oobrisk.fvcm: extracts out-of-bag risk
## predict.fvcm: prediction for 'forest.tvcm' objects
## plot.fvcm:    plot method for 'forest.tvcm' objects
## print.fvcm:   print method for 'forest.tvcm' objects
##
## To do:
## - print.fvcm: 'data' output
## --------------------------------------------------------- #

fvcolmm <- function(..., family = cumulative(), folds,
                    control = fvcm_control()) {
  mc <- match.call()
  mc[[1L]] <- as.name("fvcm")
  mc$fit <- "olmm"
  if ("weights" %in% names(mc)) mc$weights <- list(...)$weights
  return(eval.parent(mc))
}

fvcglm <- function(..., folds, control = fvcm_control()) {
  mc <- match.call()
  mc[[1L]] <- as.name("fvcm")
  mc$fit <- "glm"
  if ("weights" %in% names(mc)) mc$weights <- list(...)$weights
  return(eval.parent(mc))
}

fvcm <- function(..., folds, control = fvcm_control()) {
  
  mc <- match.call()
  args <- list(...)

  ## get fitting function (avoids error in tvcm)
  if (!is.null(args$fit) && is.function(args$fit))
    args$fit <- deparse(mc$fit)
  
  ## modify control parameters temporarily for a first tree
  maxwidth <- control$maxwidth
  control$maxwidth <- 1L
  verbose <- control$verbose
  control$verbose <- FALSE
  
  ## fit a prototyp tree
  if (verbose) cat("* fitting an initial tree ... ")
  initargs <- append(args, list(control = control))
  object <- do.call("tvcm", args = initargs)
  if (verbose) cat("OK\n")
  
  ## reset the depth parameter and set the verbose parameter
  object$info$control$maxwidth <- maxwidth
  folds <-
    folds[rownames(args$data) %in% rownames(model.frame(object)),, drop = FALSE]
  
  ## compute trees for subsamples
  args$object <- object
  args$folds <- folds
  args$type <- "forest"
  args$fixed <- TRUE
  args$verbose <- verbose
  
  cv <- do.call("cvrisk", args = args)
  fails <- cv$error$which  
  
  ## add new information to info slot
  if (length(fails) > 0)
    folds <- folds[, -fails, drop = FALSE]
  
  object$info$node <- append(object$info$node, cv$node)
  object$info$coefficients <-
    append(cv$coefficients, object$info$coefficient)
  object$info$folds <- cbind(object$info$folds, folds)
  object$info$error <- cv$error
  object$info$control$verbose <- verbose
  object$info$call <- mc
  
  ## modifiy class attribute to allow methods
  class(object) <- append("fvcm", class(object))
  return(object)
}

fvcm_control <- function(alpha = 1.0, maxwidth = 10L,
                         minbucket = 50L, mtry = 5L, ...)
  return(tvcm_control(alpha = alpha, maxwidth = maxwidth,
                      minbucket = minbucket, mtry = mtry, ...))

fitted.fvcm <- function(object, ...) predict(object, ...)

oobrisk.fvcm <- function(object, fun = NULL, ranef = FALSE, ...) {

  ## estimate error for an observation with tree in which the observation
  ## wasn't included
  
  if (is.null(fun)) {
    fun <- function(y, mu, wt)
      sum(object$info$family$dev.resids(y, mu, wt), na.rm = TRUE)
  }
  weights <- weights(object)
  yMat <- model.matrix(~ -1 + object$fitted[,"(response)"])
  if (object$info$family$family == "binomial" && nrow(yMat) > 1L)
    yMat <- yMat[,2L,drop = FALSE]
  mu <- predict(object, type = "response", ranef = ranef,
                na.action = na.pass, oob = TRUE)  
  if (any(is.na(mu))) {
    warning("some observations could not be predicted out of bag. ",
            "The oob error is reweighted.")

    if (is.matrix(mu)) {
      subs <- apply(mu, 1L, function(x) !all(is.na(x)))
      mu <- mu[subs,,drop=FALSE]
    } else {
      subs <- !is.na(mu)
      mu <- mu[subs]
    }
    sumOfWeights <- sum(weights)
    weights <- weights[subs]
    weights <- weights / sum(weights) * sumOfWeights
    yMat <- yMat[subs,,drop=FALSE]
  }
  return(fun(yMat, mu, weights))
}


plot.fvcm <- function(x, type = c("default", "coef", 
                           "simple", "terms"),
                      which = 1L, ask = TRUE, ...) {

  type <- match.arg(type)
  dotargs <- list(...)

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  ## set call
  call <- call(name = "plot.tvcm", x = quote(x), type = type)  
  for (arg in names(dotargs)) call[[arg]] <- dotargs[[arg]]
    
  ## modify model object if coefficients plots are called

  if (type == "terms") {

    call$ask <- ask
    eval(call)

  } else {

    if (type == "coef" | (type == "default" & depth(x) > 3L)) {
      
      ## modify default arguments
      if (!is.null(dotargs$conf.int) && dotargs$conf.int)
        warning("'conf.int' is not available for 'fvcm' objects")
      call$conf.int <- FALSE
      if (!is.null(dotargs$mean) && dotargs$mean)
        warning("'mean' is not available for 'fvcm' objects")
      call$mean <- FALSE
    }
    
    which <- intersect(which, 1:length(x$info$node))

    ## loop through all trees
    for (i in which) {

      ## set nodes
      x$node <- x$info$node[[i]]

      ## set coefficients (this is a ugly hack)
      if (type == "coef" | (type == "default" & depth(x) > 3L)) {
        if (isS4(x$info$model)) {
          slot(x$info$model, "coefficients") <- x$info$coefficients[[i]]
        } else {
          x$info$model$coefficients <- x$info$coefficients[[i]]
        }
      }
      eval(call)
    }
  }
}

predict.fvcm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class", "coef",
                           "ranef"),
                         ranef = FALSE, na.action = na.pass,
                         verbose = FALSE, ...) {
  
  type <- match.arg(type)

  ## check newdata
  if (!is.null(newdata) && !class(newdata) == "data.frame")
    stop("'newdata' must be a 'data.frame'.")
  
  ## resolve conflicts with the 'ranef' argument
  if (!class(ranef) %in% c("logical", "matrix"))
    stop("'ranef' must be a 'logical' or a 'matrix'.")
  if (!is.null(newdata) && is.logical(ranef) && ranef)
    stop("'ranef' should be 'FALSE' or a 'matrix' if 'newdata' is not 'NULL'.")
  if (type == "ranef" & (!is.logical(ranef) | is.logical(ranef) && ranef))
    stop("for 'type = 'ranef'' the argument 'ranef' must be 'FALSE'.")
  if (type == "ranef" & !is.null(newdata))
    stop("prediction for random effects for 'newdata' is not implemented.")
  
  ## get and set hidden arguments
  oob <- if (!is.null(list(...)$oob)) list(...)$oob else FALSE
  
  ## check hidden arguments
  if (oob && !is.null(newdata))
    stop("'oob' should be 'FALSE' if 'newdata' is not 'NULL'")
  
  ## modifiy class to apply 'tvcm' methods
  class(object) <- class(object)[-1L]
  
  ## get training data
  mf <- model.frame(object)
  
  ## save the original model
  dummymodel <- object$info$model
  
  ## set newdata
  if (is.null(newdata)) newdata <- mf

  ## set oob folds for option 'oob'
  folds <- if (oob) {
    object$info$folds
  } else {
    matrix(1L, nrow(newdata), length(object$info$node))
  }
  
  ## get formulas
  formList <- vcrpart_formula(object$info$formula$root)

  ## extract the name and levels of the response
  yName <- all.vars(object$info$formula$root)[1]
  yLevs <- if (object$info$fit == "olmm") levels(mf[, yName]) else yName
  nYLevs <- length(yLevs)

  if (type != "coef") {
    
    ## check and set response
    if (!yName %in% colnames(newdata))
      newdata[, yName] <- sample(mf[, yName], nrow(newdata), replace = TRUE)

    ## check fixed effects predictors
    feVars <- unlist(lapply(formList$fe$eta, all.vars))
    if (!all(subs <- feVars %in% colnames(newdata)))
      stop("variable(s) ", paste("'", feVars[!subs], "'", collapse = ", "),
           " are not available in 'newdata'.")

    ## check and set random effect predictors
    reVars <- unlist(lapply(unlist(formList$re[c("eta", "cond")]), all.vars))
    if (length(reVars) > 0L) {
      if (is.logical(ranef) && ranef | is.matrix(ranef)) {
        if (!all(subs <- reVars %in% colnames(newdata)))
          stop("variables ", feVars[!subs], " are not available in 'newdata'.")
      } else {
        for (var in reVars)
          newdata[, var] <- sample(mf[, var], nrow(newdata), replace = TRUE)
      }
    }

    ## check and set newdata
    Terms <- attr(mf, "terms")
    xlevels <- .getXlevels(attr(mf, "terms"), mf)
    if (is.matrix(ranef)) {
      subjectName <- slot(dummymodel, "subjectName")
      xlevels <- xlevels[names(xlevels) != subjectName]
    }
    
    newdata <- as.data.frame(model.frame(Terms, newdata,
                                         na.action = na.pass,
                                         xlev = xlevels))
  }
  
  if (verbose) cat("* predicting the coefficient function ... ")

  ## ------------------------------------------------------- #
  ## Step 1: predict the coefficients for each observation
  ## ------------------------------------------------------- #

  nEta <- if (object$info$fit == "olmm") nYLevs - 1L else 1L
  etaLabs <- paste("Eta", 1L:nEta, sep = "")
  coef <- count <- 0 * predict(object, newdata, type = "coef")
  rownames(coef) <- rownames(newdata)
  subs <- matrix(TRUE, nrow(coef), ncol(coef))
  
  for (i in seq_along(object$info$node)) {

    if (verbose) cat(".")

    ## set node
    object$node <- object$info$node[[i]]

    ## set coefficients
    if (object$info$fit == "olmm") {     
        slot(object$info$model, "coefficients") <- object$info$coefficients[[i]]
    } else {
        object$info$model[["coefficients"]] <- object$info$coefficients[[i]]
    }

    coefi <- predict(object, newdata = newdata, type = "coef",
                     ranef = FALSE, na.action = na.pass, ...)
    if (!is.matrix(coefi)) coefi <- matrix(coefi, nrow = nrow(newdata))

    ## index matrix for valid entries
    subsi <- subs
    
    ## remove observations that appear in the learning sample
    if (oob) subsi[folds[,i] > 0L, ] <- FALSE

    ## acount for skipped categories
    if (object$info$fit == "olmm" && ncol(coefi) < ncol(coef)) {        
        subsiCols <- table(mf[folds[, i] > 0, yName]) > 0L
        subsiCols <- subsiCols[-length(subsiCols)]
        etaLabsShould <- etaLabs[subsiCols]
        colnamesi <- colnames(coefi)
        etaLabsIs <- grep("Eta[1-9]+:", colnamesi, value = TRUE)
        etaLabsIs <- unique(sapply(strsplit(etaLabsIs, ":"), function(x) x[1]))
        colnamesi <- strsplit(colnamesi, ":")
        for (j in rev(seq_along(etaLabsIs))) {
            colnamesi <- lapply(colnamesi, function(x) {
                if (x[1L] == etaLabsIs[j]) x[1L] <- etaLabsShould[j]
                return(x)
            })
        }
        colnamesi <- sapply(colnamesi, function(x) paste(x, collapse = ":"))
        colnames(coefi) <- colnamesi          
    } else {
        subsiCols <- rep(TRUE, ncol(coefi))
    }

    subsi <- subs
    subsi[, !colnames(coef) %in% colnames(coefi)] <- FALSE
    if (oob) subsi[folds[,i] > 0L,] <- FALSE

    subsiC <- intersect(colnames(coef), colnames(coefi))
    subsiR <- if (oob) folds[,i] == 0L else subs[,1]
    coef[subsi] <- coef[subsi] + c(coefi[subsiR, subsiC])
    count <- count + 1 * subsi
  }
  
  if (verbose) cat(" OK\n")
  
  coef <- coef / count
  coef[apply(count, 1, function(x) any(x == 0)), ] <- NA
  coef <- coef[, names(coef(dummymodel))]
  
  ## ------------------------------------------------------- #
  ## Step 2: predict the linear predictor for each observation
  ## ------------------------------------------------------- #

  if (type == "coef") return(na.action(coef))
    
  ## create a model matrix 'X'
  if (object$info$fit == "olmm") {
    X <- olmm_merge_mm(model.matrix(terms(formList$fe$eta$ce, keep.order = TRUE),
                                    newdata, attr(slot(object$info$model, "X"),
                                                  "contrasts")),
                       model.matrix(terms(formList$fe$eta$ge, keep.order = TRUE),
                                    newdata, attr(slot(object$info$model, "X"),
                                                  "contrasts")), TRUE)
  } else {
    X <- model.matrix(terms(object$info$formula$root), newdata,
                      object$info$model$contrasts)
  }
  
  ## compute the linear predictor 'eta' based on 'coef' and 'X'
  if (object$info$fit == "olmm") {
    coef <- coef[, substr(colnames(coef), 1,12) != "ranefCholFac", drop = FALSE]
    dims <- slot(dummymodel, "dims")
    fixefMat <- function(fixef) {
      return(rbind(matrix(fixef[1:(dims["pCe"] * dims["nEta"])], dims["pCe"], dims["nEta"], byrow = FALSE), if (dims["pGe"] > 0) matrix(rep(fixef[(dims["pCe"] * dims["nEta"] + 1):dims["p"]], each = dims["nEta"]), dims["pGe"], dims["nEta"], byrow = TRUE) else NULL))
    }
    eta <- t(sapply(1:nrow(newdata), function(i) {
      X[i,,drop = FALSE] %*% fixefMat(coef[i,])
    }))
  } else {     
    eta <- t(sapply(1:nrow(newdata), function(i) {
      X[i,,drop = FALSE] %*% coef[i, ]
    }))
    eta <- matrix(eta, ncol = 1L)
  }
  colnames(eta) <- etaLabs
  rownames(eta) <- rownames(newdata)

  if (type == "link") return(na.action(eta))

  ## ------------------------------------------------------- #
  ## Step 3: create a new 'empty' model and predict the outcomes
  ## ------------------------------------------------------- #
      
  ## set the formula 'form' and the intial values 'start'
  start <- NULL
  if (object$info$fit == "olmm") { # the formula for olmms
    terms <- "fe(intercept=FALSE)"
    
    ## add random effect terms
    mTerms <- terms(object$info$formula$root, specials = "re")
    if (length(subs <- attr(mTerms, "specials")$re) > 0L) {
      
      ## the random effect term
      terms <- c(terms, rownames(attr(mTerms, "factors"))[subs])
      reTerms <-
        grep("ranefCholFac", names(object$info$coefficients[[1L]]), value = TRUE)
      
      ## compute the mean estimated random effect variance
      start <- sapply(seq_along(object$info$coefficients),
                      function(i) object$info$coefficients[[i]][reTerms])
      start <- apply(matrix(start, ncol = length(reTerms)), 2L, mean)
      names(start) <- reTerms  
    }
  } else { # the formula for glms
    terms <- "-1" 
  }
  form <- as.formula(paste(yName, "~", paste(terms, collapse = "+")))
  
  ## ensure for olmms that each response category is available
  if (is.factor(newdata[, yName]) && length(unique(newdata[, yName])) < nYLevs) {
    subs <- nrow(newdata) + 1L:nYLevs
    newdata <- rbind(newdata, newdata[rep(1L, nYLevs),,drop = FALSE])
    newdata[subs, yName] <- yLevs
    if (object$info$fit == "olmm") {
      sN <- slot(object$info$model, "subjectName")
      levs <- c(levels(newdata[,sN]), "RetoBuergin") 
      newdata[sN] <- factor(newdata[,sN], levels = levs)
      newdata[subs, sN] <- "RetoBuergin"
    }
    eta <- rbind(eta, matrix(0, nYLevs, ncol(eta)))
    folds <- rbind(folds, matrix(-1L, nYLevs, length(object$info$node)))
  }

  ## set the offset of the model as predicted linear predictors
  offset <- eta
  offset[is.na(offset)] <- 0 # just to get the calls working
  
  ## create a call for the 'empty' model
  oobCall <- call(name = object$info$fit,
                  form = quote(form),
                  data = quote(newdata),
                  offset = quote(offset),
                  start = quote(start),
                  na.action = na.pass)
  for (arg in names(object$info$dotargs))
    oobCall[[arg]] <- object$info$dotargs[[arg]]
  oobCall <- oobCall[!duplicated(names(oobCall))]
  
  if (object$info$fit == "olmm") oobCall$doFit <- FALSE
  
  ## fit the 'empty model'
  model <- suppressWarnings(eval(oobCall))
 
  if (type == "ranef") {

    ## predict random effects
    ranef <- ranef(model) 
    ranef <- ranef[rownames(ranef) != "RetoBuergin",,drop=FALSE]
    return(na.action(ranef))
    
  } else {

    ## predict outcomes
    if (is.matrix(ranef)) {
      ranefMat <- ranef(model)
      ranefMat[rownames(ranef), ] <- ranef
      ranef <- ranefMat
    }
    pred <- predict(model, type = type, ranef = ranef, ...)
    if (!is.matrix(pred)) pred <- matrix(pred, nrow = nrow(newdata))
    pred[apply(eta, 1, function(x) any(is.na(x))),] <- NA
    pred <- pred[folds[,1] >= 0L,, drop = FALSE] # drop added observations
  }
  folds <- folds[folds[,1] >= 0L,,drop = FALSE] # drop added folds
  
  ## set the observations which appear in all trees to NA
  if (oob) pred[apply(folds, 1L, function(x) all(x == 0L)),] <- NA

  ## return predictions
  return(na.action(pred))
}


print.fvcm <- function(x, ...) {
  cat(if (x$info$control$mtry < Inf) "Random forest" else "Bagging",
      "based varying-coefficients model\n\n")
  if (length(x$info$family$family) > 0L)
    cat("  Family:", x$info$family$family, x$info$family$link, "\n")
  if (length(x$info$formula$original) > 0L)
    cat(" Formula:", paste(deparse(x$info$formula$original), collapse = "\n"), "\n")
  if (length(str <- deparseCall(x$info$call$data)) > 0L)
    cat("    Data: ", str, "\n", sep = "")
  if (length(str <- deparseCall(x$info$call$subset)) > 0L)
    cat("  Subset: ", str, "\n", sep = "")
  if (length(x$info$control) > 0L)
    cat(paste("  Method: ", x$info$control$method,
              " (maxwidth = ", x$info$control$maxwidth,
              ", minbucket = ", x$info$control$minbucket,
              ", ntrees = ", length(x$info$node),
              ")\n", sep = ""))
  if (nzchar(mess <- naprint(attr(x$data, "na.action")))) 
    cat("\n(", mess, ")\n", sep = "")
  return(invisible(x))
}

ranef.fvcm <- function(object, ...) {
  return(predict(object, type = "ranef", ...))
}
