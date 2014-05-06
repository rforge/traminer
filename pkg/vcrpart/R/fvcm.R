## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-05-02
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
## print.fvcm:   print method for 'forest.tvcm' objects
##
## To do:
## - print.fvcm: 'data' output
## --------------------------------------------------------- #


fvcm <- function(..., folds = cvfolds(object, "subsampling", 50L),
                    control = fvcm_control()) {
  
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

fvcm_control <- function(alpha = 1.0, maxwidth = 10L, mtry = 5L, ...) {
  return(tvcm_control(alpha = alpha, maxwidth = maxwidth, mtry = mtry, ...))
}

fitted.fvcm <- function(object, ...) predict(object, ...)

oobrisk.fvcm <- function(object, fun = NULL, ...) {

  ## estimate error for an observation with tree in which the observation
  ## wasn't included
  
  if (is.null(fun)) {
    fun <- function(y, mu, wt)
      sum(object$info$family$dev.resids(y, mu, wt))
  }
  weights <- weights(object)
  yMat <- model.matrix(~ -1 + object$fitted[,"(response)"])
  if (object$info$family$family == "binomial" && nrow(yMat) > 1L)
    yMat <- yMat[,2L,drop = FALSE]
  mu <- predict(object, type = "response", oob = TRUE, ...)
  rval <- fun(yMat, mu, weights)
  return(rval)
}

predict.fvcm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class", "coef"),
                         ranef = FALSE, 
                         verbose = FALSE, ...) {
  
  type <- match.arg(type)
  tp <- if (type == "coef") type else "link"
  family <- object$info$family
  
  if (!"oob" %in% list(...)) oob <- FALSE
  if (!"which" %in% list(...)) which <- 1:length(object$info$node)
  
  ## modifiy class to apply 'tvcm' methods
  class(object) <- class(object)[-1L]
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
  args <- object$info$dotargs
  args$data <- model.frame(object$info$model)
  args$family <- family
  if (object$info$fit == "olmm") args$doFit <- FALSE
  if (object$info$fit == "glm") 
    args$control <- glm.control(epsilon = 10.0, maxit = 1)

  ranefMat <- ranef
  
  if (verbose) cat("* predicting ")
  
  ## extract each tree and aggregate predicted coefficients
  for (i in 1:K) {

    if (verbose) cat(".")
    
    ## build model
    object$node <- object$info$node[[i]]
    args$start <- object$info$coefficients[[i]]
    if (depth(object$node) < 1L) {
      args$formula <- object$info$formula$root
    } else {
      args$formula <- object$info$formula$tree
      args$data$Node <-
        tvcm_get_node(object, object$data, TRUE, object$fitted[, "(weights)"])
      args$contrasts$Node <- contrasts(args$data$Node)
    }
    object$info$model <- suppressWarnings(do.call(object$info$fit, args))

    ## compute random effects if 'ranef = TRUE'
    if (is.logical(ranef) && ranef && inherits(family, "family.olmm")) {
      .Call("olmm_update_u", object$info$model, PACKAGE = "vcrpart")
      ranefMat <- ranef(object$info$model)
    }

    ## overwrite coefficients
    if (!inherits(family, "family.olmm"))
      object$info$model$coefficients <- object$info$coefficients[[i]]

    ## predict 'eta'
    predi <- predict(object, newdata = newdata, type = tp, ranef = ranefMat, ...)
    
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
  if (!is.matrix(pred)) pred <- matrix(pred, ncol = 1L)
  
  n <- nrow(pred)

  if (type %in% c("coef", "link")) {
    
    rval <- pred

  } else {
    
    yName <- all.vars(object$info$formula$root)[1]
    yLevs <- if (inherits(family, "family.olmm")) {
      levels(args$data[, yName])
    } else {
      yName
    }
          
    rval <- matrix(0, nrow(pred), length(yLevs),
                   dimnames = list(rownames(pred), yLevs))
    colnames(rval) <- yLevs
    rownames(rval) <- rownames(pred)

    ## create an 'empty' formula
    start <- NULL
    if (inherits(family, "family.olmm")) {

      terms <- "fe(intercept=FALSE)"
      
      ## add random effect terms to 'form' for marginal prediction
      mTerms <- terms(object$info$formula$root, specials = "re")
      if (!ranef && length(subs <- attr(mTerms, "specials")$re) > 0L) {
        terms <- c(terms, rownames(attr(mTerms, "factors"))[subs])
        reTerms <-
          grep("ranefCholFac", names(object$info$coefficients[[1L]]), value = TRUE)
        start <- sapply(1:K, function(i) object$info$coefficients[[i]][reTerms])
        start <- apply(matrix(start, ncol = length(reTerms)), 2L, mean)
        names(start) <- reTerms
      } else {
        start <- NULL
      }
    } else {
      terms <- "-1"
    }
    form <- as.formula(paste(yName, "~", paste(terms, collapse = "+")))
    
    ## create an 'olmm' object
    newdata[, yName] <- rep(args$data[, yName], length.out = nrow(newdata))
    args <- appendDefArgs(list(formula = form,
                               data = newdata,
                               offset = pred,
                               start = start),
                          object$info$dotargs)
    if (inherits(family, "family.olmm")) args$doFit <- FALSE    
    model <- do.call(object$info$fit, args)

    ## predict
    rval <- predict(model, type = type, ranef = ranefMat, ...)
  }
  
  return(rval)
}

print.fvcm <- function(x, ...) {
  cat("Random forest-based varying-coefficients model\n\n")
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
              ", number of trees = ", length(x$info$node),
              ")\n", sep = ""))
  if (nzchar(mess <- naprint(attr(x$data, "na.action")))) 
    cat("\n(", mess, ")\n", sep = "")
  return(invisible(x))
}
