## --------------------------------------------------------- #
## Author:          Reto Buergin
## E-Mail:          reto.buergin@unige.ch, rbuergin@gmx.ch
## Date:            2014-17-03
##
## Description:
## S3 methods for tvcolmm objects
##
## References:
## party:           http://CRAN.R-project.org/package=party
## partykit:        http://CRAN.R-project.org/package=partykit
##
## Methods:
## coef.tvcolmm:    coefficients.tvcolmm
##                  tree structure)
## getCall:         extract original call
## logLik.tvcolmm:  extract log Likelihood
## model.frame:
## neglogLik.tvcolmm: extract negative log Likelihood
## prune:
## predict.tvcolmm: predict responses (see prediction of 'olmm' class)
## print.tvcolmm:   print tvcolmm objects
## resid.tvcolmm, residuals.tvcolmm: extract residuals
##
## Todo:
## - anova.tvcolmm
## --------------------------------------------------------- #

coef.tvcolmm <- function(object, ...) tvcolmm_get_estimates(object, ...)

coefficients.tvcolmm <- coef.tvcolmm 

extract.tvcolmm <- function(object, what = c("control", "model", 
                                      "sctest", "p.value",
                                      "lossgrid", "selected", 
                                      "coef", "sd", "var"),
                            ids = nodeids(object), steps = NULL, ...) {
  
  what <- match.arg(what)
  splitpath <- object$info$splitpath
  if (is.null(steps))
    steps <- which(sapply(splitpath, function(x) !is.null(unlist(x))))
  sctest <- object$info$control$sctest
  
  if (what == "control") {

    return(object$info$control)
    
  } else if (sctest && what == "sctest") {
    
    splitpath <- object$info$splitpath
    rval <- lapply(splitpath[steps], function(x) x$sctest)
    return(rval)
    
  } else if (what == "lossgrid") {
    
    splitpath <- object$info$splitpath
    rval <- lapply(splitpath[steps], function(x) x$lossgrid)
    return(rval)
    
  } else if (what == "model") {
    
    return(object$info$model)
    
  } else if (what == "selected") {

    splitpath <- object$info$splitpath
    rval <- unique(unlist(lapply(splitpath[steps], function(x) x$varid)))
    if (length(rval) > 0L) rval <- colnames(object$data)[rval]
    return(rval)
    
  } else if (sctest && what == "p.value") {
    
    rval <- unlist(sapply(splitpath[steps],
                          function(x) {
                            if (is.null(x$sctest)) return(NA)
                            rval <- na.omit(c(x$sctest$p.value))
                            if (length(rval) == 0) return(NA)
                            return(min(rval, na.rm = TRUE))
                          }))
    return(rval)
    
  } else if (what %in% c("coef", "sd", "var")){

    return(tvcolmm_get_estimates(object, what = what))
    
  }
  return(NULL)
}

fitted.tvcolmm <- function(object, ...) {
  args <- append(list(object = object), list(...))
  args$newdata <- NULL # delete the newdata argument 
  return(do.call(predict, args = args)) # ... and call predict
}

getCall.tvcolmm <- function(x, ...) return(x$info$call)

logLik.tvcolmm <- function(object, cv = FALSE,
                           folds = tvcolmm_folds(object, "kfold"),
                           ...) {
  if (cv) {
    cv <- cv.tvcolmm(object, folds, alpha.max = object$info$control$alpha, 
                     fixed = TRUE, ...)
    weights <- weights(object)
    oobWeights <- matrix(weights, nrow(folds), ncol(folds))
    rval <- -sum(cv$loss[, 1] * colSums((folds == 0) * oobWeights))
    if (sum((folds == 0) * oobWeights) !=  sum(weights))
      rval <- rval * sum(weights) / sum((folds == 0) * oobWeights)
    dims <- extract(object, "model")@dims
    attr(rval, "nall") <- attr(rval, "nobs") <- dims[["n"]]
    attr(rval, "df") <- dims[["nPar"]]
    class(rval) <- "logLik"
  }
  else {
    rval <- logLik(object$info$model)
  }
  return(rval)
}

model.frame.tvcolmm <- function(formula, ...)
  return(cbind(model.frame(formula$info$model), formula$data))

neglogLik.tvcolmm <- function(object, ...)
    return(-as.numeric(logLik(object, ...)))

predict.tvcolmm <- function(object, newdata = NULL,
                           type = c("prob", "class", "node",
                             "link", "terms"), ...) {

  ## match type
  type <- match.arg(type)

  ## the terminal node identifiers
  ids <- nodeids(object, terminal = TRUE)
  
  ## extract node ids
  if (is.null(newdata)) {
    fitted <- object$fitted[["(fitted)"]]
    names(fitted) <- rownames(object$info$model@frame)
  } else {
    fitted <- fitted_node(object$node, newdata[, colnames(object$data), drop = FALSE])
    newdata$Part <- tvcolmm_get_part(object, newdata)
  }

  ## return fitted node ids
  if (type == "node") {

    ## return fitted nodes
    return(fitted)
  } else if (type == "terms") {
      
    ## extract individual effects
    what <- list(...)$what
    if (is.null(what)) what <- "coef"
    coef <- extract(object, what)
    varcoef <- lapply(fitted, function(x) coef$varying[as.character(x), ])
    varcoef <- matrix(unlist(varcoef), nrow = length(fitted), byrow = TRUE,
                      dimnames = list(names(fitted),
                        colnames(coef$varying)))

    rescoef <- matrix(coef$restricted, nrow = length(fitted),
                      ncol = length(coef$restricted), byrow = TRUE,
                      dimnames = list(names(fitted),
                        names(coef$restricted)))

    rval <- cbind(varcoef, rescoef)
    
    cn <- names(coef(extract(object, "model")))
    cn <- sub("Part[0-9]+", "Part", cn)
    cn <- sub("Part:", "", cn)
    cn <- unique(cn)

    rval <- rval[, cn]
    
    return(rval)
  } else {

    ## call predict.olmm
    return(predict(object$info$model, newdata = newdata, type = type, ...))
  }
  return(fitted)
}

print.tvcolmm <- function(x, ...) {
  
  so <- summary(x$info$model)
  etaVar <- so$FEmatEtaVar
  etaInv <- so$FEmatEtaInv
  
  header_panel <- function(x) {
    rval <- paste(x$info$title, "\n\n", sep = "")
    if (!is.null(so$family))
      rval <- paste(rval, " Family: ", so$family, "\n", sep = "")
    if (!is.null(x$info$call$formula))
      rval <- paste(rval, "Formula: ",
                    deparse(x$info$call$formula)[1], "\n", sep = "")
    if (!is.null(x$info$call$data))
      rval <- paste(rval, "   Data: ",
                    deparse(x$info$call$data)[[1]], "\n", sep = "")
    if (!is.null(x$info$call$subset))
      rval <- paste(rval, " Subset: ",
                    deparse(x$info$call$subset),"\n", sep = "")
    rval <- paste(rval, "\nTuning-parameters:",
                  "\nalpha = ", format(x$info$control$alpha, ...),
                  if (x$info$control$bonferroni) " (Bonferroni corrected)",
                  "\nminsplit = ", x$info$control$minsplit, "\n", sep = "")

    if (any(length(so$REmat) > 0 && nrow(so$REmat) > 0 |
            sum(!grepl("Part", rownames(etaInv))) > 0 |
            sum(!grepl("Part", rownames(etaVar))) > 0))
      rval <- paste(rval, "\nRestricted effects:\n", sep = "")
    
    if (length(so$REmat) > 0  && nrow(so$REmat) > 0) {
        rval <- paste(rval, "\nRandom effects:\n", sep = "")
        rval <- paste(rval, formatMatrix(so$REmat[, 1:2, drop = FALSE], ...), sep = "")
      }
    
    if (length(etaInv) > 0  && nrow(etaInv) > 0) {
      subs <- which(!grepl("Part", rownames(etaInv)))
      if (length(subs) > 0) {
        rval <-
          paste(rval, "\nPredictor-invariant fixed effects:\n", sep = "")
        rval <- paste(rval, formatMatrix(etaInv[subs, , drop = FALSE], ...), sep = "")
      }
    }
    
    if (length(etaVar) > 0  && nrow(etaVar) > 0) {
      subs <- which(!grepl("Part", rownames(etaVar)))
      if (length(subs) > 0) {
        rval <-
          paste(rval, "\nPredictor-variable fixed effects:\n", sep = "")
        rval <- paste(rval, formatMatrix(etaVar[subs, , drop = FALSE], ...), sep = "")
      }
    }
    
    rval <- paste(rval, "\nVarying effects:\n\n", sep = "")
    return(rval)
  }
  
  terminal_panel <- function(node) {

    rval <- "\n"
    
    if (length(etaInv) > 0  && nrow(etaInv) > 0) {
      subs <- rownames(etaInv) == paste("Part", node$id, sep = "") | grepl(paste("Part", node$id, ":", sep = ""), rownames(etaInv))
      
      if (sum(subs) > 0) {
        rval <-
          paste(rval, "\nPredictor-invariant fixed effects:\n", sep = "")
        rval <- paste(rval, formatMatrix(etaInv[subs, , drop = FALSE], ...), sep = "")
      }
    }
    
    if (length(so$FEmatEtaVar) > 0  && nrow(so$FEmatEtaVar) > 0) {

      subs <- rownames(etaVar) %in% paste("Eta", 1:x$info$model@dims["J"], ":Part", node$id, sep = "") | grepl(paste("Part", node$id, ":", sep = ""), rownames(etaVar))
      
      if (sum(subs) > 0) {
        rval <-
          paste(rval, "\nPredictor-variable fixed effects:\n", sep = "")
        rval <- paste(rval, formatMatrix(etaVar[subs, , drop = FALSE], ...), sep = "")
      }
    }
    
    return(rval)
  }
  
  print.party(x, header_panel = header_panel,
              terminal_panel = terminal_panel, ...)
}

prune.tvcolmm <- function(tree, alpha = NULL,
                          depth = NULL, width = NULL,
                          minsplit = NULL, minbucket = NULL, 
                          nselect = NULL, step = NULL, ...) {
  
  call <- getCall(tree$info$model)
  call$formula <- tree$info$formula$tree
  call$control <- tree$info$control
  call$data <- model.frame(tree)

  ## prune the tree structure
  node <- tvcolmm_prune_node(tree, alpha, depth, width,
                             minsplit, minbucket,
                             nselect, step)

  ## if something changes ...
  if (!identical(node, tree$node)) {
    tree$node <- node
    if (depth(tree$node) > 0L) {
      tree$fitted[, "(fitted)"] <-
        fitted_node(tree$node, data = tree$data)
      call$data$Part <-
        tvcolmm_get_part(tree, tree$data, tree$fitted[,"(weights)"])
      call$contrasts <- appendDefArgs(list(Part = contrasts(call$data$Part)),
                                      call$contrasts)
    } else {
      call$formula <- tree$info$formula$root
    }    
    tree$info$model <- try(eval(call), silent = TRUE)
    if (inherits(tree$info$model, "try-error"))
      stop("tree model fitting failed")
    tree$info$control$terms
    control <- tree$info$control
    control$info$terms <-
      setdiff(names(fixef(tree$info$model)), control$restricted)
    tree$info$control <- control
    nsteps <- width(tree)
    tree$info$nsteps <- nsteps
    splitpath <- tree$info$splitpath[1:tree$info$nsteps]
    splitpath[[nsteps]]$varid <- NULL
    splitpath[[nsteps]]$partid <- NULL
    splitpath[[nsteps]]$cutid <- NULL
    splitpath[[nsteps]]$var <- NULL
    splitpath[[nsteps]]$cutpoint <- NULL
    if (control$sctest) splitpath[[nsteps]]$lossgrid <- NULL
    tree$info$splitpath <-
      tvcolmm_modify_splitpath(splitpath, tree$node, tree$data, tree$info$control)
  }
  if (!is.null(alpha)) 
    tree$info$control$alpha <- alpha
  if (!is.null(depth))
    tree$info$control$depth <- depth
  if (!is.null(minsplit))
    tree$info$control$minsplit <- minsplit
  return(tree)
}

ranef.tvcolmm <- function(object, ...)
  return(ranef(object$info$model, ...))

resid.tvcolmm <- function(object, ...)
  return(resid(object = object$info$model, ...))

residuals.tvcolmm <- resid.tvcolmm

splitpath.tvcolmm <- function(tree, ...) tree$info$splitpath

print.splitpath.tvcolmm <- function(x, steps = NULL, ...) {

  if (!is.null(steps)) steps <- intersect(steps, 1:length(x))
  if (is.null(steps)) steps <- 1:length(x)
  
  for (step in steps) {
    if (step != steps[1]) cat("\n")
    cat("Step:", step)
    if (is.null(unlist(x[[step]]$varid))) {
      cat(" (no splitting processed)\n")
    } else {
      cat("\nSplitting variable:", x[[step]]$var)
      cat("\nPartition:", x[[step]]$partid)
      cat("\nCutpoint: ")
      cat(paste("{",paste(x[[step]]$cutpoint, collapse = "}, {"), "}\n", sep = ""))
    }
     if (!is.null(x[[step]]$sctest) | !is.null(x[[step]]$lossgrid))
       cat("\nDetails:\n")
    if (!is.null(x[[step]]$sctest)) {
      cat("\nCoefficient constancy tests (p-value):\n")
      print(x[[step]]$sctest$p.value, ...)
    }
    if (!is.null(x[[step]]$lossgrid)) {
      cat("\nLoss-minimizing grid search:\n") 
      for (j1 in 1:length(x[[step]]$lossgrid)) {
        for (j2 in 1:length(x[[step]]$lossgrid[[j1]])) {
          cat("Variable:", names(x[[step]]$lossgrid)[j1],
              "Node:", sub("Part", "", names(x[[step]]$lossgrid[[j1]])[j2]), "\n\n")
          print(x[[step]]$lossgrid[[j1]][[j2]], ...)
        }
      }
    }
  }
}

weights.tvcolmm <- function(object, ...) {
  weights(extract(object, "model"))
}
