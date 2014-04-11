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

tvcolmm_get_estimates <- function(object, what = c("coef", "sd", "var"), ...) {

  what <- match.arg(what)
  
  ids <- nodeids(object, terminal = TRUE)
  control <- extract(object, "control")

  rval <- list()
  
  ## extract coefficients
  coef <- switch(what,
                 coef = coef(extract(object, "model")),
                 sd = diag(vcov(extract(object, "model"))),
                 var = diag(vcov(extract(object, "model"))))

  ## restricted coefficients
  rval$restricted <- coef[!grepl("Part", names(coef))]
  
  ## varying coefficients
  if (depth(object) > 0L) {
    
    ## create object with all possible coefficients for each node
    terms <- tvcolmm_get_terms(object)
    varcoef <- rep(NA, length(terms) * length(ids))
    for (i in 1:length(terms))
      for (j in 1:length(ids))
        names(varcoef)[length(ids) * (i - 1) + j] <-
          sub("Part", paste("Part", ids[j], sep = ""), terms[i])
    subs <- intersect(names(varcoef), names(coef))
    varcoef[subs] <- coef[subs]
    
    if ((subs <- paste("Part", max(ids), sep = "")) %in%
        names(varcoef)) {
      con <- extract(object, "model")@contrasts$Part

      if (what == "coef") {
        varcoef[subs] <-
          sum(con[as.character(max(ids)),] *
              coef[paste("Part", setdiff(ids, max(ids)), sep = "")])
        
      } else if (what %in% c("sd", "var")) {
        varcoef[subs] <-
          sum((con[as.character(max(ids)),])^2 *
              coef[paste("Part", setdiff(ids, max(ids)), sep = "")])
      }
    }
    
    ## create a matrix of coefficients
    FUN <- function(i) {
      name <- paste("Part", i, sep = "")
      parts <- strsplit(names(varcoef), ":")
      subs <- sapply(parts, function(x) sum(x == name) > 0)
      rval <- varcoef[subs]
      names(rval) <- sub(name, "Part", names(rval))
      names(rval) <- sub("Part:", "", names(rval))
      return(rval)
    }
    
    rval$varying <- lapply(as.character(ids), FUN)
    rval$varying <-
      matrix(unlist(rval$varying), nrow = length(ids), byrow = TRUE,
             dimnames = list(ids, names(rval$varying[[1]])))
    
  } else {
    rval$varying <- NULL
  }

  if (what == "sd") {
    rval$restricted <- sqrt(rval$restricted)
    rval$varying <- sqrt(rval$varying)
  }
  
  return(rval)
}

coef.tvcolmm <- function(object, ...) tvcolmm_get_estimates(object, ...)

coefficients.tvcolmm <- coef.tvcolmm 

extract.tvcolmm <- function(object, what = c("control", "sctest",
                                      "model", "selected", "p.value",
                                      "coef", "sd", "var"),
                            ids = nodeids(object), ...) {
  
  what <- match.arg(what)

  if (what == "control") {

    return(object$info$control)
    
  } else if (what == "sctest") {

    if (depth(object) == 0L) return(list("1" = object$info$sctest))
    
    getSctest <- function(node) {
      if (is.terminal(node)) {
        rval <- info_node(node)$sctest
      } else {
        rval <- split_node(node)$info$sctest
      }
      return(rval)
    }
    
    return(nodeapply(object, ids, getSctest))
    
    
  } else if (what == "model") {

    return(object$info$model)
    
  } else if (what == "selected") {

    ids <- setdiff(nodeids(object,, FALSE), nodeids(object,, TRUE))
    getSelected <- function(node) node$split$varid
    rval <- unique(unlist(nodeapply(object, ids, getSelected)))
    rval <- colnames(object$data)[rval]
    return(rval)

  } else if (what == "p.value") {

    getPval <- function(node) {
      if (is.terminal(node)) {
        test <- info_node(node)$sctest
      } else {
        test <- split_node(node)$info$sctest
      }
      return(min(test$p.value, na.rm = TRUE))
    }
    rval <- unlist(nodeapply(object, ids, getPval))
    return(rval)
  } else if (what %in% c("coef", "sd", "var")){

    return(tvcolmm_get_estimates(object, what = what))

  }
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
    weights <- matrix(weights(extract(object, "model")), 
                      nrow(folds), ncol(folds))
    rval <- -sum(cv$loss[, 1] * colSums(folds == 0 * weights))
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

prune.tvcolmm <- function(tree, alpha = NULL, depth = NULL,
                          minsplit = NULL, nselect = NULL, ...) {
  
  call <- getCall(tree$info$model)
  call$formula <- tree$info$formula$tree
  call$control <- tree$info$control
  call$data <- model.frame(tree)

  ## prune the tree structure
  node <- tvcolmm_prune_node(tree, alpha, depth, minsplit, nselect)

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
    tree$info$test <-
      tvcolmm_fit_sctest(tree$info$model, tree$nodes, tree$data, control)
  }
  if (!is.null(alpha)) 
    tree$info$control$alpha <- alpha
  if (!is.null(depth))
    tree$info$control$depth <- depth
  if (!is.null(minsplit))
    tree$info$control$minsplit <- minsplit
  return(tree)
}


print.tvcolmm <- function(x, ...) {

  so <- summary(x$info$model)
  etaVar <- so$FEmatEtaVar
  etaInv <- so$FEmatEtaInv
  
  header_panel <- function(x) {
    rval <- paste(x$info$title, "\n", sep = "")
    if (!is.null(x$info$call$formula))
      rval <- paste(rval, " Family: ", x$family, "\n", sep = "")
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
                  "\n\t alpha = ", format(x$info$control$alpha, ...),
                  if (x$info$control$bonferroni) " (Bonferroni corrected)",
                  "\n\t minsplit = ", x$info$control$minsplit, "\n", sep = "")

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
  
  partykit::print.party(x, header_panel = header_panel,
                        terminal_panel = terminal_panel, ...)
  
}

ranef.tvcolmm <- function(object, ...)
  return(ranef(object$info$model, ...))

resid.tvcolmm <- function(object, ...)
  return(resid(object = object$info$model, ...))

residuals.tvcolmm <- resid.tvcolmm

weights.tvcolmm <- function(object, ...) {
  weights(extract(object, "model"))
}
